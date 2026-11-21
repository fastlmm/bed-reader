//! GPU reading and matrix storage, re-exported by [`gpu`](super).
use core::{marker::PhantomData, mem::size_of};
use std::{
    fs::File,
    io::{Read, Seek, SeekFrom},
};

use cubecl::{prelude::*, std::tensor::TensorHandle};
use ndarray as nd;
use ndarray::ShapeBuilder;

use super::{sealed::Value, Error};
use crate::{
    open_and_check, path_ref_to_string, resolve_index, set_up_two_bits_to_value, try_div_4, Bed,
    BedError, BedVal, Hold, ReadOptions, CB_HEADER_U64,
};

// The standardizer is for the Python bindings (and this module's tests), like the CPU one.
#[cfg(any(feature = "extension-module", test))]
pub(crate) mod standardize;

/// Units per cube.
#[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
const UNITS: u32 = 256;
/// Units per cube. `CubeCL`'s CPU runtime (internal tests) runs each unit of a `sync_cube` kernel on
/// its own spinning thread, so a cube must have no more units than the machine has cores, and
/// fewer is faster (measured on 32 cores, standardizing 500 x 10,000: 256 units hung, 32 units
/// took 2.7 s, 8 took 0.14 s, 4 took 0.05 s). Small cubes also exercise the kernels' strided loops.
#[cfg(all(test, feature = "_cubecl-cpu-test"))]
const UNITS: u32 = 4;
/// Cubes along a grid's y axis are limited to this many; kernels stride over the rest.
const MAX_GRID_Y: usize = 65_535;

/// The `CubeCL` client for GPU device 0, which GPU reads use.
///
/// A [`CubeClMatrix`] carries its own client ([`CubeClMatrix::client`]); use this one for work
/// that is not tied to a matrix.
///
/// # Errors
/// [`Error::Unavailable`] if there is no GPU.
///
/// # Example
///
/// ```rust,no_run
/// use bed_reader::gpu;
///
/// let client = gpu::client()?;
/// println!("{}", client.name());
/// # Ok::<(), Box<gpu::Error>>(())
/// ```
pub fn client() -> Result<Client, Box<Error>> {
    // This crate's own tests can run the same kernels on CubeCL's CPU runtime (internal
    // feature `_cubecl-cpu-test`), for CI without a GPU. Library builds always use CUDA.
    #[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
    let result = cubecl::Device::cuda(0);
    #[cfg(all(test, feature = "_cubecl-cpu-test"))]
    let result = cubecl::Device::cpu();
    let device = result.map_err(Error::from)?;
    Ok(device.client())
}

/// Return all GPU memory that `CubeCL` holds but no longer uses to CUDA, now.
///
/// A dropped matrix's memory otherwise returns at `CubeCL`'s next allocation. This also frees
/// `CubeCL`'s pooled small temporaries, so other GPU libraries in the process can use them.
///
/// # Errors
/// [`Error::Unavailable`] if there is no GPU; [`Error::Server`] if the GPU fails.
///
/// # Example
///
/// ```rust,no_run
/// use bed_reader::{gpu, Bed, sample_bed_file};
///
/// let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
/// let cubecl_matrix = bed.read_cubecl::<f32>()?;
/// drop(cubecl_matrix);
/// gpu::memory_cleanup()?;
/// # Ok::<(), Box<gpu::Error>>(())
/// ```
pub fn memory_cleanup() -> Result<(), Box<Error>> {
    let client = client()?;
    // Finish queued work first, so buffers that it still uses are not counted as in use later.
    pollster::block_on(client.sync()).map_err(Error::from)?;
    client.memory_cleanup().map_err(Error::from)?;
    Ok(())
}

/// A value type the `CubeCL` reader supports: `i8`, `f32`, or `f64`.
///
/// This trait is sealed: downstream crates cannot add element types. The supported
/// types are controlled by bed-reader so kernel and memory-layout assumptions stay
/// consistent with the reader's dtype support.
///
/// See [`CubeClMatrix`] for a complete read example using these types.
pub trait CubeClVal: Value + BedVal + Numeric + CubeElement + bytemuck::Pod {}
impl CubeClVal for i8 {}
impl CubeClVal for f32 {}
impl CubeClVal for f64 {}

/// A genotype matrix (individuals x SNPs) in GPU memory, from a `read_cubecl` call.
///
/// It is a thin, typed wrapper over `CubeCL`'s own types: [`CubeClMatrix::tensor_handle`] and
/// [`CubeClMatrix::client`] give what is needed to run `CubeCL` kernels on it.
///
/// # Example
///
/// ```rust,no_run
/// use bed_reader::{Bed, ReadOptions, sample_bed_file};
///
/// let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
/// let cubecl_matrix = ReadOptions::builder().f32().f().read_cubecl(&mut bed)?;
/// assert_eq!(cubecl_matrix.dim(), (3, 4));
/// assert!(cubecl_matrix.is_f());
/// let tensor_handle = cubecl_matrix.tensor_handle(); // CubeCL's own type, for use with ...
/// let client = cubecl_matrix.client(); // ... its client
/// assert_eq!(tensor_handle.shape().as_slice(), &[3, 4]);
/// println!("{}", client.name());
/// // Avoid this in real work: copying to the CPU gives up the point of a GPU read. It is here
/// // for checking small results. `to_cpu` returns an `ndarray::Array2<f32>`.
/// let cpu_array2 = cubecl_matrix.to_cpu()?;
/// assert_eq!(cpu_array2.dim(), (3, 4));
/// // Select while decoding, choose another dtype and C order. No full-matrix decode is needed.
/// let selected = ReadOptions::builder().i8().c().iid_index([2, 0]).sid_index(1..4)
///     .read_cubecl(&mut bed)?;
/// assert_eq!(selected.dim(), (2, 3));
/// assert!(!selected.is_f());
/// # use bed_reader::gpu::Error;
/// # Ok::<(), Box<Error>>(())
/// ```
pub struct CubeClMatrix<TVal> {
    // Invariant: `iid_count`, `sid_count`, and `iid_count * sid_count` fit u32, as the kernels
    // index. Every constructor (`read_cubecl`, `empty`) enforces it through `element_count`, and
    // the unchecked kernel launches rely on it.
    client: Client,
    handle: cubecl::server::Handle,
    iid_count: usize,
    sid_count: usize,
    is_f: bool,
    phantom: PhantomData<TVal>,
}

/// Reads into GPU memory: the packed bytes covering the selection go to the GPU, which decodes them.
pub(crate) fn read_cubecl<TVal: CubeClVal>(
    bed: &mut Bed,
    read_options: &ReadOptions<TVal>,
) -> Result<CubeClMatrix<TVal>, Box<Error>> {
    let iid_count = bed.iid_count()?;
    let sid_count = bed.sid_count()?;
    // Indexes are `isize`, so each count must fit `isize`. A count can exceed that only when the
    // caller gives `Bed::builder` counts that no file could have.
    if isize::try_from(iid_count).is_err() || isize::try_from(sid_count).is_err() {
        Err(BedError::IndexesTooBigForFiles(iid_count, sid_count))?;
    }
    let iid_hold = Hold::new(&read_options.iid_index, iid_count)?;
    let sid_hold = Hold::new(&read_options.sid_index, sid_count)?;
    let iid_index = resolve(iid_hold.as_ref(), iid_count, BedError::IidIndexTooBig)?;
    let sid_index = resolve(sid_hold.as_ref(), sid_count, BedError::SidIndexTooBig)?;
    let element_count = element_count(iid_index.len(), sid_index.len())?;
    let client = client()?;

    let (buf_reader, header) = open_and_check(&bed.path)?;
    let mut file = buf_reader.into_inner();
    // Mode 1 stores each SNP's individuals together (SNP-major); mode 0 stores each individual's
    // SNPs together. For mode 0, swap the roles of individuals and SNPs and flip the memory
    // order: the transposed result has the requested layout (as the CPU reader does). From here
    // on, "rows" and "columns" are the file's: the output has one row per selected column.
    let (row_count, rows, column_count, columns, kernel_is_f) = match header[2] {
        1 => (
            sid_count,
            &sid_index,
            iid_count,
            &iid_index,
            read_options.is_f,
        ),
        0 => (
            iid_count,
            &iid_index,
            sid_count,
            &sid_index,
            !read_options.is_f,
        ),
        _ => Err(BedError::BadMode(path_ref_to_string(&bed.path)))?,
    };
    // The same file-length check, and error, as the CPU reader, which also swaps counts for mode 0.
    let row_bytes_u64 = try_div_4(column_count, row_count)?;
    // `try_div_4` checked that the whole file length fits u64.
    let expected_len = row_bytes_u64 * row_count as u64 + CB_HEADER_U64;
    if file.metadata()?.len() != expected_len {
        Err(BedError::IllFormed(path_ref_to_string(&bed.path)))?;
    }

    // Within each stored row, read only the bytes spanning the selected columns, packed into
    // u32 words. Selected column `k`'s 2-bit code sits in word `column_word[k]` at bit
    // `column_shift[k]` of its row's words.
    let (span_start, span_len) = match (columns.iter().min(), columns.iter().max()) {
        (Some(&min), Some(&max)) => (min / 4, max / 4 + 1 - min / 4),
        _ => (0, 0),
    };
    let span_words = span_len.div_ceil(4);
    // The kernel indexes `packed` with u32. A few far-apart columns can make it much larger than
    // the result, so check it separately from `element_count`, before allocating or reading.
    let packed_len = self::element_count(rows.len(), span_words)?;
    let column_word: Vec<u32> = columns
        .iter()
        .map(|&column| to_u32((column / 4 - span_start) / 4))
        .collect::<Result<_, _>>()?;
    let column_shift: Vec<u32> = columns
        .iter()
        .map(|&column| to_u32(((column / 4 - span_start) % 4) * 8 + (column % 4) * 2))
        .collect::<Result<_, _>>()?;
    let packed = read_packed(
        &mut file,
        rows,
        row_bytes_u64,
        span_start,
        span_len,
        span_words,
        packed_len,
    )?;

    // Empty kernel arguments advertise one dummy element; it must fit even for f64.
    let handle = empty_dedicated(&client, (element_count.max(1) * size_of::<TVal>()).max(4));
    if element_count > 0 {
        #[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
        let packed_handle = upload_u32(&client, &packed)?;
        // The CPU runtime (internal test feature) has no CUDA device to pin-copy to.
        #[cfg(all(test, feature = "_cubecl-cpu-test"))]
        let packed_handle = client.create_from_slice(u32::as_bytes(&packed));
        let column_word_handle = client.create_from_slice(u32::as_bytes(&column_word));
        let column_shift_handle = client.create_from_slice(u32::as_bytes(&column_shift));
        let lookup_handle = genotype_lookup(&client, read_options);
        // The output is (selected file columns) x (selected file rows).
        let (fast_count, slow_count) = memory_dimensions(columns.len(), rows.len(), kernel_is_f);
        let element_count_u32 = to_u32(element_count)?;
        // SAFETY: `decode` writes only `slow * fast_count + fast` for `fast < fast_count` and
        // `slow < slow_count`, so below `element_count` (which fits u32), the length of
        // `handle`; it reads `column_word`/`column_shift` at column indexes < `columns.len()`,
        // `packed` at word indexes < `rows.len() * span_words` = `packed_len` (checked above to
        // fit u32), and `lookup` at 2-bit codes <= 3 < 4, the lengths passed.
        unsafe {
            decode::launch_unchecked::<TVal>(
                &client,
                CubeCount::Static(
                    to_u32(fast_count)?.div_ceil(UNITS),
                    to_u32(slow_count.min(MAX_GRID_Y))?,
                    1,
                ),
                CubeDim::new_1d(UNITS),
                BufferArg::from_raw_parts(packed_handle, packed.len()),
                to_u32(span_words)?,
                BufferArg::from_raw_parts(column_word_handle, column_word.len()),
                BufferArg::from_raw_parts(column_shift_handle, column_shift.len()),
                to_u32(columns.len())?,
                to_u32(rows.len())?,
                BufferArg::from_raw_parts(lookup_handle, 4),
                BufferArg::from_raw_parts(handle.clone(), element_count_u32 as usize),
                kernel_is_f,
            );
        }
    }
    Ok(CubeClMatrix {
        client,
        handle,
        iid_count: iid_index.len(),
        sid_count: sid_index.len(),
        is_f: read_options.is_f,
        phantom: PhantomData,
    })
}

/// Reads, from each selected stored row, the `span_len` bytes starting at byte `span_start`, into
/// `span_words` u32 words per row. The caller has checked that the file holds `row_bytes` bytes for
/// each of its rows, every `row` is one of them, and `packed_len` = `rows.len() * span_words`.
fn read_packed(
    file: &mut File,
    rows: &[usize],
    row_bytes: u64,
    span_start: usize,
    span_len: usize,
    span_words: usize,
    packed_len: usize,
) -> std::io::Result<Vec<u32>> {
    let mut packed = vec![0u32; packed_len.max(1)];
    // Native-endian u32 words: the kernel's bit shifts assume little-endian bytes, as on every CUDA
    // host (x86-64 and ARM64).
    let bytes: &mut [u8] = bytemuck::cast_slice_mut(&mut packed);
    for (out_row, &row) in rows.iter().enumerate() {
        // As in the CPU reader. The file's length (which fits u64) bounds `offset`, so it cannot
        // overflow. `span_len <= span_words * 4`, so `destination + span_len <= packed_len * 4`, the
        // length of `bytes`.
        let offset = row as u64 * row_bytes + span_start as u64 + CB_HEADER_U64;
        file.seek(SeekFrom::Start(offset))?;
        let destination = out_row * span_words * 4;
        file.read_exact(&mut bytes[destination..destination + span_len])?;
    }
    Ok(packed)
}

/// A matrix's contiguous and strided dimensions (`fast_count`, `slow_count`) for its memory order.
fn memory_dimensions(row_count: usize, column_count: usize, is_f: bool) -> (usize, usize) {
    if is_f {
        (row_count, column_count)
    } else {
        (column_count, row_count)
    }
}

fn genotype_lookup<TVal: CubeClVal>(
    client: &Client,
    read_options: &ReadOptions<TVal>,
) -> cubecl::server::Handle {
    let values = set_up_two_bits_to_value(read_options.is_a1_counted, read_options.missing_value);
    client.create_from_slice(bytemuck::cast_slice(&values))
}

impl<TVal: CubeClVal> CubeClMatrix<TVal> {
    #[cfg(feature = "extension-module")]
    pub(crate) fn empty(rows: usize, columns: usize, is_f: bool) -> Result<Self, Box<Error>> {
        let element_count = element_count(rows, columns)?;
        // At most `u32::MAX` elements of at most 8 bytes, so this cannot overflow usize on the
        // 64-bit targets that `cubecl-cuda` supports.
        let byte_count = element_count * size_of::<TVal>();
        let client = client()?;
        // Keep a real pointer and one valid dummy element for empty kernel arguments.
        let handle = empty_dedicated(&client, byte_count.max(size_of::<TVal>()).max(4));
        Ok(Self {
            client,
            handle,
            iid_count: rows,
            sid_count: columns,
            is_f,
            phantom: PhantomData,
        })
    }

    /// The number of individuals and SNPs.
    ///
    /// See the [`CubeClMatrix`] example.
    #[must_use]
    pub fn dim(&self) -> (usize, usize) {
        (self.iid_count, self.sid_count)
    }

    /// True for F (column-major) order, false for C (row-major) order.
    ///
    /// See the [`CubeClMatrix`] example.
    #[must_use]
    pub fn is_f(&self) -> bool {
        self.is_f
    }

    /// Copy the matrix to the CPU, as an `ndarray` array in the same memory order.
    ///
    /// Copying a large matrix back gives up the point of a GPU read; keep the work on the GPU and
    /// copy only small results or checks.
    ///
    /// See the [`CubeClMatrix`] example.
    ///
    /// # Errors
    /// [`Error::Server`] if the GPU fails.
    pub fn to_cpu(&self) -> Result<nd::Array2<TVal>, Box<Error>> {
        // Every matrix comes from `read_cubecl` or `empty`, which check `element_count` (at most
        // `u32::MAX`) and allocate at least `element_count * size_of::<TVal>()` bytes. So this
        // product cannot overflow, and the slice of `bytes` below stays in bounds.
        let element_count = self.iid_count * self.sid_count;
        let shape = (self.iid_count, self.sid_count).set_f(self.is_f);
        if element_count == 0 {
            return Ok(nd::Array2::default(shape));
        }
        let bytes = self
            .client
            .read_one(self.handle.clone())
            .map_err(Error::from)?;
        let values = bytemuck::cast_slice::<u8, TVal>(
            &bytes[..element_count * core::mem::size_of::<TVal>()],
        );
        Ok(
            nd::Array2::from_shape_vec(shape, values.to_vec())
                .expect("shape matches element count"),
        )
    }

    /// The `CubeCL` client that owns this matrix's memory, for running `CubeCL` kernels on it.
    ///
    /// See the [`CubeClMatrix`] example.
    #[must_use]
    pub fn client(&self) -> &Client {
        &self.client
    }

    /// The matrix as `CubeCL`'s own tensor type: the same GPU memory with its shape, strides (in
    /// elements), and element type. Use it with [`CubeClMatrix::client`].
    ///
    /// The memory stays valid while this matrix or any clone of the handle exists. Kernels that
    /// write through the handle change this matrix's values, even though this method takes `&self`.
    ///
    /// See the [`CubeClMatrix`] example.
    #[must_use]
    pub fn tensor_handle(&self) -> TensorHandle {
        let strides = if self.is_f {
            [1, self.iid_count]
        } else {
            [self.sid_count, 1]
        };
        TensorHandle::new(
            self.handle.clone(),
            [self.iid_count, self.sid_count],
            strides,
            TVal::cube_type(),
        )
    }

    /// CUDA only: the matrix's device pointer, for handing it to another CUDA library (for
    /// example, `CuPy`) without a copy. Finishes this matrix's pending GPU work first.
    ///
    /// Use the pointer only while both the returned [`CudaPointer`] and this matrix exist; for that
    /// time, `CubeCL` does not move the allocation. Before dropping either one, or using the
    /// matrix again through `CubeCL`, finish all foreign CUDA work that accesses the pointer. The
    /// [`CudaPointer`] pins memory; it does not track other libraries' streams.
    ///
    /// # Errors
    /// [`Error::Server`] if the GPU fails.
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    ///
    /// let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
    /// let cubecl_matrix = ReadOptions::builder().f32().read_cubecl(&mut bed)?;
    /// let cuda_pointer = cubecl_matrix.cuda_pointer()?;
    /// assert_ne!(cuda_pointer.ptr(), 0);
    /// # use bed_reader::gpu::Error;
    /// # Ok::<(), Box<Error>>(())
    /// ```
    pub fn cuda_pointer(&self) -> Result<CudaPointer, Box<Error>> {
        use cubecl::cuda::CudaRuntime;
        use cubecl_runtime::runtime::Runtime;
        pollster::block_on(self.client.sync()).map_err(Error::from)?;
        let resource = self
            .client
            .get_resource::<<CudaRuntime as Runtime>::Server>(self.handle.clone())
            .map_err(Error::from)?;
        let ptr = resource.resource().ptr;
        Ok(CudaPointer {
            ptr,
            _resource: Box::new(resource),
        })
    }
}

/// CUDA only: a matrix's device pointer, from [`CubeClMatrix::cuda_pointer`].
///
/// Use it only while both it and its matrix exist; for that time, `CubeCL` does not move the
/// allocation. See [`CubeClMatrix::cuda_pointer`] for synchronizing with other CUDA libraries.
pub struct CudaPointer {
    ptr: u64,
    /// `CubeCL's` `ManagedResource`, which pins the allocation in place.
    _resource: Box<dyn core::any::Any + Send>,
}

impl<TVal> core::fmt::Debug for CubeClMatrix<TVal> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter
            .debug_struct("CubeClMatrix")
            .field("dim", &(self.iid_count, self.sid_count))
            .field("is_f", &self.is_f)
            .field("element", &core::any::type_name::<TVal>())
            .finish_non_exhaustive()
    }
}

impl core::fmt::Debug for CudaPointer {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter
            .debug_struct("CudaPointer")
            .field("ptr", &format_args!("{:#x}", self.ptr))
            .finish_non_exhaustive()
    }
}

impl CudaPointer {
    /// The device pointer.
    ///
    /// See [`CubeClMatrix::cuda_pointer`].
    #[must_use]
    pub fn ptr(&self) -> u64 {
        self.ptr
    }
}

/// Resolves each index with [`resolve_index`], the resolver every reader shares.
fn resolve(
    index: &[isize],
    count: usize,
    too_big: fn(isize) -> BedError,
) -> Result<Vec<usize>, BedError> {
    index
        .iter()
        .map(|&index| resolve_index(index, count, too_big))
        .collect()
}

fn to_u32(value: usize) -> Result<u32, Error> {
    u32::try_from(value).map_err(|_| Error::TooLarge(value as u64))
}

/// `rows * columns`, requiring both dimensions and the product to fit `u32`, as the kernels index.
fn element_count(rows: usize, columns: usize) -> Result<usize, Error> {
    let count = u64::from(to_u32(rows)?) * u64::from(to_u32(columns)?);
    let count = u32::try_from(count).map_err(|_| Error::TooLarge(count))?;
    Ok(count as usize)
}

/// Allocates `size` bytes outside `CubeCL`'s pools. Freed memory returns to CUDA at `CubeCL`'s
/// next allocation or cleanup instead of staying cached, so other GPU libraries in the process (for
/// example, `CuPy`) can use it. Read results and packed inputs can be large; pooling them would also
/// size the pool's later pages to them. Small temporaries stay pooled.
fn empty_dedicated(client: &Client, size: usize) -> cubecl::server::Handle {
    client.memory_dedicated_allocation(size, |size| client.empty(size))
}

/// Copies `words` to a new GPU buffer (through reused pinned host buffers on CUDA, which is about 4x
/// faster than `CubeCL's` own upload under `WSL2`).
#[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
fn upload_u32(client: &Client, words: &[u32]) -> Result<cubecl::server::Handle, Box<Error>> {
    let handle = empty_dedicated(client, core::mem::size_of_val(words));
    upload_pinned_cuda(client, &handle, bytemuck::cast_slice(words))?;
    Ok(handle)
}

/// Decodes one output element per unit, on a 2D grid: x covers the axis that is contiguous in the
/// output (its rows for F order, its columns for C order), so writes are coalesced and no division
/// is needed; y covers the other axis, striding by the grid height when it exceeds the grid limit.
///
/// The output has one row per selected file column and one column per selected file row (for mode
/// 1 files, individuals x SNPs). Element (`column`, `row`)'s 2-bit code sits in word
/// `row * span_words + column_word[column]` at bit `column_shift[column]`. `is_f` is compile-time,
/// so each memory order gets its own compiled kernel.
#[cube(launch_unchecked)]
fn decode<F: Numeric>(
    packed: &[u32],
    span_words: u32,
    column_word: &[u32],
    column_shift: &[u32],
    column_count: u32,
    row_count: u32,
    lookup: &[F],
    out: &mut [F],
    #[comptime] is_f: bool,
) {
    let fast = CUBE_POS_X * UNITS + UNIT_POS_X;
    let mut slow = CUBE_POS_Y;
    let mut fast_count = row_count;
    let mut slow_count = column_count;
    if is_f {
        fast_count = column_count;
        slow_count = row_count;
    }
    if fast < fast_count {
        while slow < slow_count {
            let mut column = slow;
            let mut row = fast;
            if is_f {
                column = fast;
                row = slow;
            }
            let word = packed[(row * span_words + column_word[column as usize]) as usize];
            let code = (word >> column_shift[column as usize]) & 3u32;
            out[(slow * fast_count + fast) as usize] = lookup[code as usize];
            // Avoid wrapping the u32 loop index at the largest supported dimensions.
            if slow_count - slow <= CUBE_COUNT_Y {
                break;
            }
            slow += CUBE_COUNT_Y;
        }
    }
}

/// CUDA only: copies `bytes` into the `CubeCL` buffer `handle` through two reused pinned host buffers,
/// so filling one overlaps copying the other.
#[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
fn upload_pinned_cuda(
    client: &Client,
    handle: &cubecl::server::Handle,
    bytes: &[u8],
) -> Result<(), Box<Error>> {
    use cubecl::cuda::CudaRuntime;
    use cubecl_runtime::runtime::Runtime;
    use cudarc::driver::result;

    // CubeCL can recycle allocations before their queued work finishes because
    // its own stream orders reuse. Our upload stream is external to that ordering:
    // finish earlier CubeCL uses before overwriting the recycled allocation.
    pollster::block_on(client.sync()).map_err(Error::from)?;
    let resource = client
        .get_resource::<<CudaRuntime as Runtime>::Server>(handle.clone())
        .map_err(Error::from)?;
    // `ptr` already includes the handle's offset within its CubeCL page.
    let device_ptr = resource.resource().ptr;
    assert!(
        bytes.len() as u64 <= resource.resource().size,
        "upload must fit its destination buffer"
    );
    let (staging_guard, staging) = cuda_staging()?;
    let copy_result = (|| -> Result<(), Error> {
        for (chunk_index, chunk) in bytes.chunks(CHUNK_BYTES).enumerate() {
            let slot = chunk_index % 2;
            // SAFETY: the event guards the pinned buffer until its previous copy finished; the
            // destination range lies inside the CubeCL buffer (asserted above), which `resource`
            // keeps alive and in place until the stream is synchronized below.
            unsafe {
                result::event::synchronize(staging.events[slot])?;
                let buffer = core::slice::from_raw_parts_mut(staging.buffers[slot], chunk.len());
                buffer.copy_from_slice(chunk);
                result::memcpy_htod_async(
                    device_ptr + (chunk_index * CHUNK_BYTES) as u64,
                    buffer,
                    staging.stream,
                )?;
                result::event::record(staging.events[slot], staging.stream)?;
            }
        }
        Ok(())
    })();
    // Wait even when a copy failed: copies already queued still write into the CubeCL buffer and
    // read the pinned buffers, so neither may be released or reused until they finish.
    // SAFETY: the stream is valid; waiting also makes the data visible to CubeCL's later launches.
    let sync_result = unsafe { result::stream::synchronize(staging.stream) };
    drop(staging_guard);
    drop(resource);
    copy_result?;
    sync_result.map_err(Error::from)?;
    Ok(())
}

#[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
const CHUNK_BYTES: usize = 16 << 20;

/// Two reused pinned host buffers, their events, and a stream, for chunked copies. A `Copy` of
/// the handles in [`StagingHandles`], valid only while its lock is held.
#[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
#[derive(Clone, Copy)]
struct CudaStaging {
    stream: cudarc::driver::sys::CUstream,
    buffers: [*mut u8; 2],
    events: [cudarc::driver::sys::CUevent; 2],
}

/// The process-wide staging handles, each created on first use and never released. Each one is
/// stored as soon as it exists, so a failed initialization keeps what it created and the next
/// call resumes from there instead of leaking it and starting over.
#[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
struct StagingHandles {
    /// Device 0's primary context, retained once for the life of the process.
    context: Option<cudarc::driver::sys::CUcontext>,
    stream: Option<cudarc::driver::sys::CUstream>,
    buffers: [Option<*mut u8>; 2],
    events: [Option<cudarc::driver::sys::CUevent>; 2],
}
// SAFETY: the raw handles are only used while holding the mutex in `cuda_staging`.
#[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
unsafe impl Send for StagingHandles {}

/// Locks the shared staging handles, creating any that are missing, and returns the lock with a
/// copy of the handles. On error, handles created so far stay stored for the next call.
#[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
fn cuda_staging() -> Result<(std::sync::MutexGuard<'static, StagingHandles>, CudaStaging), Error> {
    use cudarc::driver::{result, sys, DriverError};

    /// Returns the stored handle, or creates and stores it.
    fn stored<T: Copy>(
        handle: &mut Option<T>,
        create: impl FnOnce() -> Result<T, DriverError>,
    ) -> Result<T, DriverError> {
        if let Some(handle) = *handle {
            return Ok(handle);
        }
        let created = create()?;
        *handle = Some(created);
        Ok(created)
    }

    static STAGING: std::sync::Mutex<StagingHandles> = std::sync::Mutex::new(StagingHandles {
        context: None,
        stream: None,
        buffers: [None; 2],
        events: [None; 2],
    });
    // Nothing panics while this lock is held except a bug, so poisoning is not expected.
    let mut guard = STAGING.lock().expect("staging lock");
    let handles = &mut *guard;
    // SAFETY: device 0 exists (CubeCL already uses it); its primary context is the one CubeCL
    // uses, retained once here and never released. The context must be current on whichever
    // thread copies, so it is made current on every call before creating or using the stream,
    // buffers, and events, which are likewise created once and never released.
    let staging = unsafe {
        let context = stored(&mut handles.context, || {
            result::init()?;
            result::primary_ctx::retain(result::device::get(0)?)
        })?;
        result::ctx::set_current(context)?;
        let stream = stored(&mut handles.stream, || {
            result::stream::create(result::stream::StreamKind::NonBlocking)
        })?;
        let mut buffers = [core::ptr::null_mut(); 2];
        for (buffer, handle) in buffers.iter_mut().zip(&mut handles.buffers) {
            *buffer = stored(handle, || {
                Ok(result::malloc_host(CHUNK_BYTES, 0)?.cast::<u8>())
            })?;
        }
        let mut events = [core::ptr::null_mut(); 2];
        for (event, handle) in events.iter_mut().zip(&mut handles.events) {
            *event = stored(handle, || {
                result::event::create(sys::CUevent_flags::CU_EVENT_DISABLE_TIMING)
            })?;
        }
        CudaStaging {
            stream,
            buffers,
            events,
        }
    };
    Ok((guard, staging))
}

#[cfg(test)]
mod size_tests {
    //! Size limits, which need no GPU.
    use super::{element_count, Error};

    #[test]
    fn element_count_requires_u32_dimensions_and_product() {
        assert_eq!(element_count(3, 4).ok(), Some(12));
        assert_eq!(element_count(0, 5).ok(), Some(0));
        assert!(matches!(
            element_count(65_536, 65_536),
            Err(Error::TooLarge(4_294_967_296))
        ));
        assert!(matches!(
            element_count(1 << 32, 0),
            Err(Error::TooLarge(4_294_967_296))
        ));
    }
}

#[cfg(test)]
mod tests {
    //! Compare the GPU reader and standardizer with the CPU code, exactly (NaN matches NaN). These
    //! need a GPU; without one, they fail (they do not skip).
    use std::{error::Error as StandardError, marker::PhantomData};

    use ndarray as nd;

    use super::{
        client, memory_cleanup,
        standardize::{standardize_beta, standardize_unit, CubeClFloat},
        CubeClMatrix, CubeClVal,
    };
    use crate::{impute_and_zero_mean_snps, sample_bed_file, Bed, Dist, ReadOptions};

    const FILES: [&str; 6] = [
        "small.bed",
        "smallmode0.bed",
        "some_missing.bed",
        "plink_sim_10s_100v_10pmiss.bed",
        "sparse.bed",
        "toydata.5chrom.bed",
    ];

    /// All, reversed, negatives with repeats, strided, one, and none.
    fn index_cases(count: usize) -> Vec<Vec<isize>> {
        let count = isize::try_from(count).expect("fits");
        vec![
            (0..count).collect(),
            (0..count).rev().collect(),
            vec![-1, 0, count - 1, -count, 0],
            (0..count).step_by(3).collect(),
            vec![count / 2],
            vec![],
        ]
    }

    /// Bitwise equal, including the canonical NaNs produced by the BED reader.
    fn assert_same<T: CubeClVal + PartialEq>(
        got: &nd::Array2<T>,
        expected: &nd::Array2<T>,
        context: &str,
    ) {
        assert_eq!(got.dim(), expected.dim(), "{context}");
        for (gpu_value, cpu_value) in got.iter().zip(expected.iter()) {
            assert_eq!(
                bytemuck::bytes_of(gpu_value),
                bytemuck::bytes_of(cpu_value),
                "{context}: GPU {gpu_value:?} != CPU {cpu_value:?}"
            );
        }
    }

    fn check_read<T: CubeClVal>(
        bed: &mut Bed,
        iid_index: &[isize],
        sid_index: &[isize],
        is_a1_counted: bool,
        is_f: bool,
    ) -> Result<nd::Array2<T>, Box<dyn StandardError>> {
        let mut options = ReadOptions::<T>::builder();
        options
            .iid_index(iid_index)
            .sid_index(sid_index)
            .is_a1_counted(is_a1_counted)
            .is_f(is_f);
        let context = format!("iid {iid_index:?} sid {sid_index:?} a1 {is_a1_counted} f {is_f}");
        let expected = options.read(bed)?;
        let matrix = options.read_cubecl(bed)?;
        assert_eq!(matrix.is_f(), is_f);
        assert_eq!(matrix.dim(), expected.dim());
        assert_same(&matrix.to_cpu()?, &expected, &context);
        Ok(expected)
    }

    fn check<T: CubeClFloat>(
        bed: &mut Bed,
        iid_index: &[isize],
        sid_index: &[isize],
        is_a1_counted: bool,
        is_f: bool,
    ) -> Result<(), Box<dyn StandardError>> {
        let expected = check_read::<T>(bed, iid_index, sid_index, is_a1_counted, is_f)?;
        let mut options = ReadOptions::<T>::builder();
        options
            .iid_index(iid_index)
            .sid_index(sid_index)
            .is_a1_counted(is_a1_counted)
            .is_f(is_f);
        let context = format!("iid {iid_index:?} sid {sid_index:?} a1 {is_a1_counted} f {is_f}");

        for dist in [
            Dist::Unit,
            Dist::Beta { a: 1.0, b: 25.0 },
            Dist::Beta { a: 2.0, b: 3.0 },
        ] {
            let mut cpu_val = expected.clone();
            let mut cpu_stats = nd::Array2::<T>::zeros((sid_index.len(), 2));
            let cpu_result = impute_and_zero_mean_snps(
                &mut cpu_val.view_mut(),
                &dist,
                true,
                false,
                &mut cpu_stats.view_mut(),
            );
            let mut matrix = options.read_cubecl(bed)?;
            let gpu_result = match dist {
                Dist::Unit => standardize_unit(&mut matrix),
                Dist::Beta { a, b } => standardize_beta(&mut matrix, a, b),
            };
            match (cpu_result, gpu_result) {
                (Ok(()), Ok(gpu_stats)) => {
                    assert_same(&gpu_stats, &cpu_stats, &format!("{context} stats {dist:?}"));
                    assert_same(
                        &matrix.to_cpu()?,
                        &cpu_val,
                        &format!("{context} values {dist:?}"),
                    );
                }
                (Err(cpu_error), Err(gpu_error)) => {
                    assert_eq!(
                        gpu_error.to_string(),
                        cpu_error.to_string(),
                        "{context} {dist:?}"
                    );
                }
                (cpu_result, gpu_result) => panic!(
                    "{context} {dist:?}: CPU {cpu_result:?}, GPU {:?}",
                    gpu_result.map(|_| ())
                ),
            }
        }
        Ok(())
    }

    #[test]
    fn gpu_matches_cpu() -> Result<(), Box<dyn StandardError>> {
        for file in FILES {
            let mut bed = Bed::new(sample_bed_file(file)?)?;
            let (iid_count, sid_count) = bed.dim()?;
            for iid_index in index_cases(iid_count) {
                for sid_index in index_cases(sid_count) {
                    for is_a1_counted in [true, false] {
                        for is_f in [true, false] {
                            check_read::<i8>(
                                &mut bed,
                                &iid_index,
                                &sid_index,
                                is_a1_counted,
                                is_f,
                            )?;
                            check::<f32>(&mut bed, &iid_index, &sid_index, is_a1_counted, is_f)?;
                            check::<f64>(&mut bed, &iid_index, &sid_index, is_a1_counted, is_f)?;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn gpu_entry_points_match_cpu() -> Result<(), Box<dyn StandardError>> {
        let mut bed = Bed::new(sample_bed_file("small.bed")?)?;

        // Bed::read_cubecl matches Bed::read (all data, default options).
        let all = bed.read_cubecl::<f32>()?;
        assert_same(&all.to_cpu()?, &bed.read::<f32>()?, "Bed::read_cubecl");
        assert!(all.is_f());

        // The builder and built-options entry points match the CPU read with the same options.
        let expected = ReadOptions::builder()
            .i8()
            .c()
            .count_a2()
            .sid_index([2, 0])
            .read(&mut bed)?;
        let from_builder = ReadOptions::builder()
            .i8()
            .c()
            .count_a2()
            .sid_index([2, 0])
            .read_cubecl(&mut bed)?;
        let read_options = ReadOptions::builder()
            .i8()
            .c()
            .count_a2()
            .sid_index([2, 0])
            .build()?;
        let from_options = bed.read_cubecl_with_options(&read_options)?;
        assert_eq!(from_builder.to_cpu()?, expected);
        assert_eq!(from_options.to_cpu()?, expected);
        assert!(!from_builder.is_f());

        // A later option overrides an earlier one, as on the CPU.
        let count_a1 = ReadOptions::builder()
            .f64()
            .count_a2()
            .count_a1()
            .f()
            .read_cubecl(&mut bed)?;
        assert_same(
            &count_a1.to_cpu()?,
            &ReadOptions::builder().f64().read(&mut bed)?,
            "count_a1, f",
        );
        Ok(())
    }

    #[test]
    fn gpu_read_options_match_cpu() -> Result<(), Box<dyn StandardError>> {
        fn check_options<TVal: CubeClVal + PartialEq>(
            missing_value: TVal,
        ) -> Result<(), Box<dyn StandardError>> {
            let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
            for is_f in [true, false] {
                for is_a1_counted in [true, false] {
                    let mut options = ReadOptions::<TVal>::builder();
                    options
                        .missing_value(missing_value)
                        .iid_index([true, false, true])
                        .sid_index(1..4)
                        .is_f(is_f)
                        .is_a1_counted(is_a1_counted);
                    let expected = options.read(&mut bed)?;
                    let matrix = options
                        .num_threads(1)
                        .max_concurrent_requests(1)
                        .max_chunk_bytes(1)
                        .read_cubecl(&mut bed)?;
                    let actual = matrix.to_cpu()?;
                    assert_same(&actual, &expected, "custom missing value and read options");
                    assert_eq!(actual[[0, 1]], missing_value);
                    assert_eq!(matrix.is_f(), is_f);
                }
            }
            Ok(())
        }
        check_options(-7i8)?;
        check_options(-3.25f32)?;
        check_options(-9.5f64)?;
        Ok(())
    }

    #[test]
    fn gpu_tensor_handle_describes_the_matrix() -> Result<(), Box<dyn StandardError>> {
        let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
        for is_f in [true, false] {
            let matrix = ReadOptions::builder()
                .f32()
                .is_f(is_f)
                .read_cubecl(&mut bed)?;
            let tensor_handle = matrix.tensor_handle();
            assert_eq!(tensor_handle.shape().as_slice(), &[3, 4]);
            let expected_strides: [usize; 2] = if is_f { [1, 3] } else { [4, 1] };
            assert_eq!(&tensor_handle.strides()[..], &expected_strides);
            // The handle is the matrix's own memory: reading it through the client gives the values.
            let bytes = matrix
                .client()
                .read_one(tensor_handle.handle)
                .map_err(super::Error::from)?;
            let values: &[f32] = bytemuck::cast_slice(&bytes[..12 * 4]);
            let host = matrix.to_cpu()?;
            let expected: Vec<f32> = if is_f {
                host.t().iter().copied().collect()
            } else {
                host.iter().copied().collect()
            };
            assert_same(
                &nd::Array2::from_shape_vec((1, 12), values.to_vec()).expect("12 values"),
                &nd::Array2::from_shape_vec((1, 12), expected).expect("12 values"),
                "tensor handle values",
            );
        }
        Ok(())
    }

    #[test]
    fn gpu_debug_describes_the_matrix() -> Result<(), Box<dyn StandardError>> {
        let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
        let matrix = ReadOptions::builder()
            .i8()
            .c()
            .sid_index(1..3)
            .read_cubecl(&mut bed)?;
        assert_eq!(
            format!("{matrix:?}"),
            r#"CubeClMatrix { dim: (3, 2), is_f: false, element: "i8", .. }"#
        );
        Ok(())
    }

    // CUDA only: CubeCL's CPU runtime (internal feature `_cubecl-cpu-test`) has no device pointers.
    #[cfg(not(feature = "_cubecl-cpu-test"))]
    #[test]
    fn gpu_cuda_pointer_debug_shows_the_address() -> Result<(), Box<dyn StandardError>> {
        let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
        let matrix = bed.read_cubecl::<f32>()?;
        let cuda_pointer = matrix.cuda_pointer()?;
        assert_eq!(
            format!("{cuda_pointer:?}"),
            format!("CudaPointer {{ ptr: {:#x}, .. }}", cuda_pointer.ptr())
        );
        Ok(())
    }

    // CUDA only: the CPU runtime (internal feature `_cubecl-cpu-test`) uploads without staging.
    #[cfg(not(feature = "_cubecl-cpu-test"))]
    #[test]
    fn gpu_upload_crosses_staging_chunks() -> Result<(), Box<dyn StandardError>> {
        // The sample files fit one staging chunk. Two and a half chunks exercise both pinned
        // buffers, reuse of the first, and a short final chunk.
        let word_count = (5 * super::CHUNK_BYTES / 2) / size_of::<u32>() + 3;
        let words: Vec<u32> = (0..u32::try_from(word_count)?)
            .map(|index| index.wrapping_mul(2_654_435_761))
            .collect();
        let client = super::client()?;
        let handle = super::upload_u32(&client, &words)?;
        let bytes = client.read_one(handle).map_err(super::Error::from)?;
        let uploaded: &[u32] = bytemuck::cast_slice(&bytes[..word_count * 4]);
        // Report the first wrong word rather than printing millions of them.
        let first_mismatch = uploaded
            .iter()
            .zip(&words)
            .position(|(got, expected)| got != expected);
        assert_eq!(first_mismatch, None);
        Ok(())
    }

    #[test]
    fn gpu_memory_cleanup_runs() -> Result<(), Box<dyn StandardError>> {
        let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
        let matrix = ReadOptions::builder().f32().read_cubecl(&mut bed)?;
        let expected = matrix.to_cpu()?;
        memory_cleanup()?; // must not free memory that a live matrix uses
        assert_same(&matrix.to_cpu()?, &expected, "after cleanup");
        drop(matrix);
        memory_cleanup()?;
        Ok(())
    }

    #[test]
    fn gpu_standardize_accepts_non_genotype_values() -> Result<(), Box<dyn StandardError>> {
        // A custom missing value (-1) is observed data, not NaN. It must take
        // the general floating-point path rather than the genotype optimization.
        let mut bed = Bed::new(sample_bed_file("some_missing.bed")?)?;
        let mut matrix = ReadOptions::builder()
            .missing_value(-1.0f32)
            .read_cubecl(&mut bed)?;
        let original = matrix.to_cpu()?;
        let observed_count =
            f64::from(u32::try_from(original.nrows()).expect("GPU dimensions fit u32"));
        let expected_mean = original
            .column(0)
            .iter()
            .map(|&value| f64::from(value))
            .sum::<f64>()
            / observed_count;
        let expected_variance = original
            .column(0)
            .iter()
            .map(|&value| (f64::from(value) - expected_mean).powi(2))
            .sum::<f64>()
            / observed_count;
        let statistics = standardize_unit(&mut matrix)?;
        assert!(
            (f64::from(statistics[[0, 0]]) - expected_mean).abs() < 4.0 * f64::from(f32::EPSILON)
        );
        assert!(
            (f64::from(statistics[[0, 1]]) - expected_variance.sqrt()).abs()
                < 4.0 * f64::from(f32::EPSILON)
        );
        Ok(())
    }

    #[test]
    fn gpu_standardize_constant_fractional_column_is_zero() -> Result<(), Box<dyn StandardError>> {
        // Three f64 0.1s average to 0.10000000000000002, so without exact constant detection
        // the column got a tiny positive variance and standardized to -1 instead of 0.
        constant_fractional_column_is_zero(0.1f64)?;
        constant_fractional_column_is_zero(0.1f32)?;
        Ok(())
    }

    fn constant_fractional_column_is_zero<TVal: CubeClFloat>(
        value: TVal,
    ) -> Result<(), Box<dyn StandardError>> {
        // Four individuals by two SNPs, in F order: a constant column with one missing value,
        // then a varying fractional column.
        let [quarter, half] = [0.25, 0.5].map(|number| TVal::from_f64(number).expect("converts"));
        let values = [
            value,
            value,
            TVal::nan(),
            value,
            quarter,
            half,
            half,
            quarter,
        ];
        let client = client()?;
        let handle = client.create_from_slice(bytemuck::cast_slice(&values));
        let mut matrix = CubeClMatrix::<TVal> {
            client,
            handle,
            iid_count: 4,
            sid_count: 2,
            is_f: true,
            phantom: PhantomData,
        };
        let statistics = standardize_unit(&mut matrix)?;
        assert_eq!(statistics[[0, 0]], value);
        assert!(num_traits::Float::is_infinite(statistics[[0, 1]]));
        let standardized = matrix.to_cpu()?;
        assert_eq!(standardized.column(0).to_vec(), vec![TVal::zero(); 4]);
        Ok(())
    }

    #[test]
    fn gpu_index_errors_match_cpu() -> Result<(), Box<dyn StandardError>> {
        // Mode 0 (individual-major) files swap axes internally; errors must not.
        for name in ["small.bed", "smallmode0.bed"] {
            let mut bed = Bed::new(sample_bed_file(name)?)?;
            // The two files have different shapes, so indexes just out of range differ.
            let iid_count = isize::try_from(bed.iid_count()?)?;
            let sid_count = isize::try_from(bed.sid_count()?)?;
            for (iid_index, sid_index, axis) in [
                (vec![iid_count], vec![0isize], "individual"),
                (vec![-iid_count - 1], vec![0], "individual"),
                (vec![0], vec![sid_count], "SNP"),
                (vec![0], vec![-sid_count - 1], "SNP"),
            ] {
                let mut options = ReadOptions::<f32>::builder();
                options.iid_index(&iid_index).sid_index(&sid_index);
                let cpu_error = options.read(&mut bed).expect_err("CPU should fail");
                let gpu_error = options.read_cubecl(&mut bed).expect_err("GPU should fail");
                assert_eq!(gpu_error.to_string(), cpu_error.to_string());
                assert!(
                    cpu_error
                        .to_string()
                        .starts_with(&format!("Index to {axis} ")),
                    "{name}: {cpu_error}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn gpu_rejects_counts_too_big_for_any_file() -> Result<(), Box<dyn StandardError>> {
        // Only a caller-given count can exceed isize::MAX. It fails before any GPU work.
        let too_big = usize::try_from(isize::MAX)? + 1;
        let mut bed = Bed::builder(sample_bed_file("small.bed")?)
            .iid_count(too_big)
            .sid_count(4)
            .build()?;
        let error = bed.read_cubecl::<f32>().expect_err("count too big");
        assert_eq!(
            error.to_string(),
            format!("Indexes ({too_big},4) too big for files")
        );
        Ok(())
    }
}
