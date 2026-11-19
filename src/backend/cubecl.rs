//! The `CubeCL` GPU backend (Cargo feature `cubecl-cuda`), re-exported by [`backend`](super).
use core::{fmt::Debug, marker::PhantomData, mem::size_of};
use std::io::{Read, Seek, SeekFrom};

use cubecl::{prelude::*, std::tensor::TensorHandle};
use ndarray as nd;
use ndarray::ShapeBuilder;

use super::{
    sealed::{Sealed, Value},
    Backend, Error,
};
use crate::{
    open_and_check, path_ref_to_string, Bed, BedError, BedErrorPlus, BedVal, Hold, ReadOptions,
    CB_HEADER_U64,
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
// TODO0 experiment (2026-11-19): decide go/no-go before release; remove if no-go.
#[cfg(all(test, feature = "_cubecl-cpu-test"))]
const UNITS: u32 = 4;
/// Cubes along a grid's y axis are limited to this many; kernels stride over the rest.
const MAX_GRID_Y: usize = 65_535;

/// GPU device 0, through this build's selected `CubeCL` runtime: `read` returns a
/// [`CubeClMatrix`].
///
/// # Example
///
/// See [`CubeClMatrix`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CubeCl;

impl Sealed for CubeCl {}

impl Backend for CubeCl {
    type Matrix<TVal: BedVal> = CubeClMatrix<TVal>;
}

impl CubeCl {
    /// The `CubeCL` client of the GPU that reads go to.
    ///
    /// A [`CubeClMatrix`] carries its own client ([`CubeClMatrix::client`]); use this one for work
    /// that is not tied to a matrix, such as [`CubeCl::memory_cleanup`].
    ///
    /// # Errors
    /// [`Error::Unavailable`] if there is no GPU.
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use bed_reader::backend::CubeCl;
    ///
    /// let client = CubeCl::client()?;
    /// println!("{}", client.name());
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn client() -> Result<Client, Box<BedErrorPlus>> {
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
    /// use bed_reader::{backend::CubeCl, Bed, ReadOptions, sample_bed_file};
    ///
    /// let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
    /// let cubecl_matrix = ReadOptions::builder().f32().cubecl().read(&mut bed)?;
    /// drop(cubecl_matrix);
    /// CubeCl::memory_cleanup()?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn memory_cleanup() -> Result<(), Box<BedErrorPlus>> {
        let client = Self::client()?;
        // Finish queued work first, so buffers that it still uses are not counted as in use later.
        pollster::block_on(client.sync()).map_err(Error::from)?;
        client.memory_cleanup().map_err(Error::from)?;
        Ok(())
    }
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

/// A genotype matrix (individuals x SNPs) in GPU memory, from a `.cubecl()` read.
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
/// let cubecl_matrix = ReadOptions::builder().f32().f().cubecl().read(&mut bed)?;
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
///     .cubecl().read(&mut bed)?;
/// assert_eq!(selected.dim(), (2, 3));
/// assert!(!selected.is_f());
/// # use bed_reader::BedErrorPlus;
/// # Ok::<(), Box<BedErrorPlus>>(())
/// ```
pub struct CubeClMatrix<TVal> {
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
    read_options: &ReadOptions<TVal, CubeCl>,
) -> Result<CubeClMatrix<TVal>, Box<BedErrorPlus>> {
    let iid_count = bed.iid_count()?;
    let sid_count = bed.sid_count()?;
    let iid_hold = Hold::new(&read_options.iid_index, iid_count)?;
    let sid_hold = Hold::new(&read_options.sid_index, sid_count)?;
    let iid_index = resolve(iid_hold.as_ref(), iid_count, BedError::IidIndexTooBig)?;
    let sid_index = resolve(sid_hold.as_ref(), sid_count, BedError::SidIndexTooBig)?;
    let element_count = element_count(iid_index.len(), sid_index.len())?;
    let client = CubeCl::client()?;

    let (buf_reader, header) = open_and_check(&bed.path)?;
    let mut file = buf_reader.into_inner();
    // Mode 1 stores each SNP's individuals together (SNP-major); mode 0 stores each individual's
    // SNPs together. For mode 0, swap the roles of individuals and SNPs and flip the memory
    // order: the transposed result has the requested layout (as the CPU reader does).
    let (row_count, rows, col_count, cols, kernel_is_f) = match header[2] {
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
    let row_bytes = col_count.div_ceil(4);
    let expected_len = (row_bytes as u64)
        .checked_mul(row_count as u64)
        .and_then(|len| len.checked_add(CB_HEADER_U64))
        .ok_or(BedError::IndexesTooBigForFiles(iid_count, sid_count))?;
    if file.metadata()?.len() != expected_len {
        Err(BedError::IllFormed(path_ref_to_string(&bed.path)))?;
    }

    // Within each stored row, read only the bytes spanning the selected columns, packed into
    // u32 words. Column `c`'s 2-bit code sits in word `col_word[c]` at bit `col_shift[c]`.
    let (span_start, span_len) = match (cols.iter().min(), cols.iter().max()) {
        (Some(&min), Some(&max)) => (min / 4, max / 4 + 1 - min / 4),
        _ => (0, 0),
    };
    let span_words = span_len.div_ceil(4);
    let col_word: Vec<u32> = cols
        .iter()
        .map(|&c| to_u32((c / 4 - span_start) / 4))
        .collect::<Result<_, _>>()?;
    let col_shift: Vec<u32> = cols
        .iter()
        .map(|&c| to_u32(((c / 4 - span_start) % 4) * 8 + (c % 4) * 2))
        .collect::<Result<_, _>>()?;
    let mut packed = vec![0u32; (rows.len() * span_words).max(1)];
    {
        let bytes: &mut [u8] = bytemuck::cast_slice_mut(&mut packed);
        for (out_row, &row) in rows.iter().enumerate() {
            let pos = CB_HEADER_U64 + (row * row_bytes + span_start) as u64;
            file.seek(SeekFrom::Start(pos))?;
            let dst = out_row * span_words * 4;
            file.read_exact(&mut bytes[dst..dst + span_len])?;
        }
    }

    // Empty kernel arguments advertise one dummy element; it must fit even for f64.
    let handle = empty_dedicated(&client, (element_count.max(1) * size_of::<TVal>()).max(4));
    if element_count > 0 {
        #[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
        let packed_handle = upload_u32(&client, &packed)?;
        // The CPU runtime (internal test feature) has no CUDA device to pin-copy to.
        #[cfg(all(test, feature = "_cubecl-cpu-test"))]
        let packed_handle = client.create_from_slice(u32::as_bytes(&packed));
        let col_word_handle = client.create_from_slice(u32::as_bytes(&col_word));
        let col_shift_handle = client.create_from_slice(u32::as_bytes(&col_shift));
        let lut_handle = genotype_lookup(&client, read_options);
        let (fast_count, slow_count) = memory_dimensions(cols.len(), rows.len(), kernel_is_f);
        let element_count_u32 = to_u32(element_count)?;
        // The kernel indexes `packed` with u32. A few far-apart columns can make it much larger
        // than the result, so check it separately from `element_count`.
        to_u32(packed.len())?;
        // SAFETY: `decode` writes only `slow * fast_count + fast` for `fast < fast_count` and
        // `slow < slow_count`, so below `element_count` (which fits u32), the length of
        // `handle`; it reads `col_word`/`col_shift` at column indexes < `cols.len()` and
        // `packed` at word indexes < `rows.len() * span_words` = `packed.len()` (checked above to
        // fit u32), the lengths passed.
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
                BufferArg::from_raw_parts(col_word_handle, col_word.len()),
                BufferArg::from_raw_parts(col_shift_handle, col_shift.len()),
                to_u32(cols.len())?,
                to_u32(rows.len())?,
                u32::from(!kernel_is_f),
                BufferArg::from_raw_parts(lut_handle, 4),
                BufferArg::from_raw_parts(handle.clone(), element_count_u32 as usize),
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

/// The contiguous and strided dimensions for the selected memory order.
fn memory_dimensions(iid_count: usize, sid_count: usize, is_f: bool) -> (usize, usize) {
    if is_f {
        (iid_count, sid_count)
    } else {
        (sid_count, iid_count)
    }
}

fn genotype_lookup<TVal: CubeClVal>(
    client: &Client,
    read_options: &ReadOptions<TVal, CubeCl>,
) -> cubecl::server::Handle {
    let missing = read_options.missing_value;
    let values: [TVal; 4] = if read_options.is_a1_counted {
        [genotype(2), missing, genotype(1), genotype(0)]
    } else {
        [genotype(0), missing, genotype(1), genotype(2)]
    };
    client.create_from_slice(bytemuck::cast_slice(&values))
}

impl<TVal: CubeClVal> CubeClMatrix<TVal> {
    #[cfg(all(feature = "extension-module", feature = "cubecl-cuda"))]
    pub(crate) fn empty(
        rows: usize,
        columns: usize,
        is_f: bool,
    ) -> Result<Self, Box<BedErrorPlus>> {
        let element_count = element_count(rows, columns)?;
        // At most `u32::MAX` elements of at most 8 bytes, so this cannot overflow usize on the
        // 64-bit targets that `cubecl-cuda` supports.
        let byte_count = element_count * size_of::<TVal>();
        let client = CubeCl::client()?;
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
    pub fn to_cpu(&self) -> Result<nd::Array2<TVal>, Box<BedErrorPlus>> {
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
    /// The memory stays valid while this matrix or any clone of the handle exists.
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
    /// The returned [`CudaPointer`] keeps the allocation alive and in place; use the pointer only
    /// while both it and this matrix exist. Before releasing them, or using the
    /// buffer again through `CubeCL`, finish all foreign CUDA work that accesses
    /// the pointer. This guard pins memory; it does not track other libraries' streams.
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
    /// let cubecl_matrix = ReadOptions::builder().f32().cubecl().read(&mut bed)?;
    /// let cuda_pointer = cubecl_matrix.cuda_pointer()?;
    /// assert_ne!(cuda_pointer.ptr(), 0);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn cuda_pointer(&self) -> Result<CudaPointer, Box<BedErrorPlus>> {
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

/// CUDA only: a matrix's device pointer, from [`CubeClMatrix::cuda_pointer`]. While it exists, `CubeCL`
/// keeps the allocation alive and never relocates it.
pub struct CudaPointer {
    ptr: u64,
    /// `CubeCL's` `ManagedResource`, which pins the allocation in place.
    _resource: Box<dyn core::any::Any + Send>,
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

/// Resolves Python-style indexes (negative counts from the end), with the CPU reader's errors.
fn resolve(
    index: &[isize],
    count: usize,
    too_big: fn(isize) -> BedError,
) -> Result<Vec<usize>, BedError> {
    let count_isize =
        isize::try_from(count).map_err(|_| BedError::IndexesTooBigForFiles(count, 0))?;
    index
        .iter()
        .map(|&i| {
            if (0..count_isize).contains(&i) {
                Ok(i.unsigned_abs())
            } else if (-count_isize..0).contains(&i) {
                Ok(count - i.unsigned_abs())
            } else {
                Err(too_big(i))
            }
        })
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

fn genotype<TVal: CubeClVal>(value: i8) -> TVal {
    <TVal as From<i8>>::from(value)
}

/// Copies `words` to a new GPU buffer (through reused pinned host buffers on CUDA, which is about 4x
/// faster than `CubeCL's` own upload under `WSL2`).
#[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
fn upload_u32(client: &Client, words: &[u32]) -> Result<cubecl::server::Handle, Box<BedErrorPlus>> {
    let handle = empty_dedicated(client, core::mem::size_of_val(words));
    upload_pinned_cuda(client, &handle, bytemuck::cast_slice(words))?;
    Ok(handle)
}

/// Decodes one output element per unit, on a 2D grid: x covers the axis that is contiguous in the
/// output (individuals for F order, SNPs for C order), so writes are coalesced and no division is
/// needed; y covers the other axis, striding by the grid height when it exceeds the grid limit.
/// Element (`iid_i`, `sid_j`)'s 2-bit code sits in word `sid_j * span_words + iid_word[iid_i]` at
/// bit `iid_shift[iid_i]`.
#[cube(launch_unchecked)]
fn decode<F: Numeric>(
    packed: &[u32],
    span_words: u32,
    iid_word: &[u32],
    iid_shift: &[u32],
    iid_count: u32,
    sid_count: u32,
    c_order: u32,
    lut: &[F],
    out: &mut [F],
) {
    let fast = CUBE_POS_X * UNITS + UNIT_POS_X;
    let mut slow = CUBE_POS_Y;
    let mut fast_count = iid_count;
    let mut slow_count = sid_count;
    if c_order == 1u32 {
        fast_count = sid_count;
        slow_count = iid_count;
    }
    if fast < fast_count {
        while slow < slow_count {
            let mut iid_i = fast;
            let mut sid_j = slow;
            if c_order == 1u32 {
                iid_i = slow;
                sid_j = fast;
            }
            let word = packed[(sid_j * span_words + iid_word[iid_i as usize]) as usize];
            let code = (word >> iid_shift[iid_i as usize]) & 3u32;
            out[(slow * fast_count + fast) as usize] = lut[code as usize];
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
) -> Result<(), Box<BedErrorPlus>> {
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
    let device_ptr = resource.resource().ptr;
    let staging_guard = cuda_staging()?;
    let staging = staging_guard
        .as_ref()
        .expect("cuda_staging returns only after initializing the staging buffers");
    for (chunk_index, chunk) in bytes.chunks(CHUNK_BYTES).enumerate() {
        let slot = chunk_index % 2;
        // SAFETY: the event guards the pinned buffer until its previous copy finished; the
        // destination range lies inside the CubeCL buffer, which `resource` keeps alive and in place.
        unsafe {
            result::event::synchronize(staging.events[slot]).map_err(Error::from)?;
            let buffer = core::slice::from_raw_parts_mut(staging.buffers[slot], chunk.len());
            buffer.copy_from_slice(chunk);
            result::memcpy_htod_async(
                device_ptr + (chunk_index * CHUNK_BYTES) as u64,
                buffer,
                staging.stream,
            )
            .map_err(Error::from)?;
            result::event::record(staging.events[slot], staging.stream).map_err(Error::from)?;
        }
    }
    // SAFETY: the stream is valid; waiting makes the data visible to CubeCL's later launches.
    unsafe { result::stream::synchronize(staging.stream).map_err(Error::from)? };
    drop(staging_guard);
    drop(resource);
    Ok(())
}

#[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
const CHUNK_BYTES: usize = 16 << 20;

/// Two reused pinned host buffers, their events, and a stream, for chunked copies.
#[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
struct CudaStaging {
    /// Device 0's primary context, retained once for the life of the process.
    context: cudarc::driver::sys::CUcontext,
    stream: cudarc::driver::sys::CUstream,
    buffers: [*mut u8; 2],
    events: [cudarc::driver::sys::CUevent; 2],
}
// SAFETY: the raw handles are only used while holding the mutex in `cuda_staging`.
#[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
unsafe impl Send for CudaStaging {}

/// Locks the shared staging state, creating it on first use. On success the `Option` is `Some`.
/// A failed first initialization leaves it `None`, so a later call retries instead of panicking.
#[cfg(not(all(test, feature = "_cubecl-cpu-test")))]
fn cuda_staging() -> Result<std::sync::MutexGuard<'static, Option<CudaStaging>>, Error> {
    use cudarc::driver::{result, sys};
    static STAGING: std::sync::Mutex<Option<CudaStaging>> = std::sync::Mutex::new(None);
    // Nothing panics while this lock is held except a bug, so poisoning is not expected.
    let mut guard = STAGING.lock().expect("staging lock");
    if let Some(staging) = guard.as_ref() {
        // SAFETY: the context was retained when the staging state was created and is never
        // released; it must be current on whichever thread copies.
        unsafe { result::ctx::set_current(staging.context)? };
    } else {
        // SAFETY: device 0 exists (CubeCL already uses it); its primary context is the one CubeCL
        // uses. It is retained once here, kept in the staging state, and never released. If a
        // later step fails, handles already created leak (a bounded cost, at most once per
        // failed attempt).
        unsafe {
            result::init()?;
            let device = result::device::get(0)?;
            let context = result::primary_ctx::retain(device)?;
            result::ctx::set_current(context)?;
            let stream = result::stream::create(result::stream::StreamKind::NonBlocking)?;
            let buffers = [
                result::malloc_host(CHUNK_BYTES, 0)?.cast::<u8>(),
                result::malloc_host(CHUNK_BYTES, 0)?.cast::<u8>(),
            ];
            let events = [
                result::event::create(sys::CUevent_flags::CU_EVENT_DISABLE_TIMING)?,
                result::event::create(sys::CUevent_flags::CU_EVENT_DISABLE_TIMING)?,
            ];
            *guard = Some(CudaStaging {
                context,
                stream,
                buffers,
                events,
            });
        }
    }
    Ok(guard)
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
    use ndarray as nd;

    use super::{
        standardize::{standardize_beta, standardize_unit, CubeClFloat},
        CubeCl, CubeClVal,
    };
    use crate::{impute_and_zero_mean_snps, sample_bed_file, Bed, BedErrorPlus, Dist, ReadOptions};

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
    ) -> Result<nd::Array2<T>, Box<BedErrorPlus>> {
        let mut options = ReadOptions::<T>::builder();
        options
            .iid_index(iid_index)
            .sid_index(sid_index)
            .is_a1_counted(is_a1_counted)
            .is_f(is_f);
        let context = format!("iid {iid_index:?} sid {sid_index:?} a1 {is_a1_counted} f {is_f}");
        let expected = options.read(bed)?;
        let matrix = options.cubecl().read(bed)?;
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
    ) -> Result<(), Box<BedErrorPlus>> {
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
            let mut matrix = options.cubecl().read(bed)?;
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
    fn gpu_matches_cpu() -> Result<(), Box<BedErrorPlus>> {
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
    fn gpu_builder_chain_order() -> Result<(), Box<BedErrorPlus>> {
        // `.cubecl()` can go anywhere in the chain, and the value-type and order methods work after it.
        let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
        let expected = ReadOptions::builder()
            .i8()
            .c()
            .count_a2()
            .sid_index([2, 0])
            .read(&mut bed)?;
        let first = ReadOptions::builder()
            .cubecl()
            .i8()
            .c()
            .count_a2()
            .sid_index([2, 0])
            .read(&mut bed)?;
        let last = ReadOptions::builder()
            .i8()
            .c()
            .count_a2()
            .sid_index([2, 0])
            .cubecl()
            .read(&mut bed)?;
        let cpu = ReadOptions::builder()
            .i8()
            .cpu()
            .c()
            .count_a2()
            .sid_index([2, 0])
            .read(&mut bed)?;
        assert_eq!(first.to_cpu()?, expected);
        assert_eq!(last.to_cpu()?, expected);
        assert_eq!(cpu, expected);
        assert!(!first.is_f());
        let count_a1 = ReadOptions::builder()
            .f64()
            .cubecl()
            .count_a2()
            .count_a1()
            .f()
            .read(&mut bed)?;
        assert_same(
            &count_a1.to_cpu()?,
            &ReadOptions::builder().f64().read(&mut bed)?,
            "count_a1, f",
        );
        Ok(())
    }

    #[test]
    fn gpu_read_options_match_cpu() -> Result<(), Box<BedErrorPlus>> {
        fn check_options<TVal: CubeClVal + PartialEq>(
            missing_value: TVal,
        ) -> Result<(), Box<BedErrorPlus>> {
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
                        .cubecl()
                        .num_threads(1)
                        .max_concurrent_requests(1)
                        .max_chunk_bytes(1)
                        .read(&mut bed)?;
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
    fn gpu_tensor_handle_describes_the_matrix() -> Result<(), Box<BedErrorPlus>> {
        let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
        for is_f in [true, false] {
            let matrix = ReadOptions::builder()
                .f32()
                .is_f(is_f)
                .cubecl()
                .read(&mut bed)?;
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
    fn gpu_memory_cleanup_runs() -> Result<(), Box<BedErrorPlus>> {
        let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
        let matrix = ReadOptions::builder().f32().cubecl().read(&mut bed)?;
        let expected = matrix.to_cpu()?;
        CubeCl::memory_cleanup()?; // must not free memory that a live matrix uses
        assert_same(&matrix.to_cpu()?, &expected, "after cleanup");
        drop(matrix);
        CubeCl::memory_cleanup()?;
        Ok(())
    }

    #[test]
    fn gpu_standardize_accepts_non_genotype_values() -> Result<(), Box<BedErrorPlus>> {
        // A custom missing value (-1) is observed data, not NaN. It must take
        // the general floating-point path rather than the genotype optimization.
        let mut bed = Bed::new(sample_bed_file("some_missing.bed")?)?;
        let mut matrix = ReadOptions::builder()
            .missing_value(-1.0f32)
            .cubecl()
            .read(&mut bed)?;
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
        let stats = standardize_unit(&mut matrix)?;
        assert!((f64::from(stats[[0, 0]]) - expected_mean).abs() < 4.0 * f64::from(f32::EPSILON));
        assert!(
            (f64::from(stats[[0, 1]]) - expected_variance.sqrt()).abs()
                < 4.0 * f64::from(f32::EPSILON)
        );
        Ok(())
    }

    #[test]
    fn gpu_index_errors_match_cpu() -> Result<(), Box<BedErrorPlus>> {
        let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
        for (iid_index, sid_index) in [
            (vec![3isize], vec![0isize]),
            (vec![-4], vec![0]),
            (vec![0], vec![4]),
            (vec![0], vec![-5]),
        ] {
            let mut options = ReadOptions::<f32>::builder();
            options.iid_index(&iid_index).sid_index(&sid_index);
            let cpu_error = options.read(&mut bed).expect_err("CPU should fail");
            let gpu_error = options
                .cubecl()
                .read(&mut bed)
                .err()
                .expect("GPU should fail");
            assert_eq!(gpu_error.to_string(), cpu_error.to_string());
        }
        Ok(())
    }
}
