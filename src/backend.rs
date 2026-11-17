//! Read genotype data into GPU memory (Cargo feature `cubecl-cuda` or `cubecl-hip`). **Experimental**: this API
//! may change in a minor release.
//!
//! A read's backend is a second type parameter, beside its value type. Just as `.i8()`, `.f32()`, and
//! `.f64()` pick the value type, `.cpu()` (the default) and `.cubecl()` pick the backend, and the two
//! together decide what `read` returns: an `ndarray` array from the CPU, a [`CubeClMatrix`] from the
//! GPU. [`Cpu`] and [`CubeCl`] are the backends; they share the [`Backend`] trait.
//!
//! A GPU read copies the bytes covering the selected individuals of each selected SNP (2 bits per
//! genotype) to the GPU and decodes them there into an `i8`, `f32`, or `f64` [`CubeClMatrix`], in C
//! or F order.
//!
//! # Precision and throughput
//!
//! On many consumer GPUs, `f64` arithmetic is much slower than `f32`; some data-center
//! GPUs have high `f64` throughput. This is a computation tradeoff, not a claim that
//! allocating or decoding `f64` is always slow. `f64` storage also takes twice the
//! memory. Choose `.f32().cubecl()` when its precision is sufficient; keep `f64` when
//! numerical accuracy requires it. For example:
//!
//! ```rust,no_run
//! use bed_reader::{Bed, ReadOptions, sample_bed_file};
//! let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
//! let values = ReadOptions::builder().f32().cubecl().read(&mut bed)?;
//! assert_eq!(values.dim(), (3, 4));
//! # use bed_reader::BedErrorPlus;
//! # Ok::<(), Box<BedErrorPlus>>(())
//! ```
//!
//! Output dtype does not specify every intermediate: standardization may still
//! use `f64` reductions for accuracy. Rust does not emit Python's `GpuFloat64Warning`.
//!
//! The GPU work runs through [CubeCL](https://github.com/tracel-ai/cubecl), so one code base runs on
//! NVIDIA (CUDA) and AMD (HIP) GPUs.
//!
//! # Enabling GPU reads
//!
//! These GPU changes are currently development-only. Use this checkout as a path
//! dependency until a GPU-enabled release is available. From an application
//! directory next to the bed-reader checkout, enable one backend Cargo feature:
//!
//! ```bash
//! # NVIDIA: requires a compatible driver and CUDA toolkit, including nvcc.
//! cargo add bed-reader --path ../bed-reader --features cubecl-cuda
//! # AMD alternative: requires a compatible ROCm/HIP driver and toolkit.
//! cargo add bed-reader --path ../bed-reader --features cubecl-hip
//! ```
//!
//! The internal `cubecl` Cargo feature cannot be enabled alone. CPU reads are still the
//! default; see [`CubeClMatrix`] for a complete typed read, accessor, and host-copy
//! example. Missing values default to -127 for `i8` and NaN for floating-point
//! types; `missing_value` overrides the sentinel on either device. Dtype,
//! `iid_index`, `sid_index`, layout (`f`, `c`, `is_f`), and allele-counting
//! (`count_a1`, `count_a2`, `is_a1_counted`) options all apply to GPU reads.
//! `.cubecl()` is an explicit Rust API choice; Python's `ARRAY_MODULE` setting
//! does not select the Rust read device.
//!
//! Device 0 is used. If CUDA is compiled in it is selected; HIP is selected only
//! when compiled without CUDA. There is no runtime CUDA-to-HIP fallback.
//! GPU reads require local files. Cloud reads and `read_and_fill` remain CPU-only;
//! CPU thread and cloud-request options do not configure the GPU. A decoded
//! matrix has at most `u32::MAX` elements, subject to available memory.
//!
//! [`CubeClMatrix::to_host`] copies data to the CPU. For further GPU computation,
//! use [`CubeClMatrix::tensor_handle`] and [`CubeClMatrix::client`] instead.
//! Standardization helpers are internal Python plumbing, not a public Rust API.
//!
//! # Memory and CUDA interoperation
//!
//! `CubeCL` pools freed allocations. Drop matrices and exported handles before
//! calling [`CubeCl::memory_cleanup`] if another library needs that memory.
//! On CUDA builds, `CubeClMatrix::cuda_pointer` finishes pending `CubeCL` work
//! and pins the allocation. Rust callers must finish foreign CUDA work before
//! releasing the pointer guard or accessing the buffer again through `CubeCL`;
//! the guard does not synchronize foreign streams. Python `CuPy` ownership and
//! synchronization are managed by the bindings.
//!
use core::{fmt::Debug, marker::PhantomData, mem::size_of};
use std::io::{Read, Seek, SeekFrom};

use cubecl::{prelude::*, std::tensor::TensorHandle};
use ndarray as nd;
use ndarray::ShapeBuilder;

use crate::{
    open_and_check, path_ref_to_string, Bed, BedError, BedErrorPlus, BedVal, Hold, ReadOptions,
    ReadOptionsBuilder, CB_HEADER_U64,
};

// The standardizer is for the Python bindings (and this module's tests), like the CPU one.
#[cfg(any(feature = "extension-module", test))]
pub(crate) mod standardize;

/// Units per cube.
const UNITS: u32 = 256;
/// Cubes along a grid's y axis are limited to this many; kernels stride over the rest.
const MAX_GRID_Y: usize = 65_535;

/// GPU errors.
///
/// This enum is non-exhaustive so additional backend and runtime failures can be
/// represented in minor releases. Downstream matches must include a fallback arm.
#[derive(thiserror::Error, Debug)]
#[non_exhaustive]
pub enum Error {
    /// No usable GPU, or this build does not include its runtime.
    #[error(transparent)]
    Unavailable(#[from] cubecl::device::DeviceUnavailable),

    /// A GPU runtime error.
    #[error(transparent)]
    Server(#[from] cubecl_runtime::server::ServerError),

    /// Legacy genotype-only error, retained for compatibility. General floating-point
    /// standardization now handles columns outside the integer optimization.
    #[error("GPU standardization needs values 0, 1, 2, or NaN")]
    NotGenotype,

    /// Centered floating-point statistics overflowed or were otherwise invalid.
    #[error("Invalid GPU variance: {0}")]
    InvalidVariance(f64),

    /// The GPU code currently indexes with `u32`.
    #[error("Too large for the GPU code ({0} exceeds u32::MAX)")]
    TooLarge(usize),
}

impl From<Error> for Box<BedErrorPlus> {
    fn from(error: Error) -> Self {
        Box::new(BedErrorPlus::GpuError(error))
    }
}

mod sealed {
    /// Keeps [`Backend`](super::Backend) to the backends this crate defines.
    pub trait Sealed {}

    /// Restricts `CubeCL` element types to the types supported by this crate.
    pub trait Value {}

    impl Value for i8 {}
    impl Value for f32 {}
    impl Value for f64 {}
}

/// Where a read puts its result: [`Cpu`] or [`CubeCl`].
///
/// It is the second type parameter of [`ReadOptions`] and [`ReadOptionsBuilder`], as
/// [`BedVal`] (`i8`, `f32`, `f64`) is the first. Together they decide the type `read` returns.
///
/// This trait is sealed: downstream crates cannot implement it. bed-reader can
/// add backend types, builder methods, and opt-in Cargo features in minor releases
/// while preserving existing backends' entry points and result types.
///
/// # Example
///
/// See [`CubeClMatrix`].
pub trait Backend: sealed::Sealed + Clone + Copy + Debug + Default + Send + Sync + 'static {
    /// What `read` returns for value type `TVal` with this backend.
    type Matrix<TVal: BedVal>;
}

/// The CPU, the default [`Backend`]: `read` returns an `ndarray` array.
///
/// # Example
///
/// See [`ReadOptionsBuilder::cpu`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Cpu;

/// GPU device 0, through this build's selected `CubeCL` runtime: `read` returns a
/// [`CubeClMatrix`]. CUDA takes compile-time precedence if both vendor Cargo features
/// are enabled; there is no runtime fallback to HIP.
///
/// # Example
///
/// See [`CubeClMatrix`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CubeCl;

impl sealed::Sealed for Cpu {}
impl sealed::Sealed for CubeCl {}

impl Backend for Cpu {
    type Matrix<TVal: BedVal> = nd::Array2<TVal>;
}

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
        #[cfg(feature = "cubecl-cuda")]
        let result = cubecl::Device::cuda(0);
        #[cfg(all(feature = "cubecl-hip", not(feature = "cubecl-cuda")))]
        let result = cubecl::Device::rocm(0);
        let device = result.map_err(Error::from)?;
        Ok(device.client())
    }

    /// Ask `CubeCL` to return GPU memory that it holds but no longer uses.
    ///
    /// `CubeCL` keeps the memory of dropped matrices in a pool for its own later use. Other GPU
    /// libraries in the process (for example, `CuPy`) cannot use that memory until it is returned.
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
    /// let matrix = ReadOptions::builder().f32().cubecl().read(&mut bed)?;
    /// drop(matrix);
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
pub trait CubeClVal: sealed::Value + BedVal + Numeric + CubeElement + bytemuck::Pod {}
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
/// let matrix = ReadOptions::builder().f32().f().cubecl().read(&mut bed)?;
/// assert_eq!(matrix.dim(), (3, 4));
/// assert!(matrix.is_f());
/// let tensor_handle = matrix.tensor_handle(); // CubeCL's own type, for use with ...
/// let client = matrix.client(); // ... its client
/// assert_eq!(tensor_handle.shape().as_slice(), &[3, 4]);
/// println!("{}", client.name());
/// let val = matrix.to_host()?; // explicit CPU copy: ndarray::Array2<f32>
/// assert_eq!(val.dim(), (3, 4));
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

impl<TVal: BedVal> ReadOptionsBuilder<TVal> {
    /// Read on the CPU, into an `ndarray` array. This is the default, so it changes nothing; it
    /// is the counterpart of [`ReadOptionsBuilder::cubecl`].
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    ///
    /// let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
    /// let val = ReadOptions::builder().f32().cpu().read(&mut bed)?;
    /// assert_eq!(val.dim(), (3, 4));
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn cpu(&mut self) -> &mut Self {
        self
    }

    /// Read on the GPU, through `CubeCL`, into a [`CubeClMatrix`].
    ///
    /// Unlike `.f32()`, this changes the builder's type (the CPU is where every builder starts),
    /// so it returns a new builder. It can go anywhere in the chain.
    ///
    /// # Example
    ///
    /// See [`CubeClMatrix`].
    #[must_use]
    pub fn cubecl(&self) -> ReadOptionsBuilder<TVal, CubeCl> {
        ReadOptionsBuilder {
            missing_value: self.missing_value,
            iid_index: self.iid_index.clone(),
            sid_index: self.sid_index.clone(),
            is_f: self.is_f,
            is_a1_counted: self.is_a1_counted,
            num_threads: self.num_threads,
            max_concurrent_requests: self.max_concurrent_requests,
            max_chunk_bytes: self.max_chunk_bytes,
            _backend: PhantomData,
        }
    }
}

impl<TVal: CubeClVal> ReadOptionsBuilder<TVal, CubeCl> {
    /// Read genotype data into GPU memory. The options that set CPU threads and cloud chunking
    /// do not apply.
    ///
    /// # Errors
    /// See [`BedError`] and [`BedErrorPlus`] (including its `GpuError` variant).
    ///
    /// # Example
    ///
    /// See [`CubeClMatrix`].
    pub fn read(&self, bed: &mut Bed) -> Result<CubeClMatrix<TVal>, Box<BedErrorPlus>> {
        let read_options = self.build()?;
        read_cubecl(bed, &read_options)
    }
}

impl<TVal: BedVal> ReadOptionsBuilder<TVal, CubeCl> {
    /// Order of the output matrix, Fortran-style (default). As [`ReadOptionsBuilder::f`] on the CPU.
    ///
    /// # Example
    ///
    /// See [`CubeClMatrix`].
    pub fn f(&mut self) -> &mut Self {
        self.is_f(true);
        self
    }

    /// Order of the output matrix, C. As [`ReadOptionsBuilder::c`] on the CPU.
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    ///
    /// let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
    /// let matrix = ReadOptions::builder().f64().cubecl().c().count_a2().read(&mut bed)?;
    /// assert!(!matrix.is_f());
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn c(&mut self) -> &mut Self {
        self.is_f(false);
        self
    }

    /// Count the number of A1 alleles (default). As [`ReadOptionsBuilder::count_a1`] on the CPU.
    ///
    /// # Example
    ///
    /// See [`ReadOptionsBuilder::c`](#method.c-1), which uses the similar `count_a2`.
    pub fn count_a1(&mut self) -> &mut Self {
        self.is_a1_counted = Some(true);
        self
    }

    /// Count the number of A2 alleles. As [`ReadOptionsBuilder::count_a2`] on the CPU.
    ///
    /// # Example
    ///
    /// See [`ReadOptionsBuilder::c`](#method.c-1).
    pub fn count_a2(&mut self) -> &mut Self {
        self.is_a1_counted = Some(false);
        self
    }
}

impl ReadOptionsBuilder<i8, CubeCl> {
    /// Output a [`CubeClMatrix`] of `i8`. As [`ReadOptionsBuilder::i8`] on the CPU.
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    ///
    /// let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
    /// let matrix = ReadOptions::builder().cubecl().i8().read(&mut bed)?;
    /// assert_eq!(matrix.to_host()?[[0, 0]], 1i8);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn i8(&mut self) -> &mut Self {
        self
    }
}

impl ReadOptionsBuilder<f32, CubeCl> {
    /// Output a [`CubeClMatrix`] of `f32`. As [`ReadOptionsBuilder::f32`] on the CPU.
    ///
    /// # Example
    ///
    /// See [`ReadOptionsBuilder::i8`](#method.i8-1), which is similar.
    pub fn f32(&mut self) -> &mut Self {
        self
    }
}

impl ReadOptionsBuilder<f64, CubeCl> {
    /// Output a [`CubeClMatrix`] of `f64`. As [`ReadOptionsBuilder::f64`] on the CPU.
    /// On many consumer GPUs, `f64` arithmetic is substantially slower than `f32`.
    /// Consider [`ReadOptionsBuilder::f32`] if its precision is sufficient. Storage
    /// also takes twice the memory; this method never silently changes precision.
    ///
    /// # Example
    ///
    /// See [`ReadOptionsBuilder::i8`](#method.i8-1), which is similar.
    pub fn f64(&mut self) -> &mut Self {
        self
    }
}

/// Reads into GPU memory: the packed bytes covering the selection go to the GPU, which decodes them.
fn read_cubecl<TVal: CubeClVal>(
    bed: &mut Bed,
    read_options: &ReadOptions<TVal, CubeCl>,
) -> Result<CubeClMatrix<TVal>, Box<BedErrorPlus>> {
    let iid_count = bed.iid_count()?;
    let sid_count = bed.sid_count()?;
    let iid_hold = Hold::new(&read_options.iid_index, iid_count)?;
    let sid_hold = Hold::new(&read_options.sid_index, sid_count)?;
    let iid_index = resolve(iid_hold.as_ref(), iid_count, BedError::IidIndexTooBig)?;
    let sid_index = resolve(sid_hold.as_ref(), sid_count, BedError::SidIndexTooBig)?;
    let element_count = iid_index
        .len()
        .checked_mul(sid_index.len())
        .ok_or(Error::TooLarge(usize::MAX))?;
    to_u32(element_count)?;
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
    let handle = client.empty((element_count.max(1) * size_of::<TVal>()).max(4));
    if element_count > 0 {
        let packed_handle = upload_u32(&client, &packed)?;
        let col_word_handle = client.create_from_slice(u32::as_bytes(&col_word));
        let col_shift_handle = client.create_from_slice(u32::as_bytes(&col_shift));
        let lut_handle = genotype_lookup(&client, read_options);
        let (fast_count, slow_count) = memory_dimensions(cols.len(), rows.len(), kernel_is_f);
        let element_count_u32 = to_u32(element_count)?;
        // SAFETY: `decode` writes only `slow * fast_count + fast` for `fast < fast_count` and
        // `slow < slow_count`, so below `element_count` (which fits u32), the length of
        // `handle`; it reads `col_word`/`col_shift` at column indexes < `cols.len()` and
        // `packed` at word indexes < `rows.len() * span_words`, the lengths passed.
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
    pub(crate) fn empty(rows: usize, columns: usize, is_f: bool) -> Result<Self, Box<BedErrorPlus>> {
        to_u32(rows)?;
        to_u32(columns)?;
        let element_count = rows
            .checked_mul(columns)
            .ok_or(Error::TooLarge(usize::MAX))?;
        to_u32(element_count)?;
        let byte_count = element_count
            .checked_mul(size_of::<TVal>())
            .ok_or(Error::TooLarge(usize::MAX))?;
        let client = CubeCl::client()?;
        // Keep a real pointer and one valid dummy element for empty kernel arguments.
        let handle = client.empty(byte_count.max(size_of::<TVal>()).max(4));
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

    /// Copy the matrix to host memory, in its memory order.
    ///
    /// See the [`CubeClMatrix`] example.
    ///
    /// # Errors
    /// [`Error::Server`] if the GPU fails.
    pub fn to_host(&self) -> Result<nd::Array2<TVal>, Box<BedErrorPlus>> {
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
    /// let matrix = ReadOptions::builder().f32().cubecl().read(&mut bed)?;
    /// let cuda_pointer = matrix.cuda_pointer()?;
    /// assert_ne!(cuda_pointer.ptr(), 0);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[cfg(feature = "cubecl-cuda")]
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
#[cfg(feature = "cubecl-cuda")]
pub struct CudaPointer {
    ptr: u64,
    /// `CubeCL's` `ManagedResource`, which pins the allocation in place.
    _resource: Box<dyn core::any::Any + Send>,
}

#[cfg(feature = "cubecl-cuda")]
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
    u32::try_from(value).map_err(|_| Error::TooLarge(value))
}

fn genotype<TVal: CubeClVal>(value: i8) -> TVal {
    <TVal as From<i8>>::from(value)
}

/// Copies `words` to a new GPU buffer (through reused pinned host buffers on CUDA, which is about 4x
/// faster than `CubeCL's` own upload under `WSL2`).
fn upload_u32(client: &Client, words: &[u32]) -> Result<cubecl::server::Handle, Box<BedErrorPlus>> {
    #[cfg(feature = "cubecl-cuda")]
    {
        let handle = client.empty(core::mem::size_of_val(words));
        upload_pinned_cuda(client, &handle, bytemuck::cast_slice(words))?;
        Ok(handle)
    }
    #[cfg(not(feature = "cubecl-cuda"))]
    {
        Ok(client.create_from_slice(u32::as_bytes(words)))
    }
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
#[cfg(feature = "cubecl-cuda")]
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
    let staging = cuda_staging();
    for (chunk_index, chunk) in bytes.chunks(CHUNK_BYTES).enumerate() {
        let slot = chunk_index % 2;
        // SAFETY: the event guards the pinned buffer until its previous copy finished; the
        // destination range lies inside the CubeCL buffer, which `resource` keeps alive and in place.
        unsafe {
            result::event::synchronize(staging.events[slot]).expect("CUDA event sync");
            let buffer = core::slice::from_raw_parts_mut(staging.buffers[slot], chunk.len());
            buffer.copy_from_slice(chunk);
            result::memcpy_htod_async(
                device_ptr + (chunk_index * CHUNK_BYTES) as u64,
                buffer,
                staging.stream,
            )
            .expect("CUDA copy to GPU");
            result::event::record(staging.events[slot], staging.stream).expect("CUDA event record");
        }
    }
    // SAFETY: the stream is valid; waiting makes the data visible to CubeCL's later launches.
    unsafe { result::stream::synchronize(staging.stream).expect("CUDA stream sync") };
    drop(resource);
    Ok(())
}

#[cfg(feature = "cubecl-cuda")]
const CHUNK_BYTES: usize = 16 << 20;

/// Two reused pinned host buffers, their events, and a stream, for chunked copies.
#[cfg(feature = "cubecl-cuda")]
struct CudaStaging {
    stream: cudarc::driver::sys::CUstream,
    buffers: [*mut u8; 2],
    events: [cudarc::driver::sys::CUevent; 2],
}
// SAFETY: the raw handles are only used while holding the mutex in `cuda_staging`.
#[cfg(feature = "cubecl-cuda")]
unsafe impl Send for CudaStaging {}

#[cfg(feature = "cubecl-cuda")]
fn cuda_staging() -> std::sync::MutexGuard<'static, CudaStaging> {
    use cudarc::driver::{result, sys};
    static STAGING: std::sync::OnceLock<std::sync::Mutex<CudaStaging>> = std::sync::OnceLock::new();
    let staging = STAGING.get_or_init(|| {
        // SAFETY: device 0 exists (CubeCL already uses it); the primary context is the one CubeCL
        // uses.
        unsafe {
            result::init().expect("cuInit");
            let device = result::device::get(0).expect("CUDA device 0");
            let context = result::primary_ctx::retain(device).expect("CUDA primary context");
            result::ctx::set_current(context).expect("CUDA set context");
            let stream = result::stream::create(result::stream::StreamKind::NonBlocking)
                .expect("CUDA stream");
            let buffers = [0, 1].map(|_| {
                result::malloc_host(CHUNK_BYTES, 0)
                    .expect("pinned buffer")
                    .cast::<u8>()
            });
            let events = [0, 1].map(|_| {
                result::event::create(sys::CUevent_flags::CU_EVENT_DISABLE_TIMING)
                    .expect("CUDA event")
            });
            std::sync::Mutex::new(CudaStaging {
                stream,
                buffers,
                events,
            })
        }
    });
    let guard = staging.lock().expect("staging lock");
    // SAFETY: device 0's primary context is valid for the life of the process; it must be current
    // on whichever thread copies.
    unsafe {
        let device = result::device::get(0).expect("CUDA device 0");
        let context = result::primary_ctx::retain(device).expect("CUDA primary context");
        result::ctx::set_current(context).expect("CUDA set context");
    }
    guard
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
        assert_same(&matrix.to_host()?, &expected, &context);
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
                        &matrix.to_host()?,
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
        assert_eq!(first.to_host()?, expected);
        assert_eq!(last.to_host()?, expected);
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
            &count_a1.to_host()?,
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
                    let actual = matrix.to_host()?;
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
            let host = matrix.to_host()?;
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
        let expected = matrix.to_host()?;
        CubeCl::memory_cleanup()?; // must not free memory that a live matrix uses
        assert_same(&matrix.to_host()?, &expected, "after cleanup");
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
        let original = matrix.to_host()?;
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
