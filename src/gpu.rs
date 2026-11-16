//! Read genotype data into GPU memory and standardize it there (`gpu-cuda` or `gpu-hip` feature).
//!
//! [`Bed::read_gpu_with_options`] (or [`ReadOptionsBuilder::read_gpu`]) reads the bytes covering the
//! selected individuals of each selected SNP, copies those packed bytes (2 bits per genotype) to the
//! GPU, and decodes them there into an `f32`/`f64` [`Matrix`] in C or F order. [`Matrix::standardize_unit`]
//! and [`Matrix::standardize_beta`] standardize it in place, matching the CPU standardization
//! exactly: the GPU computes exact integer sums per SNP, the host finishes each SNP's statistics with
//! the same function the CPU code uses, and the GPU applies them.
//!
//! The GPU work runs through [CubeCL](https://github.com/tracel-ai/cubecl), so one code base runs on
//! NVIDIA (CUDA) and AMD (HIP) GPUs.
use core::marker::PhantomData;
use std::io::{Read, Seek, SeekFrom};

use cubecl::prelude::*;
use ndarray as nd;
use ndarray::ShapeBuilder;

use crate::{
    find_factor, open_and_check, path_ref_to_string, snp_mean_std, Bed, BedError, BedErrorPlus,
    BedVal, Dist, Hold, ReadOptions, ReadOptionsBuilder, CB_HEADER_U64,
};

/// Units per cube.
const UNITS: u32 = 256;
/// Cubes along a grid's y axis are limited to this many; kernels stride over the rest.
const MAX_GRID_Y: usize = 65_535;

/// GPU errors.
#[derive(thiserror::Error, Debug)]
#[non_exhaustive]
pub enum Error {
    /// No usable GPU, or this build does not include its runtime.
    #[error(transparent)]
    Unavailable(#[from] cubecl::device::DeviceUnavailable),

    /// A GPU runtime error.
    #[error(transparent)]
    Server(#[from] cubecl_runtime::server::ServerError),

    /// GPU standardization needs genotype values (0, 1, 2, or missing), so its sums are exact.
    #[error("GPU standardization needs values 0, 1, 2, or NaN")]
    NotGenotype,

    /// The GPU code currently indexes with `u32`.
    #[error("Too large for the GPU code ({0} exceeds u32::MAX)")]
    TooLarge(usize),
}

impl From<Error> for Box<BedErrorPlus> {
    fn from(error: Error) -> Self {
        Box::new(BedErrorPlus::GpuError(error))
    }
}

/// A GPU to read into.
///
/// # Example
///
/// See [`Matrix`].
pub struct Device {
    client: Client,
}

impl Device {
    /// The first GPU of a runtime this build includes (CUDA, then HIP).
    ///
    /// # Errors
    /// [`Error::Unavailable`] if there is none.
    pub fn new() -> Result<Self, Box<BedErrorPlus>> {
        #[cfg(feature = "gpu-cuda")]
        let result = cubecl::Device::cuda(0);
        #[cfg(all(feature = "gpu-hip", not(feature = "gpu-cuda")))]
        let result = cubecl::Device::rocm(0);
        let device = result.map_err(Error::from)?;
        Ok(Self { client: device.client() })
    }
}

/// A float type the GPU code reads into: `f32` or `f64`.
pub trait GpuVal:
    BedVal
    + Float
    + CubeElement
    + bytemuck::Pod
    + num_traits::Float
    + num_traits::ToPrimitive
    + num_traits::FromPrimitive
{
}
impl GpuVal for f32 {}
impl GpuVal for f64 {}

/// A genotype matrix (individuals x SNPs) in GPU memory.
///
/// # Example
///
/// ```no_run
/// use bed_reader::{gpu, Bed, ReadOptions, sample_bed_file};
///
/// let device = gpu::Device::new()?;
/// let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
/// let mut matrix = ReadOptions::builder().f32().f().read_gpu(&mut bed, &device)?;
/// let stats = matrix.standardize_unit()?; // in place; each SNP's (mean, std)
/// let val = matrix.to_host()?;
/// assert_eq!(val.dim(), (3, 4));
/// assert_eq!(stats.dim(), (4, 2));
/// # use bed_reader::BedErrorPlus;
/// # Ok::<(), Box<BedErrorPlus>>(())
/// ```
pub struct Matrix<TVal: GpuVal> {
    client: Client,
    handle: cubecl::server::Handle,
    iid_count: usize,
    sid_count: usize,
    is_f: bool,
    phantom: PhantomData<TVal>,
}

impl Bed {
    /// Read genotype data into GPU memory with options. Only floats (`f32`, `f64`) are supported.
    ///
    /// > Also see [`ReadOptionsBuilder::read_gpu`], which this backs, and the [`Matrix`] example.
    ///
    /// # Errors
    /// See [`BedError`] and [`BedErrorPlus`] (including its `GpuError` variant).
    pub fn read_gpu_with_options<TVal: GpuVal>(
        &mut self,
        read_options: &ReadOptions<TVal>,
        device: &Device,
    ) -> Result<Matrix<TVal>, Box<BedErrorPlus>> {
        let iid_count = self.iid_count()?;
        let sid_count = self.sid_count()?;
        let iid_hold = Hold::new(&read_options.iid_index, iid_count)?;
        let sid_hold = Hold::new(&read_options.sid_index, sid_count)?;
        let iid_index = resolve(iid_hold.as_ref(), iid_count, BedError::IidIndexTooBig)?;
        let sid_index = resolve(sid_hold.as_ref(), sid_count, BedError::SidIndexTooBig)?;
        let client = device.client.clone();

        let (buf_reader, header) = open_and_check(&self.path)?;
        let mut file = buf_reader.into_inner();
        // Mode 1 stores each SNP's individuals together (SNP-major); mode 0 stores each individual's
        // SNPs together. For mode 0, swap the roles of individuals and SNPs and flip the memory
        // order: the transposed result has the requested layout (as the CPU reader does).
        let (row_count, rows, col_count, cols, kernel_is_f) = match header[2] {
            1 => (sid_count, &sid_index, iid_count, &iid_index, read_options.is_f),
            0 => (iid_count, &iid_index, sid_count, &sid_index, !read_options.is_f),
            _ => Err(BedError::BadMode(path_ref_to_string(&self.path)))?,
        };
        let row_bytes = col_count.div_ceil(4);
        let expected_len = (row_bytes as u64)
            .checked_mul(row_count as u64)
            .and_then(|len| len.checked_add(CB_HEADER_U64))
            .ok_or(BedError::IndexesTooBigForFiles(iid_count, sid_count))?;
        if file.metadata()?.len() != expected_len {
            Err(BedError::IllFormed(path_ref_to_string(&self.path)))?;
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

        let element_count = cols.len() * rows.len();
        let handle = client.empty((element_count * core::mem::size_of::<TVal>()).max(4));
        if element_count > 0 {
            let packed_handle = upload_u32(&client, &packed)?;
            let col_word_handle = client.create_from_slice(u32::as_bytes(&col_word));
            let col_shift_handle = client.create_from_slice(u32::as_bytes(&col_shift));
            let missing = read_options.missing_value;
            let lut: [TVal; 4] = if read_options.is_a1_counted {
                [genotype(2), missing, genotype(1), genotype(0)]
            } else {
                [genotype(0), missing, genotype(1), genotype(2)]
            };
            let lut_handle = client.create_from_slice(bytemuck::cast_slice(&lut));
            let (fast_count, slow_count) = if kernel_is_f {
                (cols.len(), rows.len())
            } else {
                (rows.len(), cols.len())
            };
            let element_count_u32 = to_u32(element_count)?;
            // SAFETY: `decode` writes only `slow * fast_count + fast` for `fast < fast_count` and
            // `slow < slow_count`, so below `element_count` (which fits u32), the length of
            // `handle`; it reads `col_word`/`col_shift` at column indexes < `cols.len()` and
            // `packed` at word indexes < `rows.len() * span_words`, the lengths passed.
            unsafe {
                decode::launch_unchecked::<TVal>(
                    &client,
                    CubeCount::Static(to_u32(fast_count)?.div_ceil(UNITS), to_u32(slow_count.min(MAX_GRID_Y))?, 1),
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
        Ok(Matrix {
            client,
            handle,
            iid_count: iid_index.len(),
            sid_count: sid_index.len(),
            is_f: read_options.is_f,
            phantom: PhantomData,
        })
    }
}

impl<TVal: GpuVal> ReadOptionsBuilder<TVal> {
    /// Read genotype data into GPU memory. Only floats (`f32`, `f64`) are supported.
    ///
    /// > Also see [`Bed::read_gpu_with_options`] and the [`Matrix`] example.
    ///
    /// # Errors
    /// See [`BedError`] and [`BedErrorPlus`] (including its `GpuError` variant).
    pub fn read_gpu(&self, bed: &mut Bed, device: &Device) -> Result<Matrix<TVal>, Box<BedErrorPlus>> {
        let read_options = self.build()?;
        bed.read_gpu_with_options(&read_options, device)
    }
}

impl<TVal: GpuVal> Matrix<TVal> {
    /// The number of individuals and SNPs.
    ///
    /// See the [`Matrix`] example.
    #[must_use]
    pub fn dim(&self) -> (usize, usize) {
        (self.iid_count, self.sid_count)
    }

    /// True for F (column-major) order, false for C (row-major) order.
    ///
    /// See the [`Matrix`] example.
    #[must_use]
    pub fn is_f(&self) -> bool {
        self.is_f
    }

    /// Copy the matrix to host memory, in its memory order.
    ///
    /// See the [`Matrix`] example.
    ///
    /// # Errors
    /// [`Error::Server`] if the GPU fails.
    pub fn to_host(&self) -> Result<nd::Array2<TVal>, Box<BedErrorPlus>> {
        let element_count = self.iid_count * self.sid_count;
        let shape = (self.iid_count, self.sid_count).set_f(self.is_f);
        if element_count == 0 {
            return Ok(nd::Array2::default(shape));
        }
        let bytes = self.client.read_one(self.handle.clone()).map_err(Error::from)?;
        let values = bytemuck::cast_slice::<u8, TVal>(&bytes[..element_count * core::mem::size_of::<TVal>()]);
        Ok(nd::Array2::from_shape_vec(shape, values.to_vec()).expect("shape matches element count"))
    }

    /// Standardize in place to zero mean and unit variance (missing values become 0), matching the
    /// CPU standardization exactly. Returns each SNP's mean and standard deviation (`sid_count x 2`).
    ///
    /// See the [`Matrix`] example.
    ///
    /// # Errors
    /// [`BedError::NoIndividuals`] for a SNP with no observed values; [`Error::NotGenotype`] if a
    /// value is not 0, 1, 2, or NaN.
    pub fn standardize_unit(&mut self) -> Result<nd::Array2<TVal>, Box<BedErrorPlus>> {
        self.standardize(&Dist::Unit)
    }

    /// Standardize in place with Beta(`a`, `b`) weights (missing values become 0), matching the CPU
    /// standardization exactly. Returns each SNP's mean and standard deviation (`sid_count x 2`).
    ///
    /// Also see the [`Matrix`] example.
    ///
    /// ```no_run
    /// use bed_reader::{gpu, Bed, ReadOptions, sample_bed_file};
    ///
    /// let device = gpu::Device::new()?;
    /// let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
    /// let mut matrix = ReadOptions::builder().f64().read_gpu(&mut bed, &device)?;
    /// let stats = matrix.standardize_beta(1.0, 25.0)?;
    /// assert_eq!(stats.dim(), (4, 2));
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    ///
    /// # Errors
    /// As [`Matrix::standardize_unit`], plus [`BedError::IllegalSnpMean`] and
    /// [`BedError::CannotCreateBetaDist`].
    pub fn standardize_beta(&mut self, a: f64, b: f64) -> Result<nd::Array2<TVal>, Box<BedErrorPlus>> {
        self.standardize(&Dist::Beta { a, b })
    }

    fn standardize(&mut self, dist: &Dist) -> Result<nd::Array2<TVal>, Box<BedErrorPlus>> {
        let (iid_count, sid_count) = (self.iid_count, self.sid_count);
        let mut stats = nd::Array2::<TVal>::zeros((sid_count, 2));
        if sid_count == 0 {
            return Ok(stats);
        }
        let element_count = iid_count * sid_count;
        let c_order = u32::from(!self.is_f);
        let client = &self.client;
        let sums_handle = client.empty(4 * sid_count * 4);
        // SAFETY: one cube per SNP (`sid_j < sid_count`); it reads the matrix only at indexes below
        // `iid_count * sid_count` (= `element_count`, which fits u32) and writes
        // `sums[4 * sid_j .. 4 * sid_j + 4]`, inside the `4 * sid_count` words allocated.
        unsafe {
            snp_sums::launch_unchecked::<TVal>(
                client,
                CubeCount::Static(to_u32(sid_count)?, 1, 1),
                CubeDim::new_1d(UNITS),
                BufferArg::from_raw_parts(self.handle.clone(), to_u32(element_count.max(1))? as usize),
                to_u32(iid_count)?,
                to_u32(sid_count)?,
                c_order,
                BufferArg::from_raw_parts(sums_handle.clone(), 4 * sid_count),
            );
        }
        let sums_bytes = client.read_one(sums_handle).map_err(Error::from)?;
        let sums: &[u32] = bytemuck::cast_slice(&sums_bytes[..4 * sid_count * 4]);

        // Finish each SNP's statistics with the CPU code's functions. The integer sums convert to
        // `TVal` exactly whenever the CPU's own float sums are exact (below 2^24 for f32).
        let two = genotype::<TVal>(2);
        let mut mean = Vec::with_capacity(sid_count);
        let mut factor = Vec::with_capacity(sid_count);
        let mut snc = Vec::with_capacity(sid_count);
        for (&[count, sum, sum2, not_genotype], mut stats_row) in sums.as_chunks::<4>().0.iter().zip(stats.axis_iter_mut(nd::Axis(0))) {
            if not_genotype != 0 {
                Err(Error::NotGenotype)?;
            }
            let (mean_s, std) = snp_mean_std(from_u32::<TVal>(count), from_u32(sum), from_u32(sum2), dist, two)?;
            stats_row[0] = mean_s;
            stats_row[1] = std;
            mean.push(mean_s);
            factor.push(find_factor(dist, mean_s, std)?);
            snc.push(u32::from(num_traits::Float::is_infinite(std)));
        }

        if element_count > 0 {
            let mean_handle = client.create_from_slice(bytemuck::cast_slice(&mean));
            let factor_handle = client.create_from_slice(bytemuck::cast_slice(&factor));
            let snc_handle = client.create_from_slice(u32::as_bytes(&snc));
            let (fast_count, slow_count) = if self.is_f { (iid_count, sid_count) } else { (sid_count, iid_count) };
            // SAFETY: as for `decode`: writes only `slow * fast_count + fast` < `element_count`;
            // reads the per-SNP arrays at `sid_j < sid_count`, their length.
            unsafe {
                apply::launch_unchecked::<TVal>(
                    client,
                    CubeCount::Static(to_u32(fast_count)?.div_ceil(UNITS), to_u32(slow_count.min(MAX_GRID_Y))?, 1),
                    CubeDim::new_1d(UNITS),
                    BufferArg::from_raw_parts(self.handle.clone(), element_count),
                    to_u32(iid_count)?,
                    to_u32(sid_count)?,
                    c_order,
                    BufferArg::from_raw_parts(mean_handle, sid_count),
                    BufferArg::from_raw_parts(factor_handle, sid_count),
                    BufferArg::from_raw_parts(snc_handle, sid_count),
                );
            }
        }
        Ok(stats)
    }

    /// CUDA only: the matrix's device pointer, for handing it to another CUDA library (for
    /// example, `CuPy`) without a copy. Finishes this matrix's pending GPU work first.
    ///
    /// The returned [`CudaPointer`] keeps the allocation alive and in place; use the pointer only
    /// while both it and this matrix exist.
    ///
    /// # Errors
    /// [`Error::Server`] if the GPU fails.
    #[cfg(feature = "gpu-cuda")]
    pub fn cuda_pointer(&self) -> Result<CudaPointer, Box<BedErrorPlus>> {
        use cubecl::cuda::CudaRuntime;
        use cubecl_runtime::runtime::Runtime;
        pollster::block_on(self.client.sync()).map_err(Error::from)?;
        let resource = self
            .client
            .get_resource::<<CudaRuntime as Runtime>::Server>(self.handle.clone())
            .map_err(Error::from)?;
        let ptr = resource.resource().ptr;
        Ok(CudaPointer { ptr, _resource: Box::new(resource) })
    }
}

/// CUDA only: a matrix's device pointer, from [`Matrix::cuda_pointer`]. While it exists, `CubeCL`
/// keeps the allocation alive and never relocates it.
#[cfg(feature = "gpu-cuda")]
pub struct CudaPointer {
    ptr: u64,
    /// `CubeCL's` `ManagedResource`, which pins the allocation in place.
    _resource: Box<dyn core::any::Any + Send>,
}

#[cfg(feature = "gpu-cuda")]
impl CudaPointer {
    /// The device pointer.
    ///
    /// See [`Matrix::cuda_pointer`].
    #[must_use]
    pub fn ptr(&self) -> u64 {
        self.ptr
    }
}

/// Resolves Python-style indexes (negative counts from the end), with the CPU reader's errors.
fn resolve(index: &[isize], count: usize, too_big: fn(isize) -> BedError) -> Result<Vec<usize>, BedError> {
    let count_isize = isize::try_from(count).map_err(|_| BedError::IndexesTooBigForFiles(count, 0))?;
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

fn genotype<TVal: GpuVal>(value: i8) -> TVal {
    <TVal as From<i8>>::from(value)
}

fn from_u32<TVal: GpuVal>(value: u32) -> TVal {
    // f32/f64 from u32 is the CPU's sum of the same integers whenever that sum is exact.
    <TVal as num_traits::NumCast>::from(value).expect("u32 converts to f32/f64")
}

/// Copies `words` to a new GPU buffer (through reused pinned host buffers on CUDA, which is about 4x
/// faster than `CubeCL's` own upload under `WSL2`).
fn upload_u32(client: &Client, words: &[u32]) -> Result<cubecl::server::Handle, Box<BedErrorPlus>> {
    #[cfg(feature = "gpu-cuda")]
    {
        let handle = client.empty(core::mem::size_of_val(words));
        upload_pinned_cuda(client, &handle, bytemuck::cast_slice(words))?;
        Ok(handle)
    }
    #[cfg(not(feature = "gpu-cuda"))]
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
fn decode<F: Float>(
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
            slow += CUBE_COUNT_Y;
        }
    }
}

/// Per-SNP sums over observed (non-NaN) values: one cube per SNP. Writes `count`, `sum`,
/// `sum of squares`, and a flag set if any value is not 0, 1, or 2 to `sums[4 * sid_j ..]`. The
/// values are small integers, so these integer sums are exact in any order.
#[cube(launch_unchecked)]
fn snp_sums<F: Float>(val: &[F], iid_count: u32, sid_count: u32, c_order: u32, sums: &mut [u32]) {
    let sid_j = CUBE_POS_X;
    let mut count = 0u32;
    let mut sum = 0u32;
    let mut sum2 = 0u32;
    let mut not_genotype = 0u32;
    if sid_j < sid_count {
        let mut iid_i = UNIT_POS_X;
        while iid_i < iid_count {
            let mut index = sid_j * iid_count + iid_i;
            if c_order == 1u32 {
                index = iid_i * sid_count + sid_j;
            }
            let value = val[index as usize];
            // NaN (missing) fails every ordered comparison, so it takes neither branch and is not
            // counted. (`F::is_nan` returns a line-wide bool, which a scalar `if` can't use.)
            if value >= F::from_int(0) {
                let genotype = u32::cast_from(value);
                if genotype > 2u32 || F::cast_from(genotype) != value {
                    not_genotype = 1u32;
                }
                count += 1u32;
                sum += genotype;
                sum2 += genotype * genotype;
            } else if value < F::from_int(0) {
                not_genotype = 1u32;
            }
            iid_i += UNITS;
        }
    }
    // Reduce over the cube: within each plane, then across planes.
    let mut partials = Shared::<[u32]>::new_slice((4u32 * UNITS) as usize);
    let plane = UNIT_POS_X / PLANE_DIM;
    let count_plane = plane_sum(count);
    let sum_plane = plane_sum(sum);
    let sum2_plane = plane_sum(sum2);
    let not_genotype_plane = plane_max(not_genotype);
    if UNIT_POS_PLANE == 0u32 {
        partials[(4u32 * plane) as usize] = count_plane;
        partials[(4u32 * plane + 1u32) as usize] = sum_plane;
        partials[(4u32 * plane + 2u32) as usize] = sum2_plane;
        partials[(4u32 * plane + 3u32) as usize] = not_genotype_plane;
    }
    sync_cube();
    if UNIT_POS_X == 0u32 && sid_j < sid_count {
        let planes = UNITS / PLANE_DIM;
        let mut total_count = 0u32;
        let mut total_sum = 0u32;
        let mut total_sum2 = 0u32;
        let mut total_not_genotype = 0u32;
        let mut p = 0u32;
        while p < planes {
            total_count += partials[(4u32 * p) as usize];
            total_sum += partials[(4u32 * p + 1u32) as usize];
            total_sum2 += partials[(4u32 * p + 2u32) as usize];
            total_not_genotype = u32::max(total_not_genotype, partials[(4u32 * p + 3u32) as usize]);
            p += 1u32;
        }
        sums[(4u32 * sid_j) as usize] = total_count;
        sums[(4u32 * sid_j + 1u32) as usize] = total_sum;
        sums[(4u32 * sid_j + 2u32) as usize] = total_sum2;
        sums[(4u32 * sid_j + 3u32) as usize] = total_not_genotype;
    }
}

/// Applies per-SNP `mean` and `factor` in place, on the same 2D grid as `decode`: missing values
/// and single-value SNPs (`snc` = 1) become 0, others `(value - mean) * factor`. Call only after
/// `snp_sums` found every value to be 0, 1, 2, or NaN. (No multiply-add
/// pattern, so the GPU rounds exactly as the CPU does.)
#[cube(launch_unchecked)]
fn apply<F: Float>(val: &mut [F], iid_count: u32, sid_count: u32, c_order: u32, mean: &[F], factor: &[F], snc: &[u32]) {
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
            let mut sid_j = slow;
            if c_order == 1u32 {
                sid_j = fast;
            }
            let index = (slow * fast_count + fast) as usize;
            let value = val[index];
            // `snp_sums` already rejected values other than 0, 1, 2, and NaN; NaN (missing) fails
            // `>=`, so missing values and single-value SNPs become 0.
            if value >= F::from_int(0) && snc[sid_j as usize] == 0u32 {
                val[index] = (value - mean[sid_j as usize]) * factor[sid_j as usize];
            } else {
                val[index] = F::from_int(0);
            }
            slow += CUBE_COUNT_Y;
        }
    }
}

/// CUDA only: copies `bytes` into the `CubeCL` buffer `handle` through two reused pinned host buffers,
/// so filling one overlaps copying the other.
#[cfg(feature = "gpu-cuda")]
fn upload_pinned_cuda(client: &Client, handle: &cubecl::server::Handle, bytes: &[u8]) -> Result<(), Box<BedErrorPlus>> {
    use cubecl::cuda::CudaRuntime;
    use cubecl_runtime::runtime::Runtime;
    use cudarc::driver::result;

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
            result::memcpy_htod_async(device_ptr + (chunk_index * CHUNK_BYTES) as u64, buffer, staging.stream)
                .expect("CUDA copy to GPU");
            result::event::record(staging.events[slot], staging.stream).expect("CUDA event record");
        }
    }
    // SAFETY: the stream is valid; waiting makes the data visible to CubeCL's later launches.
    unsafe { result::stream::synchronize(staging.stream).expect("CUDA stream sync") };
    drop(resource);
    Ok(())
}

#[cfg(feature = "gpu-cuda")]
const CHUNK_BYTES: usize = 16 << 20;

/// Two reused pinned host buffers, their events, and a stream, for chunked copies.
#[cfg(feature = "gpu-cuda")]
struct CudaStaging {
    stream: cudarc::driver::sys::CUstream,
    buffers: [*mut u8; 2],
    events: [cudarc::driver::sys::CUevent; 2],
}
// SAFETY: the raw handles are only used while holding the mutex in `cuda_staging`.
#[cfg(feature = "gpu-cuda")]
unsafe impl Send for CudaStaging {}

#[cfg(feature = "gpu-cuda")]
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
            let stream = result::stream::create(result::stream::StreamKind::NonBlocking).expect("CUDA stream");
            let buffers = [0, 1].map(|_| result::malloc_host(CHUNK_BYTES, 0).expect("pinned buffer").cast::<u8>());
            let events =
                [0, 1].map(|_| result::event::create(sys::CUevent_flags::CU_EVENT_DISABLE_TIMING).expect("CUDA event"));
            std::sync::Mutex::new(CudaStaging { stream, buffers, events })
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

    use super::{Device, GpuVal};
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

    fn assert_same<T: GpuVal>(got: &nd::Array2<T>, expected: &nd::Array2<T>, context: &str) {
        assert_eq!(got.dim(), expected.dim(), "{context}");
        for (g, e) in got.iter().zip(expected.iter()) {
            let same = (num_traits::Float::is_nan(*g) && num_traits::Float::is_nan(*e))
                || bytemuck::bytes_of(g) == bytemuck::bytes_of(e);
            assert!(same, "{context}: GPU {g:?} != CPU {e:?}");
        }
    }

    fn check<T: GpuVal>(device: &Device, bed: &mut Bed, iid_index: &[isize], sid_index: &[isize], is_a1_counted: bool, is_f: bool) -> Result<(), Box<BedErrorPlus>> {
        let mut options = ReadOptions::<T>::builder();
        options.iid_index(iid_index).sid_index(sid_index).is_a1_counted(is_a1_counted).is_f(is_f);
        let context = format!("iid {iid_index:?} sid {sid_index:?} a1 {is_a1_counted} f {is_f}");
        let expected = options.read(bed)?;
        let matrix = options.read_gpu(bed, device)?;
        assert_eq!(matrix.is_f(), is_f);
        assert_same(&matrix.to_host()?, &expected, &context);

        for dist in [Dist::Unit, Dist::Beta { a: 1.0, b: 25.0 }, Dist::Beta { a: 2.0, b: 3.0 }] {
            let mut cpu_val = expected.clone();
            let mut cpu_stats = nd::Array2::<T>::zeros((sid_index.len(), 2));
            let cpu_result = impute_and_zero_mean_snps(&mut cpu_val.view_mut(), &dist, true, false, &mut cpu_stats.view_mut());
            let mut matrix = options.read_gpu(bed, device)?;
            let gpu_result = match dist {
                Dist::Unit => matrix.standardize_unit(),
                Dist::Beta { a, b } => matrix.standardize_beta(a, b),
            };
            match (cpu_result, gpu_result) {
                (Ok(()), Ok(gpu_stats)) => {
                    assert_same(&gpu_stats, &cpu_stats, &format!("{context} stats {dist:?}"));
                    assert_same(&matrix.to_host()?, &cpu_val, &format!("{context} values {dist:?}"));
                }
                (Err(cpu_error), Err(gpu_error)) => {
                    assert_eq!(gpu_error.to_string(), cpu_error.to_string(), "{context} {dist:?}");
                }
                (cpu_result, gpu_result) => panic!("{context} {dist:?}: CPU {cpu_result:?}, GPU {:?}", gpu_result.map(|_| ())),
            }
        }
        Ok(())
    }

    #[test]
    fn gpu_matches_cpu() -> Result<(), Box<BedErrorPlus>> {
        let device = Device::new()?;
        for file in FILES {
            let mut bed = Bed::new(sample_bed_file(file)?)?;
            let (iid_count, sid_count) = bed.dim()?;
            for iid_index in index_cases(iid_count) {
                for sid_index in index_cases(sid_count) {
                    for is_a1_counted in [true, false] {
                        for is_f in [true, false] {
                            check::<f32>(&device, &mut bed, &iid_index, &sid_index, is_a1_counted, is_f)?;
                            check::<f64>(&device, &mut bed, &iid_index, &sid_index, is_a1_counted, is_f)?;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn gpu_standardize_rejects_non_genotype_values() -> Result<(), Box<BedErrorPlus>> {
        // A custom missing value (-1) is not a genotype, so the exact integer sums can't be used.
        let device = Device::new()?;
        let mut bed = Bed::new(sample_bed_file("some_missing.bed")?)?;
        let mut matrix = ReadOptions::builder().missing_value(-1.0f32).read_gpu(&mut bed, &device)?;
        let error = matrix.standardize_unit().err().expect("should fail");
        assert!(matches!(*error, BedErrorPlus::GpuError(super::Error::NotGenotype)), "{error:?}");
        Ok(())
    }

    #[test]
    fn gpu_index_errors_match_cpu() -> Result<(), Box<BedErrorPlus>> {
        let device = Device::new()?;
        let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
        for (iid_index, sid_index) in [(vec![3isize], vec![0isize]), (vec![-4], vec![0]), (vec![0], vec![4]), (vec![0], vec![-5])] {
            let mut options = ReadOptions::<f32>::builder();
            options.iid_index(&iid_index).sid_index(&sid_index);
            let cpu_error = options.read(&mut bed).expect_err("CPU should fail");
            let gpu_error = options.read_gpu(&mut bed, &device).err().expect("GPU should fail");
            assert_eq!(gpu_error.to_string(), cpu_error.to_string());
        }
        Ok(())
    }
}
