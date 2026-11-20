//! Standardize a [`CubeClMatrix`] in place on the GPU. Validated genotype columns
//! retain exact CPU-compatible integer statistics; other floating-point columns
//! use two-pass, centered f64 reductions to avoid variance cancellation.
//!
//! Like the CPU standardizer, this is for the Python bindings (`PySnpTools`), not the public Rust API.
mod general;

use cubecl::prelude::*;
use ndarray as nd;

use super::{memory_dimensions, to_u32, CubeClMatrix, CubeClVal, MAX_GRID_Y, UNITS};
use crate::{find_factor, gpu, snp_mean_std, BedError, Dist};

/// Internal standardizer failures; invalid variance is not part of the public GPU API.
/// Boxed conversions are explicit for the same coherence reasons as [`gpu::Error`].
#[derive(thiserror::Error, Debug)]
pub(crate) enum Error {
    #[error(transparent)]
    Gpu(#[from] Box<gpu::Error>),
    #[error(transparent)]
    Bed(#[from] BedError),
    #[error("Invalid GPU variance: {0}")]
    InvalidVariance(f64),
}

impl From<Box<gpu::Error>> for Box<Error> {
    fn from(error: Box<gpu::Error>) -> Self {
        Box::new(Error::from(error))
    }
}

impl From<gpu::Error> for Box<Error> {
    fn from(error: gpu::Error) -> Self {
        Box::new(Error::Gpu(Box::new(error)))
    }
}

impl From<BedError> for Box<Error> {
    fn from(error: BedError) -> Self {
        Box::new(Error::from(error))
    }
}

/// A float type the GPU code standardizes: `f32` or `f64`.
pub(crate) trait CubeClFloat:
    CubeClVal + Float + num_traits::Float + num_traits::ToPrimitive + num_traits::FromPrimitive
{
}
impl CubeClFloat for f32 {}
impl CubeClFloat for f64 {}

/// Standardize in place to zero mean and unit variance (missing values become 0), matching the
/// CPU genotype standardization exactly. Fractional inputs use centered GPU statistics.
/// Returns each SNP's mean and standard deviation (`sid_count x 2`).
///
/// # Errors
/// [`BedError::NoIndividuals`] for a SNP with no observed values;
/// [`Error::InvalidVariance`] if floating-point moments overflow or become invalid.
pub(crate) fn standardize_unit<TVal: CubeClFloat>(
    matrix: &mut CubeClMatrix<TVal>,
) -> Result<nd::Array2<TVal>, Box<Error>> {
    standardize(matrix, &Dist::Unit)
}

/// Standardize in place with Beta(`a`, `b`) weights (missing values become 0). Genotype columns
/// match CPU statistics exactly; fractional columns use centered GPU statistics.
/// Returns each SNP's mean and standard deviation (`sid_count x 2`).
///
/// # Errors
/// As [`standardize_unit`], plus [`BedError::IllegalSnpMean`] and
/// [`BedError::CannotCreateBetaDist`].
pub(crate) fn standardize_beta<TVal: CubeClFloat>(
    matrix: &mut CubeClMatrix<TVal>,
    a: f64,
    b: f64,
) -> Result<nd::Array2<TVal>, Box<Error>> {
    standardize(matrix, &Dist::Beta { a, b })
}

fn standardize<TVal: CubeClFloat>(
    matrix: &mut CubeClMatrix<TVal>,
    dist: &Dist,
) -> Result<nd::Array2<TVal>, Box<Error>> {
    let (iid_count, sid_count) = (matrix.iid_count, matrix.sid_count);
    let mut stats = nd::Array2::<TVal>::zeros((sid_count, 2));
    if sid_count == 0 {
        return Ok(stats);
    }
    let element_count = iid_count * sid_count;
    // Each observed genotype contributes at most four to sum2. Both the sums
    // and the four-word-per-SNP output are indexed with u32 in the kernels.
    for count in [iid_count, sid_count] {
        let times_four = u64::from(to_u32(count)?) * 4;
        u32::try_from(times_four).map_err(|_| gpu::Error::TooLarge(times_four))?;
    }
    let c_order = u32::from(!matrix.is_f);
    let client = &matrix.client;
    // CUDA has plane operations; only CubeCL's CPU runtime (internal tests) lacks them.
    let use_plane = client
        .features()
        .plane
        .contains(cubecl::ir::features::Plane::Ops);
    let sums_handle = launch_snp_sums(matrix, use_plane)?;
    let sums_bytes = client
        .read_one(sums_handle.clone())
        .map_err(gpu::Error::from)?;
    let sums: &[u32] = bytemuck::cast_slice(&sums_bytes[..4 * sid_count * 4]);

    let general_statistics = if sums.chunks_exact(4).any(|column| column[3] != 0) {
        Some(general::statistics(matrix, &sums_handle)?)
    } else {
        None
    };

    // Finish each SNP's statistics with the CPU code's functions. The integer sums convert to
    // `TVal` exactly whenever the CPU's own float sums are exact (below 2^24 for f32).
    let two = <TVal as From<i8>>::from(2);
    let mut mean = Vec::with_capacity(sid_count);
    let mut factor = Vec::with_capacity(sid_count);
    let mut snc = Vec::with_capacity(sid_count);
    for (sid_index, (column, mut stats_row)) in sums
        .chunks_exact(4)
        .zip(stats.axis_iter_mut(nd::Axis(0)))
        .enumerate()
    {
        let &[count, sum, sum2, not_genotype] = column else {
            unreachable!("chunks_exact(4) produces four statistics per SNP");
        };
        let (mean_s, std) = if not_genotype == 0 {
            snp_mean_std(
                from_u32::<TVal>(count),
                from_u32(sum),
                from_u32(sum2),
                dist,
                two,
            )?
        } else {
            general::mean_std::<TVal>(
                general_statistics
                    .as_ref()
                    .expect("non-genotype statistics")[sid_index],
                dist,
            )?
        };
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
        let (fast_count, slow_count) = memory_dimensions(iid_count, sid_count, matrix.is_f);
        // SAFETY: as for `decode`: writes only `slow * fast_count + fast` < `element_count`;
        // reads the per-SNP arrays at `sid_j < sid_count`, their length.
        unsafe {
            apply::launch_unchecked::<TVal>(
                client,
                CubeCount::Static(
                    to_u32(fast_count)?.div_ceil(UNITS),
                    to_u32(slow_count.min(MAX_GRID_Y))?,
                    1,
                ),
                CubeDim::new_1d(UNITS),
                BufferArg::from_raw_parts(matrix.handle.clone(), element_count),
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

/// Launches `snp_sums` into a new `4 * sid_count`-word buffer and returns it. `use_plane` picks
/// the plane (warp) reduction, which needs `Plane::Ops`, or the shared-memory tree reduction;
/// both give identical results. The caller has checked that `4 * iid_count`, `4 * sid_count`,
/// and `iid_count * sid_count` fit u32.
fn launch_snp_sums<TVal: CubeClFloat>(
    matrix: &CubeClMatrix<TVal>,
    use_plane: bool,
) -> Result<cubecl::server::Handle, Box<Error>> {
    let (iid_count, sid_count) = (matrix.iid_count, matrix.sid_count);
    let element_count = iid_count * sid_count;
    let c_order = u32::from(!matrix.is_f);
    let client = &matrix.client;
    let sums_handle = client.empty(4 * sid_count * 4);
    // SAFETY: one cube per SNP (`sid_j < sid_count`); it reads the matrix only at indexes below
    // `iid_count * sid_count` (= `element_count`, which fits u32) and writes
    // `sums[4 * sid_j .. 4 * sid_j + 4]`, inside the `4 * sid_count` words allocated.
    unsafe {
        snp_sums::launch_unchecked::<TVal>(
            client,
            CubeCount::Static(to_u32(sid_count)?, 1, 1),
            CubeDim::new_1d(UNITS),
            BufferArg::from_raw_parts(
                matrix.handle.clone(),
                to_u32(element_count.max(1))? as usize,
            ),
            to_u32(iid_count)?,
            to_u32(sid_count)?,
            c_order,
            BufferArg::from_raw_parts(sums_handle.clone(), 4 * sid_count),
            use_plane,
        );
    }
    Ok(sums_handle)
}

fn from_u32<TVal: CubeClFloat>(value: u32) -> TVal {
    // f32/f64 from u32 is the CPU's sum of the same integers whenever that sum is exact.
    <TVal as num_traits::NumCast>::from(value).expect("u32 converts to f32/f64")
}

/// Per-SNP sums over observed (non-NaN) values: one cube per SNP. Writes `count`, `sum`,
/// `sum of squares`, and a flag set if any value is not 0, 1, or 2 to `sums[4 * sid_j ..]`. The
/// values are small integers, so these integer sums are exact in any order.
#[cube(launch_unchecked)]
fn snp_sums<F: Float>(
    val: &[F],
    iid_count: u32,
    sid_count: u32,
    c_order: u32,
    sums: &mut [u32],
    #[comptime] use_plane: bool,
) {
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
            if value >= F::from_int(0) && value <= F::from_int(2) {
                let genotype = u32::cast_from(value);
                if F::cast_from(genotype) == value {
                    count += 1u32;
                    sum += genotype;
                    sum2 += genotype * genotype;
                } else {
                    not_genotype = 1u32;
                }
            } else if value < F::from_int(0) || value > F::from_int(2) {
                // Never cast a foreign float to u32: it may overflow before
                // the general floating-point path gets a chance to handle it.
                not_genotype = 1u32;
            }
            if iid_count - iid_i <= UNITS {
                break;
            }
            iid_i += UNITS;
        }
    }
    let mut partials = Shared::<[u32]>::new_slice((4u32 * UNITS) as usize);
    if use_plane {
        // Reduce over the cube: within each plane, then across planes.
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
                total_not_genotype =
                    u32::max(total_not_genotype, partials[(4u32 * p + 3u32) as usize]);
                p += 1u32;
            }
            sums[(4u32 * sid_j) as usize] = total_count;
            sums[(4u32 * sid_j + 1u32) as usize] = total_sum;
            sums[(4u32 * sid_j + 2u32) as usize] = total_sum2;
            sums[(4u32 * sid_j + 3u32) as usize] = total_not_genotype;
        }
    } else {
        // Runtimes without plane operations (CubeCL's CPU runtime, used by this crate's tests):
        // a shared-memory tree reduction. Integer sums and max are exact in any order, so this
        // gives the same results as the plane path above.
        partials[(4u32 * UNIT_POS_X) as usize] = count;
        partials[(4u32 * UNIT_POS_X + 1u32) as usize] = sum;
        partials[(4u32 * UNIT_POS_X + 2u32) as usize] = sum2;
        partials[(4u32 * UNIT_POS_X + 3u32) as usize] = not_genotype;
        sync_cube();
        // A runtime variable, as in `general::centered`, so the loop can halve it.
        let mut stride = 0u32;
        stride += UNITS / 2u32;
        while stride > 0u32 {
            if UNIT_POS_X < stride {
                let mine = 4u32 * UNIT_POS_X;
                let other = 4u32 * (UNIT_POS_X + stride);
                partials[mine as usize] += partials[other as usize];
                partials[(mine + 1u32) as usize] += partials[(other + 1u32) as usize];
                partials[(mine + 2u32) as usize] += partials[(other + 2u32) as usize];
                partials[(mine + 3u32) as usize] = u32::max(
                    partials[(mine + 3u32) as usize],
                    partials[(other + 3u32) as usize],
                );
            }
            sync_cube();
            stride /= 2u32;
        }
        if UNIT_POS_X == 0u32 && sid_j < sid_count {
            sums[(4u32 * sid_j) as usize] = partials[0];
            sums[(4u32 * sid_j + 1u32) as usize] = partials[1];
            sums[(4u32 * sid_j + 2u32) as usize] = partials[2];
            sums[(4u32 * sid_j + 3u32) as usize] = partials[3];
        }
    }
}

/// Applies per-SNP `mean` and `factor` in place, on the same 2D grid as `decode`: missing values
/// and single-value SNPs (`snc` = 1) become 0, others `(value - mean) * factor`. Call only after
/// the per-SNP statistics have been validated. (No multiply-add pattern, so
/// genotype columns continue to round exactly as the CPU does.)
#[cube(launch_unchecked)]
fn apply<F: Float>(
    val: &mut [F],
    iid_count: u32,
    sid_count: u32,
    c_order: u32,
    mean: &[F],
    factor: &[F],
    snc: &[u32],
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
            let mut sid_j = slow;
            if c_order == 1u32 {
                sid_j = fast;
            }
            let index = (slow * fast_count + fast) as usize;
            let value = val[index];
            // NaN fails both ordered comparisons; negative finite values must
            // still be transformed. Single-value SNPs always become zero.
            if (value >= F::from_int(0) || value < F::from_int(0)) && snc[sid_j as usize] == 0u32 {
                val[index] = (value - mean[sid_j as usize]) * factor[sid_j as usize];
            } else {
                val[index] = F::from_int(0);
            }
            if slow_count - slow <= CUBE_COUNT_Y {
                break;
            }
            slow += CUBE_COUNT_Y;
        }
    }
}

// CUDA only: CubeCL's CPU runtime (internal feature `_cubecl-cpu-test`) cannot compile the plane
// path, so there only the tree path runs, through the standardizer tests.
#[cfg(all(test, not(feature = "_cubecl-cpu-test")))]
mod tests {
    //! These need a GPU; without one, they fail (they do not skip).
    use std::error::Error as StandardError;

    use super::{launch_snp_sums, CubeClFloat};
    use crate::{gpu, sample_bed_file, Bed, ReadOptions};

    fn sums<TVal: CubeClFloat>(
        matrix: &super::CubeClMatrix<TVal>,
        use_plane: bool,
    ) -> Result<Vec<u8>, Box<dyn StandardError>> {
        let handle = launch_snp_sums(matrix, use_plane)?;
        let bytes = matrix.client.read_one(handle).map_err(gpu::Error::from)?;
        Ok(bytes[..16 * matrix.sid_count].to_vec())
    }

    /// The plane and tree reductions of `snp_sums` give identical sums and flags, for genotype
    /// data with NaNs and for non-genotype data (missing value -1), in both memory orders.
    #[test]
    fn snp_sums_plane_and_tree_reductions_agree() -> Result<(), Box<dyn StandardError>> {
        for file in [
            "some_missing.bed",
            "plink_sim_10s_100v_10pmiss.bed",
            "toydata.5chrom.bed",
        ] {
            for is_f in [true, false] {
                let mut bed = Bed::new(sample_bed_file(file)?)?;
                for missing in [f32::NAN, -1.0] {
                    let matrix = ReadOptions::builder()
                        .missing_value(missing)
                        .is_f(is_f)
                        .read_cubecl(&mut bed)?;
                    let plane = sums(&matrix, true)?;
                    assert_eq!(
                        plane,
                        sums(&matrix, false)?,
                        "{file} is_f={is_f} f32 {missing}"
                    );
                    // Only an observed missing value (-1, not NaN) is a non-genotype value.
                    let expected = !missing.is_nan()
                        && matrix
                            .to_cpu()?
                            .iter()
                            .any(|value| value.to_bits() == missing.to_bits());
                    let has_non_genotype = plane.chunks_exact(16).any(|column| column[12] != 0);
                    assert_eq!(
                        has_non_genotype, expected,
                        "{file} is_f={is_f} non-genotype flag"
                    );
                }
                let matrix = ReadOptions::builder()
                    .f64()
                    .is_f(is_f)
                    .read_cubecl(&mut bed)?;
                assert_eq!(
                    sums(&matrix, true)?,
                    sums(&matrix, false)?,
                    "{file} is_f={is_f} f64"
                );
            }
        }
        Ok(())
    }
}
