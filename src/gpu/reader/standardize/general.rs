//! Centered statistics for floating-point columns outside the exact genotype path.
use cubecl::{prelude::*, server::Handle};
use num_traits::Float as HostFloat;

use super::{CubeClFloat, Error, UNITS};
use crate::{
    gpu::{
        self,
        reader::{to_u32, CubeClMatrix},
    },
    BedError, Dist,
};

/// Per-SNP `[count, mean, variance]` of the observed (non-NaN) values, computed on the GPU in f64
/// for the SNPs whose `genotype_sums` flag (word 3 of 4) marks a non-genotype value; all zeros for
/// the others, whose statistics come from the exact integer sums instead.
pub(super) fn statistics<TVal: CubeClFloat>(
    matrix: &CubeClMatrix<TVal>,
    genotype_sums: &Handle,
) -> Result<Vec<[f64; 3]>, Box<gpu::Error>> {
    let (iid_count, sid_count) = matrix.dim();
    let statistics_handle = matrix.client.empty(3 * sid_count * size_of::<f64>());
    // SAFETY: one cube per SNP (`sid_index < sid_count`). Matrix offsets are below
    // iid_count * sid_count, which fits u32 by `CubeClMatrix`'s invariant. The input flags have four words per
    // SNP and the output has three f64 values per SNP, written only by unit 0 of each cube.
    unsafe {
        centered::launch_unchecked::<TVal>(
            &matrix.client,
            CubeCount::Static(to_u32(sid_count)?, 1, 1),
            CubeDim::new_1d(UNITS),
            BufferArg::from_raw_parts(matrix.handle.clone(), iid_count * sid_count),
            to_u32(iid_count)?,
            to_u32(sid_count)?,
            BufferArg::from_raw_parts(genotype_sums.clone(), 4 * sid_count),
            BufferArg::from_raw_parts(statistics_handle.clone(), 3 * sid_count),
            matrix.is_f,
        );
    }
    let bytes = matrix
        .client
        .read_one(statistics_handle)
        .map_err(gpu::Error::from)?;
    Ok(bytemuck::cast_slice::<u8, [f64; 3]>(&bytes[..3 * sid_count * size_of::<f64>()]).to_vec())
}

pub(super) fn mean_std<TVal: CubeClFloat>(
    [count, mean, variance]: [f64; 3],
    dist: &Dist,
) -> Result<(TVal, TVal), Error> {
    if count == 0.0 {
        return Err(BedError::NoIndividuals.into());
    }
    let mean = TVal::from_f64(mean).expect("f64 converts to f32/f64");
    if HostFloat::is_nan(mean)
        || (matches!(dist, Dist::Beta { .. })
            && (mean < TVal::zero() || mean > TVal::from_f64(2.0).expect("two converts")))
    {
        return Err(BedError::IllegalSnpMean.into());
    }
    // A centered sum of squares cannot be negative, so a negative or non-finite variance means the
    // moments overflowed or were invalid: report it rather than treating the SNP as constant.
    if !variance.is_finite() || variance < 0.0 {
        return Err(Error::InvalidVariance(variance));
    }
    let deviation = if variance == 0.0 {
        f64::INFINITY
    } else {
        variance.sqrt()
    };
    let deviation = TVal::from_f64(deviation).expect("f64 converts to f32/f64");
    if variance > 0.0 && (!HostFloat::is_finite(deviation) || deviation <= TVal::zero()) {
        return Err(Error::InvalidVariance(variance));
    }
    Ok((mean, deviation))
}

/// One cube per SNP, like `snp_sums`. Only non-genotype columns do any work. Two passes in f64:
/// first the mean, then the sum of squared deviations from it, which avoids the cancellation in
/// E[x*x] - E[x]**2. In each pass the cube's units stride over individuals (coalesced in F
/// order) and their partial sums are combined in shared memory. Pass 1 also finds the observed
/// minimum and maximum: a rounded mean of equal values (three 0.1s average to
/// 0.10000000000000002) would otherwise give a constant SNP a tiny positive variance, so when
/// they are equal the mean is that value and every deviation is exactly zero.
#[cube(launch_unchecked)]
fn centered<F: Float>(
    values: &[F],
    iid_count: u32,
    sid_count: u32,
    genotype_sums: &[u32],
    statistics: &mut [f64],
    #[comptime] is_f: bool,
) {
    let sid_index = CUBE_POS_X;
    // Uniform across the cube, so every unit reaches the same `sync_cube` calls below.
    let mut active = false;
    if sid_index < sid_count {
        active = genotype_sums[(4u32 * sid_index + 3u32) as usize] != 0u32;
    }

    // Pass 1: count, sum, minimum, and maximum of the observed (non-NaN) values. A unit that
    // observes nothing keeps an infinite minimum and negative-infinite maximum, which merge away.
    let mut count = 0u32;
    let mut sum = 0.0f64;
    let mut minimum = f64::INFINITY.runtime();
    let mut maximum = f64::NEG_INFINITY.runtime();
    if active {
        let mut iid_index = UNIT_POS_X;
        while iid_index < iid_count {
            let value = value_at(values, iid_index, sid_index, iid_count, sid_count, is_f);
            // Ordered comparisons exclude NaN while retaining negative values.
            if value >= 0.0f64 || value < 0.0f64 {
                minimum = f64::min(minimum, value);
                maximum = f64::max(maximum, value);
                count += 1u32;
                sum += value;
            }
            if iid_count - iid_index <= UNITS {
                break;
            }
            iid_index += UNITS;
        }
    }
    let mut counts = Shared::<[u32]>::new_slice(UNITS as usize);
    let mut sums = Shared::<[f64]>::new_slice(UNITS as usize);
    let mut minimums = Shared::<[f64]>::new_slice(UNITS as usize);
    let mut maximums = Shared::<[f64]>::new_slice(UNITS as usize);
    counts[UNIT_POS_X as usize] = count;
    sums[UNIT_POS_X as usize] = sum;
    minimums[UNIT_POS_X as usize] = minimum;
    maximums[UNIT_POS_X as usize] = maximum;
    sync_cube();
    // A runtime variable (a literal start), since a value computed from the constant `UNITS`
    // stays a compile-time constant that the loop could not halve.
    let mut stride = 0u32;
    stride += UNITS / 2u32;
    while stride > 0u32 {
        if UNIT_POS_X < stride {
            let unit = UNIT_POS_X as usize;
            let other = (UNIT_POS_X + stride) as usize;
            counts[unit] += counts[other];
            sums[unit] += sums[other];
            minimums[unit] = f64::min(minimums[unit], minimums[other]);
            maximums[unit] = f64::max(maximums[unit], maximums[other]);
        }
        sync_cube();
        stride /= 2u32;
    }
    let total_count = counts[0];
    // NaN when no value is observed; the host reports that SNP as having no individuals.
    let mut mean = sums[0] / f64::cast_from(total_count);
    // With any observation the minimum never exceeds the maximum, so `>=` means they are equal.
    if total_count != 0u32 && minimums[0] >= maximums[0] {
        mean = minimums[0];
    }

    // Pass 2: sum of squared deviations from the mean.
    let mut squares = 0.0f64;
    if active {
        let mut iid_index = UNIT_POS_X;
        while iid_index < iid_count {
            let value = value_at(values, iid_index, sid_index, iid_count, sid_count, is_f);
            if value >= 0.0f64 || value < 0.0f64 {
                let difference = value - mean;
                squares += difference * difference;
            }
            if iid_count - iid_index <= UNITS {
                break;
            }
            iid_index += UNITS;
        }
    }
    let mut square_sums = Shared::<[f64]>::new_slice(UNITS as usize);
    square_sums[UNIT_POS_X as usize] = squares;
    sync_cube();
    stride = 0u32;
    stride += UNITS / 2u32;
    while stride > 0u32 {
        if UNIT_POS_X < stride {
            square_sums[UNIT_POS_X as usize] += square_sums[(UNIT_POS_X + stride) as usize];
        }
        sync_cube();
        stride /= 2u32;
    }

    if UNIT_POS_X == 0u32 && sid_index < sid_count {
        if active {
            statistics[(3u32 * sid_index) as usize] = f64::cast_from(total_count);
            statistics[(3u32 * sid_index + 1u32) as usize] = mean;
            statistics[(3u32 * sid_index + 2u32) as usize] =
                square_sums[0] / f64::cast_from(total_count);
        } else {
            statistics[(3u32 * sid_index) as usize] = 0.0f64;
            statistics[(3u32 * sid_index + 1u32) as usize] = 0.0f64;
            statistics[(3u32 * sid_index + 2u32) as usize] = 0.0f64;
        }
    }
}

/// One value, as f64, at the given individual and SNP in either memory order.
#[cube]
fn value_at<F: Float>(
    values: &[F],
    iid_index: u32,
    sid_index: u32,
    iid_count: u32,
    sid_count: u32,
    #[comptime] is_f: bool,
) -> f64 {
    let mut index = iid_index * sid_count + sid_index;
    if is_f {
        index = sid_index * iid_count + iid_index;
    }
    f64::cast_from(values[index as usize])
}
