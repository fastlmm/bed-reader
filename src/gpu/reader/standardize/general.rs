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

pub(super) fn statistics<TVal: CubeClFloat>(
    matrix: &CubeClMatrix<TVal>,
    genotype_sums: &Handle,
) -> Result<Vec<[f64; 3]>, Box<gpu::Error>> {
    let (iid_count, sid_count) = matrix.dim();
    let statistics_handle = matrix.client.empty(3 * sid_count * size_of::<f64>());
    // SAFETY: one cube per SNP (`sid_index < sid_count`). Matrix offsets are below
    // iid_count * sid_count, already checked by the reader. The input flags have four words per
    // SNP and the output has three f64 values per SNP, written only by unit 0 of each cube.
    unsafe {
        centered::launch_unchecked::<TVal>(
            &matrix.client,
            CubeCount::Static(to_u32(sid_count)?, 1, 1),
            CubeDim::new_1d(UNITS),
            BufferArg::from_raw_parts(matrix.handle.clone(), iid_count * sid_count),
            to_u32(iid_count)?,
            to_u32(sid_count)?,
            u32::from(!matrix.is_f),
            BufferArg::from_raw_parts(genotype_sums.clone(), 4 * sid_count),
            BufferArg::from_raw_parts(statistics_handle.clone(), 3 * sid_count),
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
    // A centered sum of squares cannot be negative. Do not turn invalid or
    // overflowed moments into a constant column as the legacy CPU formula does.
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
/// order) and their partial sums are combined in shared memory.
#[cube(launch_unchecked)]
fn centered<F: Float>(
    values: &[F],
    iid_count: u32,
    sid_count: u32,
    c_order: u32,
    genotype_sums: &[u32],
    statistics: &mut [f64],
) {
    let sid_index = CUBE_POS_X;
    // Uniform across the cube, so every unit reaches the same `sync_cube` calls below.
    let mut active = false;
    if sid_index < sid_count {
        active = genotype_sums[(4u32 * sid_index + 3u32) as usize] != 0u32;
    }

    // Pass 1: count and sum of the observed (non-NaN) values.
    let mut count = 0u32;
    let mut sum = 0.0f64;
    if active {
        let mut iid_index = UNIT_POS_X;
        while iid_index < iid_count {
            let mut index = sid_index * iid_count + iid_index;
            if c_order == 1u32 {
                index = iid_index * sid_count + sid_index;
            }
            let value = f64::cast_from(values[index as usize]);
            // Ordered comparisons exclude NaN while retaining negative values.
            if value >= 0.0f64 || value < 0.0f64 {
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
    counts[UNIT_POS_X as usize] = count;
    sums[UNIT_POS_X as usize] = sum;
    sync_cube();
    // A runtime variable (a literal start), since a value computed from the constant `UNITS`
    // stays a compile-time constant that the loop could not halve.
    let mut stride = 0u32;
    stride += UNITS / 2u32;
    while stride > 0u32 {
        if UNIT_POS_X < stride {
            counts[UNIT_POS_X as usize] += counts[(UNIT_POS_X + stride) as usize];
            sums[UNIT_POS_X as usize] += sums[(UNIT_POS_X + stride) as usize];
        }
        sync_cube();
        stride /= 2u32;
    }
    let total_count = counts[0];
    // NaN when no value is observed; the host reports that SNP as having no individuals.
    let mean = sums[0] / f64::cast_from(total_count);

    // Pass 2: sum of squared deviations from the mean.
    let mut squares = 0.0f64;
    if active {
        let mut iid_index = UNIT_POS_X;
        while iid_index < iid_count {
            let mut index = sid_index * iid_count + iid_index;
            if c_order == 1u32 {
                index = iid_index * sid_count + sid_index;
            }
            let value = f64::cast_from(values[index as usize]);
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
