//! Centered statistics for floating-point columns outside the exact genotype path.
use cubecl::{prelude::*, server::Handle};
use num_traits::Float as HostFloat;

use super::{CubeClFloat, UNITS};
use crate::{
    backend::{to_u32, CubeClMatrix, Error},
    BedError, BedErrorPlus, Dist,
};

pub(super) fn statistics<TVal: CubeClFloat>(
    matrix: &CubeClMatrix<TVal>,
    genotype_sums: &Handle,
) -> Result<Vec<[f64; 3]>, Box<BedErrorPlus>> {
    let (iid_count, sid_count) = matrix.dim();
    let statistics_handle = matrix.client.empty(3 * sid_count * size_of::<f64>());
    // SAFETY: one unit per SNP, guarded by sid_index < sid_count. Matrix offsets
    // are below iid_count * sid_count, already checked by the reader. The input
    // flags have four words per SNP and the output has three f64 values per SNP.
    unsafe {
        centered::launch_unchecked::<TVal>(
            &matrix.client,
            CubeCount::Static(to_u32(sid_count)?.div_ceil(UNITS), 1, 1),
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
        .map_err(Error::from)?;
    Ok(bytemuck::cast_slice::<u8, [f64; 3]>(&bytes[..3 * sid_count * size_of::<f64>()]).to_vec())
}

pub(super) fn mean_std<TVal: CubeClFloat>(
    [count, mean, variance]: [f64; 3],
    dist: &Dist,
) -> Result<(TVal, TVal), BedErrorPlus> {
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
        return Err(BedErrorPlus::GpuError(Error::InvalidVariance(variance)));
    }
    let deviation = if variance == 0.0 {
        f64::INFINITY
    } else {
        variance.sqrt()
    };
    let deviation = TVal::from_f64(deviation).expect("f64 converts to f32/f64");
    if variance > 0.0 && (!HostFloat::is_finite(deviation) || deviation <= TVal::zero()) {
        return Err(BedErrorPlus::GpuError(Error::InvalidVariance(variance)));
    }
    Ok((mean, deviation))
}

/// One GPU unit per non-genotype column. Subtract the first observed value
/// before summing, then sum squared deviations about that shifted mean in a
/// second pass, both in f64. This avoids cancellation in E[x*x] - E[x]**2 and
/// avoids summing a large common offset. Only per-SNP output is allocated.
#[cube(launch_unchecked)]
fn centered<F: Float>(
    values: &[F],
    iid_count: u32,
    sid_count: u32,
    c_order: u32,
    genotype_sums: &[u32],
    statistics: &mut [f64],
) {
    let sid_index = CUBE_POS_X * UNITS + UNIT_POS_X;
    if sid_index < sid_count && genotype_sums[(4u32 * sid_index + 3u32) as usize] != 0u32 {
        let mut count = 0u32;
        let mut anchor = 0.0f64;
        let mut shifted_sum = 0.0f64;
        let mut iid_index = 0u32;
        while iid_index < iid_count {
            let mut index = sid_index * iid_count + iid_index;
            if c_order == 1u32 {
                index = iid_index * sid_count + sid_index;
            }
            let value = f64::cast_from(values[index as usize]);
            // Ordered comparisons exclude NaN while retaining negative values.
            if value >= 0.0f64 || value < 0.0f64 {
                if count == 0u32 {
                    anchor = value;
                }
                shifted_sum += value - anchor;
                count += 1u32;
            }
            iid_index += 1u32;
        }
        let shifted_mean = shifted_sum / f64::cast_from(count);
        let mut centered_sum = 0.0f64;
        iid_index = 0u32;
        while iid_index < iid_count {
            let mut index = sid_index * iid_count + iid_index;
            if c_order == 1u32 {
                index = iid_index * sid_count + sid_index;
            }
            let value = f64::cast_from(values[index as usize]);
            if value >= 0.0f64 || value < 0.0f64 {
                let difference = (value - anchor) - shifted_mean;
                centered_sum += difference * difference;
            }
            iid_index += 1u32;
        }
        statistics[(3u32 * sid_index) as usize] = f64::cast_from(count);
        statistics[(3u32 * sid_index + 1u32) as usize] = anchor + shifted_mean;
        statistics[(3u32 * sid_index + 2u32) as usize] = centered_sum / f64::cast_from(count);
    } else if sid_index < sid_count {
        statistics[(3u32 * sid_index) as usize] = 0.0f64;
        statistics[(3u32 * sid_index + 1u32) as usize] = 0.0f64;
        statistics[(3u32 * sid_index + 2u32) as usize] = 0.0f64;
    }
}
