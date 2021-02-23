// Inspired by C++ version by Chris Widmer and Carl Kadie

use core::fmt::Debug;
use ndarray as nd;
use ndarray::ShapeBuilder;
use num_traits::{Float, FromPrimitive, ToPrimitive};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rayon::{iter::ParallelBridge, ThreadPoolBuildError};
use statrs::distribution::{Beta, Continuous};
use std::{
    convert::TryFrom,
    ops::{Div, Sub},
};
use std::{
    fs::File,
    io::{BufRead, BufWriter, Read, Write},
};
use std::{io::SeekFrom, path::PathBuf};
use std::{
    io::{BufReader, Seek},
    path::Path,
};
use thiserror::Error;

const BED_FILE_MAGIC1: u8 = 0x6C; // 0b01101100 or 'l' (lowercase 'L')
const BED_FILE_MAGIC2: u8 = 0x1B; // 0b00011011 or <esc>
const CB_HEADER_U64: u64 = 3;
const CB_HEADER_USIZE: usize = 3;

// About ndarray
//  https://docs.rs/ndarray/0.14.0/ndarray/parallel/index.html
//  https://rust-lang-nursery.github.io/rust-cookbook/concurrency/parallel.html
//  https://github.com/rust-ndarray/ndarray/blob/master/README-quick-start.md
//  https://datacrayon.com/posts/programming/rust-notebooks/multidimensional-arrays-and-operations-with-ndarray
//  https://docs.rs/ndarray/0.14.0/ndarray/doc/ndarray_for_numpy_users/index.html
//  https://docs.rs/ndarray-npy
//  https://rust-lang-nursery.github.io/rust-cookbook/science/mathematics/linear_algebra.html

/// BedErrorPlus enumerates all possible errors returned by this library.
/// Based on https://nick.groenen.me/posts/rust-error-handling/#the-library-error-type
#[derive(Error, Debug)]
pub enum BedErrorPlus {
    #[error(transparent)]
    IOError(#[from] std::io::Error),

    #[error(transparent)]
    BedError(#[from] BedError),

    #[error(transparent)]
    ThreadPoolError(#[from] ThreadPoolBuildError),
}
// https://docs.rs/thiserror/1.0.23/thiserror/
#[derive(Error, Debug, Clone)]
pub enum BedError {
    #[error("Ill-formed BED file. BED file header is incorrect or length is wrong. '{0}'")]
    IllFormed(String),

    #[error(
        "Ill-formed BED file. BED file header is incorrect. Expected mode to be 0 or 1. '{0}'"
    )]
    BadMode(String),

    #[error("Attempt to write illegal value to BED file. Only 0,1,2,missing allowed. '{0}'")]
    BadValue(String),

    #[error("No individual observed for the SNP.")]
    NoIndividuals,

    #[error("Illegal SNP mean.")]
    IllegalSnpMean,

    #[error("Index to individual larger than the number of individuals. (Index value {0})")]
    IidIndexTooBig(usize),

    #[error("Index to SNP larger than the number of SNPs. (Index value {0})")]
    SidIndexTooBig(usize),

    #[error("Length of iid_index ({0}) and sid_index ({1}) must match dimensions of output array ({2},{3}).")]
    IndexMismatch(usize, usize, usize, usize),

    #[error("Indexes ({0},{1}) too big for files")]
    IndexesTooBigForFiles(usize, usize),

    #[error("Subset: length of iid_index ({0}) and sid_index ({1}) must match dimensions of output array ({2},{3}).")]
    SubsetMismatch(usize, usize, usize, usize),

    #[error("Cannot convert beta values to/from float 64")]
    CannotConvertBetaToFromF64,

    #[error("Cannot create Beta Dist with given parameters ({0},{1})")]
    CannotCreateBetaDist(f64, f64),

    #[error("Cannot open metadata file. '{0}'")]
    CannotOpenFamOrBim(String),
}

fn read_no_alloc<TOut: Copy + Default + From<i8> + Debug + Sync + Send>(
    filename: &str,
    iid_count: usize,
    sid_count: usize,
    count_a1: bool,
    iid_index: &[usize],
    sid_index: &[usize],
    missing_value: TOut,
    val: &mut nd::ArrayViewMut2<'_, TOut>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.
) -> Result<(), BedErrorPlus> {
    let mut buf_reader = BufReader::new(File::open(filename)?);
    let mut bytes_vector: Vec<u8> = vec![0; CB_HEADER_USIZE];
    buf_reader.read_exact(&mut bytes_vector)?;
    if (BED_FILE_MAGIC1 != bytes_vector[0]) || (BED_FILE_MAGIC2 != bytes_vector[1]) {
        return Err(BedError::IllFormed(filename.to_string()).into());
    }
    match bytes_vector[2] {
        0 => {
            let mut val_t = val.view_mut().reversed_axes();
            return internal_read_no_alloc(
                buf_reader,
                filename,
                sid_count,
                iid_count,
                count_a1,
                sid_index,
                iid_index,
                missing_value,
                &mut val_t,
            );
        }
        1 => {
            return internal_read_no_alloc(
                buf_reader,
                filename,
                iid_count,
                sid_count,
                count_a1,
                iid_index,
                sid_index,
                missing_value,
                val,
            );
        }
        _ => {
            return Err(BedError::BadMode(filename.to_string()).into());
        }
    }
}

trait Max {
    fn max() -> Self;
}

impl Max for u8 {
    fn max() -> u8 {
        std::u8::MAX
    }
}

impl Max for u64 {
    fn max() -> u64 {
        std::u64::MAX
    }
}

// We make this generic instead of u64, so that we can test it via u8
fn try_div_4<T: Max + TryFrom<usize> + Sub<Output = T> + Div<Output = T> + Ord>(
    in_iid_count: usize,
    in_sid_count: usize,
    cb_header: T,
) -> Result<(usize, T), BedErrorPlus> {
    // 4 genotypes per byte so round up without overflow
    let in_iid_count_div4 = if in_iid_count > 0 {
        (in_iid_count - 1) / 4 + 1
    } else {
        0
    };
    let in_iid_count_div4_t = match T::try_from(in_iid_count_div4) {
        Ok(v) => v,
        Err(_) => return Err(BedError::IndexesTooBigForFiles(in_iid_count, in_sid_count).into()),
    };
    let in_sid_count_t = match T::try_from(in_sid_count) {
        Ok(v) => v,
        Err(_) => return Err(BedError::IndexesTooBigForFiles(in_iid_count, in_sid_count).into()),
    };

    let m: T = Max::max(); // Don't know how to move this into the next line.
    if in_sid_count > 0 && (m - cb_header) / in_sid_count_t < in_iid_count_div4_t {
        return Err(BedError::IndexesTooBigForFiles(in_iid_count, in_sid_count).into());
    }

    return Ok((in_iid_count_div4, in_iid_count_div4_t));
}

fn internal_read_no_alloc<TOut: Copy + Default + From<i8> + Debug + Sync + Send>(
    mut buf_reader: BufReader<File>,
    filename: &str,
    in_iid_count: usize,
    in_sid_count: usize,
    count_a1: bool,
    iid_index: &[usize],
    sid_index: &[usize],
    missing_value: TOut,
    out_val: &mut nd::ArrayViewMut2<'_, TOut>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.
) -> Result<(), BedErrorPlus> {
    // Find the largest in_iid_i (if any) and check its size.
    if let Some(in_max_iid_i) = iid_index.iter().max() {
        if *in_max_iid_i >= in_iid_count {
            return Err(BedError::IidIndexTooBig(*in_max_iid_i).into());
        }
    }

    let out_iid_count = iid_index.len();
    let out_sid_count = sid_index.len();

    let (in_iid_count_div4, in_iid_count_div4_u64) =
        try_div_4(in_iid_count, in_sid_count, CB_HEADER_U64)?;

    let from_two_bits_to_value = set_up_two_bits_to_value(count_a1, missing_value);

    // "as" and math is safe because of early checks
    if buf_reader.seek(SeekFrom::End(0))?
        != in_iid_count_div4_u64 * (in_sid_count as u64) + CB_HEADER_U64
    {
        return Err(BedErrorPlus::BedError(BedError::IllFormed(
            filename.to_string(),
        )));
    }

    // See https://morestina.net/blog/1432/parallel-stream-processing-with-rayon
    // Possible optimization: We could try to read only the iid info needed
    // Possible optimization: We could read snp in their input order instead of their output order
    (0..out_sid_count)
        // Read all the iid info for one snp from the disk
        .map(|out_sid_i| {
            let in_sid_i = sid_index[out_sid_i];
            if in_sid_i >= in_sid_count {
                return Err(BedErrorPlus::BedError(BedError::SidIndexTooBig(in_sid_i)));
            }
            let mut bytes_vector: Vec<u8> = vec![0; in_iid_count_div4];
            let pos: u64 = (in_sid_i as u64) * in_iid_count_div4_u64 + CB_HEADER_U64; // "as" and math is safe because of early checks
            buf_reader.seek(SeekFrom::Start(pos))?;
            buf_reader.read_exact(&mut bytes_vector)?;
            return Ok(bytes_vector);
        })
        // Zip in the column of the output array
        .zip(out_val.axis_iter_mut(nd::Axis(1)))
        // In parallel, decompress the iid info and put it in its column
        .par_bridge() // This seems faster that parallel zip
        .try_for_each(|(bytes_vector_result, mut col)| {
            match bytes_vector_result {
                Err(e) => Err(e),
                Ok(bytes_vector) => {
                    for out_iid_i in 0..out_iid_count {
                        // Possible optimization: We could pre-compute the conversion, the division, the mod, and the multiply*2
                        let in_iid_i = iid_index[out_iid_i];
                        let i_div_4 = in_iid_i / 4;
                        let i_mod_4 = in_iid_i % 4;
                        let genotype_byte: u8 = (bytes_vector[i_div_4] >> (i_mod_4 * 2)) & 0x03;
                        col[out_iid_i] = from_two_bits_to_value[genotype_byte as usize];
                    }
                    Ok(())
                }
            }
        })?;

    return Ok(());
}

fn set_up_two_bits_to_value<TOut: From<i8>>(count_a1: bool, missing_value: TOut) -> [TOut; 4] {
    let homozygous_primary_allele = TOut::from(0); // Major Allele
    let heterozygous_allele = TOut::from(1);
    let homozygous_secondary_allele = TOut::from(2); // Minor Allele

    let from_two_bits_to_value;
    if count_a1 {
        from_two_bits_to_value = [
            homozygous_secondary_allele, // look-up 0
            missing_value,               // look-up 1
            heterozygous_allele,         // look-up 2
            homozygous_primary_allele,   // look-up 3
        ];
    } else {
        from_two_bits_to_value = [
            homozygous_primary_allele,   // look-up 0
            missing_value,               // look-up 1
            heterozygous_allele,         // look-up 2
            homozygous_secondary_allele, // look-up 3
        ];
    }
    return from_two_bits_to_value;
}

// could make count_a1, etc. optional
pub fn read_with_indexes<TOut: From<i8> + Default + Copy + Debug + Sync + Send>(
    filename: &str,
    iid_index: &[usize],
    sid_index: &[usize],
    output_is_order_f: bool,
    count_a1: bool,
    missing_value: TOut,
) -> Result<nd::Array2<TOut>, BedErrorPlus> {
    let (iid_count, sid_count) = counts(filename)?;

    let shape = ShapeBuilder::set_f((iid_index.len(), sid_index.len()), output_is_order_f);
    let mut val = nd::Array2::<TOut>::default(shape);

    read_no_alloc(
        filename,
        iid_count,
        sid_count,
        count_a1,
        iid_index,
        sid_index,
        missing_value,
        &mut val.view_mut(),
    )?;

    return Ok(val);
}

pub fn read<TOut: From<i8> + Default + Copy + Debug + Sync + Send>(
    filename: &str,
    output_is_order_f: bool,
    count_a1: bool,
    missing_value: TOut,
) -> Result<nd::Array2<TOut>, BedErrorPlus> {
    let (iid_count, sid_count) = counts(filename)?;

    let iid_index: Vec<usize> = (0..iid_count).collect();
    let sid_index: Vec<usize> = (0..sid_count).collect();

    let shape = ShapeBuilder::set_f((iid_count, sid_count), output_is_order_f);
    let mut val = nd::Array2::<TOut>::default(shape);

    read_no_alloc(
        filename,
        iid_count,
        sid_count,
        count_a1,
        &iid_index,
        &sid_index,
        missing_value,
        &mut val.view_mut(),
    )?;

    return Ok(val);
}

pub fn write<T: From<i8> + Default + Copy + Debug + Sync + Send + PartialEq>(
    filename: &str,
    val: &nd::ArrayView2<'_, T>,
    count_a1: bool,
    missing: T,
) -> Result<(), BedErrorPlus> {
    let mut writer = BufWriter::new(File::create(filename)?);
    writer.write_all(&[BED_FILE_MAGIC1, BED_FILE_MAGIC2, 0x01])?;

    let zero_code = if count_a1 { 3u8 } else { 0u8 };
    let two_code = if count_a1 { 0u8 } else { 3u8 };

    let homozygous_primary_allele = T::from(0); // Major Allele
    let heterozygous_allele = T::from(1);
    let homozygous_secondary_allele = T::from(2); // Minor Allele

    let (iid_count, sid_count) = val.dim();

    // 4 genotypes per byte so round up
    let (iid_count_div4, _) = try_div_4(iid_count, sid_count, CB_HEADER_U64)?;

    let use_nan = missing != missing;
    for column in val.axis_iter(nd::Axis(1)) {
        let mut bytes_vector: Vec<u8> = vec![0; iid_count_div4]; // inits to 0
        for (iid_i, &v0) in column.iter().enumerate() {
            let genotype_byte = if v0 == homozygous_primary_allele {
                zero_code
            } else if v0 == heterozygous_allele {
                2
            } else if v0 == homozygous_secondary_allele {
                two_code
            } else if (use_nan && v0 != v0) || (!use_nan && v0 == missing) {
                1
            } else {
                return Err(BedError::BadValue(filename.to_string()).into());
            };
            // Possible optimization: We could pre-compute the conversion, the division, the mod, and the multiply*2
            let i_div_4 = iid_i / 4;
            let i_mod_4 = iid_i % 4;
            bytes_vector[i_div_4] |= genotype_byte << (i_mod_4 * 2);
        }
        writer.write_all(&bytes_vector)?;
    }
    return Ok(());
}

fn count_lines(path_buf: PathBuf) -> Result<usize, BedErrorPlus> {
    let file = match File::open(&path_buf) {
        Err(_) => {
            let string_path = path_buf.to_string_lossy().to_string();
            return Err(BedErrorPlus::BedError(BedError::CannotOpenFamOrBim(
                string_path,
            )));
        }
        Ok(file) => file,
    };
    let reader = BufReader::new(file);
    let count = reader.lines().count();
    return Ok(count);
}
pub fn counts(filename: &str) -> Result<(usize, usize), BedErrorPlus> {
    let path = Path::new(filename);
    let iid_count = count_lines(path.with_extension("fam"))?;
    let sid_count = count_lines(path.with_extension("bim"))?;
    return Ok((iid_count, sid_count));
}

pub fn matrix_subset_no_alloc<
    TIn: Copy + Default + Debug + Sync + Send + Sized,
    TOut: Copy + Default + Debug + Sync + Send + From<TIn>,
>(
    in_val: &nd::ArrayView3<'_, TIn>,
    iid_index: &[usize],
    sid_index: &[usize],
    out_val: &mut nd::ArrayViewMut3<'_, TOut>,
) -> Result<(), BedErrorPlus> {
    let out_iid_count = iid_index.len();
    let out_sid_count = sid_index.len();
    let did_count = in_val.dim().2;

    if (out_iid_count, out_sid_count, did_count) != out_val.dim() {
        return Err(BedError::SubsetMismatch(
            out_iid_count,
            out_sid_count,
            out_val.dim().0,
            out_val.dim().1,
        )
        .into());
    }

    // If output is F-order (or in general if iid stride is no more than sid_stride)
    if out_val.stride_of(nd::Axis(0)) <= out_val.stride_of(nd::Axis(1)) {
        // (No error are possible in the par_azip, so don't have to collect and check them)
        nd::par_azip!((mut out_col in out_val.axis_iter_mut(nd::Axis(1)),
                    in_sid_i_pr in sid_index) {
            let in_col = in_val.index_axis(nd::Axis(1), *in_sid_i_pr);
            for did_i in 0..did_count
            {
                for (out_iid_i, in_iid_i_ptr) in iid_index.iter().enumerate() {
                    out_col[(out_iid_i,did_i)] = in_col[(*in_iid_i_ptr,did_i)].into();
                }
            }
        });
        return Ok(());
    } else {
        //If output is C-order, transpose input and output and recurse
        let in_val_t = in_val.view().permuted_axes([1, 0, 2]);
        let mut out_val_t = out_val.view_mut().permuted_axes([1, 0, 2]);
        return matrix_subset_no_alloc(&in_val_t, &sid_index, &iid_index, &mut out_val_t);
    }
}

pub enum Dist {
    Unit,
    Beta { a: f64, b: f64 },
}

pub fn impute_and_zero_mean_snps<
    T: Default + Copy + Debug + Sync + Send + Float + ToPrimitive + FromPrimitive,
>(
    val: &mut nd::ArrayViewMut2<'_, T>,
    dist: Dist,
    apply_in_place: bool,
    use_stats: bool,
    stats: &mut nd::ArrayViewMut2<'_, T>,
) -> Result<(), BedErrorPlus> {
    let two = T::one() + T::one();

    // If output is F-order (or in general if iid stride is no more than sid_stride)
    if val.stride_of(nd::Axis(0)) <= val.stride_of(nd::Axis(1)) {
        let result_list = nd::Zip::from(val.axis_iter_mut(nd::Axis(1)))
            .and(stats.axis_iter_mut(nd::Axis(0)))
            .par_apply_collect(|mut col, mut stats_row| {
                _process_sid(
                    &mut col,
                    apply_in_place,
                    use_stats,
                    &mut stats_row,
                    &dist,
                    two,
                )
            });

        // Check the result list for errors
        result_list
            .iter()
            .par_bridge()
            .try_for_each(|x| (*x).clone())?;

        return Ok(());
    } else {
        //If C-order
        return _process_all_iids(val, apply_in_place, use_stats, stats, dist, two);
    }
}

fn find_factor<T: Default + Copy + Debug + Sync + Send + Float + ToPrimitive + FromPrimitive>(
    dist: &Dist,
    mean_s: T,
    std: T,
) -> Result<T, BedError> {
    if let Dist::Beta { a, b } = dist {
        // Try to create a beta dist
        let beta_dist = if let Ok(beta_dist) = Beta::new(*a, *b) {
            beta_dist
        } else {
            return Err(BedError::CannotCreateBetaDist(*a, *b));
        };

        // Try to an f64 maf
        let mut maf = if let Some(mean_u64) = mean_s.to_f64() {
            mean_u64 / 2.0
        } else {
            return Err(BedError::CannotConvertBetaToFromF64);
        };
        if maf > 0.5 {
            maf = 1.0 - maf;
        }

        // Try to put the maf in the beta dist
        return if let Some(b) = T::from_f64(beta_dist.pdf(maf)) {
            Ok(b)
        } else {
            Err(BedError::CannotConvertBetaToFromF64)
        };
    } else {
        return Ok(T::one() / std);
    }
}

fn _process_sid<T: Default + Copy + Debug + Sync + Send + Float + ToPrimitive + FromPrimitive>(
    col: &mut nd::ArrayViewMut1<'_, T>,
    apply_in_place: bool,
    use_stats: bool,
    stats_row: &mut nd::ArrayViewMut1<'_, T>,
    dist: &Dist,
    two: T,
) -> Result<(), BedError> {
    if !use_stats {
        let mut n_observed = T::zero();
        let mut sum_s = T::zero(); // the sum of a SNP over all observed individuals
        let mut sum2_s = T::zero(); // the sum of the squares of the SNP over all observed individuals

        for iid_i in 0..col.len() {
            let v = col[iid_i];
            if !v.is_nan() {
                sum_s = sum_s + v;
                sum2_s = sum2_s + v * v;
                n_observed = n_observed + T::one();
            }
        }
        if n_observed < T::one() {
            //LATER make it work (in some form) for n of 0
            return Err(BedError::NoIndividuals.into());
        }
        let mean_s = sum_s / n_observed; //compute the mean over observed individuals for the current SNP
        let mean2_s: T = sum2_s / n_observed; //compute the mean of the squared SNP

        if mean_s.is_nan()
            || (matches!(dist, Dist::Beta { a:_,b:_}) && ((mean_s > two) || (mean_s < T::zero())))
        {
            return Err(BedError::IllegalSnpMean.into());
        }

        let variance: T = mean2_s - mean_s * mean_s; //By the Cauchy Schwartz inequality this should always be positive

        let mut std = variance.sqrt();
        if std.is_nan() || std <= T::zero() {
            // All "SNPs" have the same value (aka SNC)
            std = T::infinity(); //SNCs are still meaning full in QQ plots because they should be thought of as SNPs without enough data.
        }

        stats_row[0] = mean_s;
        stats_row[1] = std;
    }

    if apply_in_place {
        {
            let mean_s = stats_row[0];
            let std = stats_row[1];
            let is_snc = std.is_infinite();

            let factor = find_factor(&dist, mean_s, std)?;

            for iid_i in 0..col.len() {
                //check for Missing (NAN) or SNC
                if col[iid_i].is_nan() || is_snc {
                    col[iid_i] = T::zero();
                } else {
                    col[iid_i] = (col[iid_i] - mean_s) * factor;
                }
            }
        }
    }
    return Ok(());
}

fn _process_all_iids<
    T: Default + Copy + Debug + Sync + Send + Float + ToPrimitive + FromPrimitive,
>(
    val: &mut nd::ArrayViewMut2<'_, T>,
    apply_in_place: bool,
    use_stats: bool,
    stats: &mut nd::ArrayViewMut2<'_, T>,
    dist: Dist,
    two: T,
) -> Result<(), BedErrorPlus> {
    let sid_count = val.dim().1;

    if !use_stats {
        // O(iid_count * sid_count)
        // Serial that respects C-order is 3-times faster than parallel that doesn't
        // So we parallelize the inner loop instead of the outer loop
        let mut n_observed_array = nd::Array1::<T>::zeros(sid_count);
        let mut sum_s_array = nd::Array1::<T>::zeros(sid_count); //the sum of a SNP over all observed individuals
        let mut sum2_s_array = nd::Array1::<T>::zeros(sid_count); //the sum of the squares of the SNP over all observed individuals
        for row in val.axis_iter(nd::Axis(0)) {
            nd::par_azip!((&v in row,
                n_observed_ptr in &mut n_observed_array,
                sum_s_ptr in &mut sum_s_array,
                sum2_s_ptr in &mut sum2_s_array
            )
                if !v.is_nan() {
                    *n_observed_ptr = *n_observed_ptr + T::one();
                    *sum_s_ptr = *sum_s_ptr + v;
                    *sum2_s_ptr = *sum2_s_ptr + v * v;
                }
            );
        }

        // O(sid_count)
        let mut result_list: Vec<Result<(), BedError>> = vec![Ok(()); sid_count];
        nd::par_azip!((mut stats_row in stats.axis_iter_mut(nd::Axis(0)),
                &n_observed in &n_observed_array,
                &sum_s in &sum_s_array,
                &sum2_s in &sum2_s_array,
                result_ptr in &mut result_list)
        {
            if n_observed < T::one() {
                *result_ptr = Err(BedError::NoIndividuals);
                return;
            }
            let mean_s = sum_s / n_observed; //compute the mean over observed individuals for the current SNP
            let mean2_s: T = sum2_s / n_observed; //compute the mean of the squared SNP

            if mean_s.is_nan()
                || (matches!(dist, Dist::Beta { a:_, b:_ }) && ((mean_s > two) || (mean_s < T::zero())))
            {
                *result_ptr = Err(BedError::IllegalSnpMean);
                return;
            }

            let variance: T = mean2_s - mean_s * mean_s; //By the Cauchy Schwartz inequality this should always be positive
            let mut std = variance.sqrt();
            if std.is_nan() || std <= T::zero() {
                // All "SNPs" have the same value (aka SNC)
                std = T::infinity(); //SNCs are still meaning full in QQ plots because they should be thought of as SNPs without enough data.
            }
            stats_row[0] = mean_s;
            stats_row[1] = std;
        });
        // Check the result list for errors
        result_list.par_iter().try_for_each(|x| (*x).clone())?;
    }

    if apply_in_place {
        // O(sid_count)
        let mut factor_array = nd::Array1::<T>::zeros(stats.dim().0);

        stats
            .axis_iter_mut(nd::Axis(0))
            .zip(&mut factor_array)
            .par_bridge()
            .try_for_each(|(stats_row, factor_ptr)| {
                match find_factor(&dist, stats_row[0], stats_row[1]) {
                    Err(e) => Err(e),
                    Ok(factor) => {
                        *factor_ptr = factor;
                        Ok(())
                    }
                }
            })?;

        // O(iid_count * sid_count)
        nd::par_azip!((mut row in val.axis_iter_mut(nd::Axis(0)))
        {
            for sid_i in 0..row.len() {
                //check for Missing (NAN) or SNC
                if row[sid_i].is_nan() || stats[(sid_i, 1)].is_infinite() {
                    row[sid_i] = T::zero();
                } else {
                    row[sid_i] = (row[sid_i] - stats[(sid_i, 0)]) * factor_array[sid_i];
                }
            }
        });
    }
    return Ok(());
}

pub fn create_pool(num_threads: usize) -> Result<rayon::ThreadPool, BedErrorPlus> {
    match rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
    {
        Err(e) => Err(e.into()),
        Ok(pool) => Ok(pool),
    }
}

mod python_module;
mod tests;
