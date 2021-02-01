// Inspired by C++ version by Chris Widmer and Carl Kadie

use core::fmt::Debug;
use ndarray as nd;
use ndarray::ShapeBuilder;
use num_traits::{Float, FromPrimitive, ToPrimitive};
use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;
use statrs::distribution::{Beta, Continuous};
use std::{
    fs::File,
    io::{BufRead, BufWriter, Read, Write},
    // cmk mem::{self, size_of},
    num::TryFromIntError,
    // cmk time::Instant,
    vec,
};
use std::{io::SeekFrom, path::PathBuf};
use std::{
    io::{BufReader, Seek},
    path::Path,
};
use thiserror::Error; // !!!cmk

const BED_FILE_MAGIC1: u8 = 0x6C; // 0b01101100 or 'l' (lowercase 'L')
const BED_FILE_MAGIC2: u8 = 0x1B; // 0b00011011 or <esc>
const CB_HEADER: u64 = 3;

// About ndarray
//  https://docs.rs/ndarray/0.14.0/ndarray/parallel/index.html
//  https://rust-lang-nursery.github.io/rust-cookbook/concurrency/parallel.html
//  https://github.com/rust-ndarray/ndarray/blob/master/README-quick-start.md
//  https://datacrayon.com/posts/programming/rust-notebooks/multidimensional-arrays-and-operations-with-ndarray
//  https://docs.rs/ndarray/0.14.0/ndarray/doc/ndarray_for_numpy_users/index.html
//  https://docs.rs/ndarray-npy
//  https://rust-lang-nursery.github.io/rust-cookbook/science/mathematics/linear_algebra.html

/// BedError enumerates all possible errors returned by this library.
/// Based on https://nick.groenen.me/posts/rust-error-handling/#the-library-error-type
#[derive(Error, Debug)]
pub enum BedError {
    #[error("Ill-formed BED file. BED file header is incorrect.")]
    IllFormed,

    #[error("Ill-formed BED file. BED file header is incorrect. Expected mode to be 0 or 1.")]
    BadMode,

    #[error("Attempt to write illegal value to BED file. Only 0,1,2,missing allowed.")]
    BadValue,

    #[error("No individual observed for the SNP.")]
    NoIndividuals,

    #[error("Illegal SNP mean.")]
    IllegalSnpMean,

    #[error("Output matrix dimensions doesn't match the length of the indexes.")]
    SubsetMismatch,

    #[error(transparent)]
    ConversionError(#[from] TryFromIntError),

    /// Represents all other cases of `std::io::Error`.
    #[error(transparent)]
    IOError(#[from] std::io::Error),
}

// !!!cmk "no_net_alloc"???
fn read_no_alloc<TOut: Copy + Default + From<i8> + Debug + Sync + Send>(
    filename: &str,
    iid_count: usize,
    sid_count: usize,
    count_a1: bool,
    iid_index: &[usize],
    sid_index: &[usize],
    missing_value: TOut,
    val: &mut nd::ArrayViewMut2<'_, TOut>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.
) -> Result<(), BedError> {
    let reader2 = BufReader::new(File::open(filename)?);
    let mut s = reader2.bytes();
    let rd1 = s.next().ok_or(BedError::IllFormed)??;
    let rd2 = s.next().ok_or(BedError::IllFormed)??;
    if (BED_FILE_MAGIC1 != rd1) || (BED_FILE_MAGIC2 != rd2) {
        return Err(BedError::IllFormed);
    }
    let rd3 = s.next().ok_or(BedError::IllFormed)??; // !!! cmk learn how to manipulate result codes
    match rd3 {
        0 => {
            //This option -- Sample major -- is untested because it is (almost?) never used.
            // cmk return Err(BedError::BadMode);
            let mut val_t = val.view_mut().reversed_axes();
            return _internal_read_no_alloc(
                filename,
                sid_count,
                count_a1,
                sid_index,
                iid_index,
                missing_value,
                &mut val_t,
            );
        }
        1 => {
            return _internal_read_no_alloc(
                filename,
                iid_count,
                count_a1,
                iid_index,
                sid_index,
                missing_value,
                val,
            );
        }
        _ => {
            return Err(BedError::BadMode);
        }
    }
}

fn _internal_read_no_alloc<TOut: Copy + Default + From<i8> + Debug + Sync + Send>(
    filename: &str,
    iid_count: usize,
    count_a1: bool,
    iid_index: &[usize],
    sid_index: &[usize],
    missing_value: TOut,
    val: &mut nd::ArrayViewMut2<'_, TOut>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.
) -> Result<(), BedError> {
    {}

    let iid_count_div4 = (iid_count + 3) / 4; // 4 genotypes per byte so round up
    let iid_count_div4_u64 = iid_count_div4 as u64; // !!!raise error
    let iid_out_count = iid_index.len();
    let sid_out_count = sid_index.len();

    let from_two_bits_to_value = set_up_two_bits_to_value(count_a1, missing_value);

    let mut reader = BufReader::new(File::open(filename)?);

    // See https://morestina.net/blog/1432/parallel-stream-processing-with-rayon
    // Possible optimization: We could try to read only the iid info needed
    // Possible optimization: We could read snp in their input order instead of their output order
    let _ = (0..sid_out_count)
        // Read all the iid info for one snp from the disk
        .map(|snp_out_i| {
            // if snp_out_i % 100 == 0 {
            //     // !!!cmk
            //     println!("{} of {}", snp_out_i, sid_out_count);
            // }
            let snp_in_i = sid_index[snp_out_i];
            let mut bytes_vector: Vec<u8> = vec![0; iid_count_div4];
            let pos: u64 = (snp_in_i as u64) * iid_count_div4_u64 + CB_HEADER; // !!!raise error
            reader.seek(SeekFrom::Start(pos)).unwrap(); // !!! cmk think about what it means to have unwrap here not ?
            reader.read_exact(&mut bytes_vector).unwrap();
            bytes_vector
        })
        // Zip in the column of the output array
        .zip(val.axis_iter_mut(nd::Axis(1)))
        // In parallel, decompress the iid info and put it in its column
        .par_bridge()
        .for_each(|(bytes_vector, mut col)| {
            for iid_out_i in 0..iid_out_count {
                // Possible optimization: We could pre-compute the conversion, the division, the mod, and the multiply*2
                let iid_in_i = iid_index[iid_out_i];
                let i_div_4 = iid_in_i / 4;
                let i_mod_4 = iid_in_i % 4;
                let genotype_byte: u8 = (bytes_vector[i_div_4] >> (i_mod_4 * 2)) & 0x03;
                col[iid_out_i] = from_two_bits_to_value[genotype_byte as usize];
            }
        });

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

// make count_a1 optional
pub fn read_with_indexes<TOut: From<i8> + Default + Copy + Debug + Sync + Send>(
    filename: &str,
    iid_index: &[usize],
    sid_index: &[usize],
    output_is_order_f: bool,
    count_a1: bool,
    missing_value: TOut,
) -> nd::Array2<TOut> {
    let path = Path::new(filename);
    let iid_count = count_lines(path.with_extension("fam"));
    let sid_count = count_lines(path.with_extension("bim"));

    let shape = ShapeBuilder::set_f((iid_index.len(), sid_index.len()), output_is_order_f);
    let mut val = nd::Array2::<TOut>::default(shape);

    let _ = read_no_alloc(
        filename, // !!!cmk can you call with keywords?
        iid_count,
        sid_count,
        count_a1,
        iid_index,
        sid_index,
        missing_value,
        &mut val.view_mut(),
    );

    return val;
}

// !!!cmk put in [No ignore] thing to force results to be used
// !!!cmk return good error messages from threads

pub fn read<TOut: From<i8> + Default + Copy + Debug + Sync + Send>(
    filename: &str,
    output_is_order_f: bool,
    count_a1: bool,
    missing_value: TOut,
) -> nd::Array2<TOut> {
    let path = Path::new(filename);
    let iid_count = count_lines(path.with_extension("fam"));
    let sid_count = count_lines(path.with_extension("bim"));

    let iid_index: Vec<usize> = (0..iid_count).collect();
    let sid_index: Vec<usize> = (0..sid_count).collect();

    let mut val = nd::Array2::<TOut>::default(ShapeBuilder::set_f(
        (iid_count as usize, sid_count as usize),
        output_is_order_f,
    ));

    let _ = read_no_alloc(
        filename, // !!!cmk can you call with keywords?
        iid_count,
        sid_count,
        count_a1,
        &iid_index,
        &sid_index,
        missing_value,
        &mut val.view_mut(),
    );

    return val;
}

// !!!cmk add thread control
pub fn write<T: From<i8> + Default + Copy + Debug + Sync + Send + PartialEq>(
    filename: &str,
    val: &nd::ArrayView2<'_, T>,
    count_a1: bool,
    missing: (bool, T),
) -> Result<(), BedError> {
    let mut writer = BufWriter::new(File::create(filename)?);
    writer.write_all(&[BED_FILE_MAGIC1, BED_FILE_MAGIC2, 0x01])?;

    let zero_code = if count_a1 { 3u8 } else { 0u8 };
    let two_code = if count_a1 { 0u8 } else { 3u8 };

    let homozygous_primary_allele = T::from(0); // Major Allele
    let heterozygous_allele = T::from(1);
    let homozygous_secondary_allele = T::from(2); // Minor Allele

    let iid_count = val.shape()[0];
    let iid_count_div4 = (iid_count + 3) / 4; // 4 genotypes per byte so round up
                                              // !!!cmk let iid_count_div4_u64 = iid_count_div4 as u64; // !!!raise error

    let (use_nan, other_missing_value) = missing;

    for column in val.axis_iter(nd::Axis(1)) {
        let mut bytes_vector: Vec<u8> = vec![0; iid_count_div4]; // inits to 0
        for (iid_i, &v0) in column.iter().enumerate() {
            let genotype_byte;
            if v0 == homozygous_primary_allele {
                genotype_byte = zero_code;
            } else if v0 == heterozygous_allele {
                genotype_byte = 2u8;
            } else if v0 == homozygous_secondary_allele {
                genotype_byte = two_code;
            } else if (use_nan && v0 != v0) || (!use_nan && v0 == other_missing_value) {
                genotype_byte = 1;
            } else {
                // println!("cmk bad value '{:?}'", v0);
                return Err(BedError::BadValue);
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

// !!!cmk make private
fn count_lines(path_buf: PathBuf) -> usize {
    let reader = BufReader::new(File::open(path_buf).expect("Cannot open file"));
    let count = reader.lines().count();
    return count;
}
// !!!cmk where are the return codes? Here and elsewhere
pub fn counts(filename: &str) -> (usize, usize) {
    let path = Path::new(filename);
    let iid_count = count_lines(path.with_extension("fam"));
    let sid_count = count_lines(path.with_extension("bim"));
    return (iid_count, sid_count);
}

// // #!!!cmk
// fn subset_64_64_f(
//     in_val: &nd::ArrayView2<'_, f64>,
//     iid_index: &[usize],
//     sid_index: &[usize],
//     out_val: &mut nd::ArrayViewMut2<'_, f64>,
// ) -> Result<(), BedError> {
//     let out_iid_count = iid_index.len();
//     let out_sid_count = sid_index.len();
//     if out_iid_count != out_val.shape()[0] || out_sid_count != out_val.shape()[1] {
//         return Err(BedError::SubsetMismatch);
//     }

//     assert!(out_val.stride_of(nd::Axis(0)) <= out_val.stride_of(nd::Axis(1)));
//     nd::par_azip!((mut col in out_val.axis_iter_mut(nd::Axis(1)),
//                 ptr_sid_i_in in sid_index)
//         {
//             let sid_i_in = *ptr_sid_i_in;
//             for (iid_i_out, ptr_iid_i_in) in iid_index.iter().enumerate() {
//                 let iid_i_in = *ptr_iid_i_in;
//                 col[(iid_i_out)] = in_val[(iid_i_in, sid_i_in)];
//             }
//     }
//     );
//     return Ok(());
// }

pub fn matrix_subset_no_alloc<
    TIn: Copy + Default + Debug + Sync + Send + Sized,
    TOut: Copy + Default + Debug + Sync + Send + From<TIn>,
>(
    in_val: &nd::ArrayView3<'_, TIn>,
    iid_index: &[usize],
    sid_index: &[usize],
    out_val: &mut nd::ArrayViewMut3<'_, TOut>,
) -> Result<(), BedError> {
    let out_iid_count = iid_index.len();
    let out_sid_count = sid_index.len();
    if out_iid_count != out_val.shape()[0] || out_sid_count != out_val.shape()[1] {
        return Err(BedError::SubsetMismatch);
    }

    let did_count = in_val.shape()[2];

    // If output is F-order (or in general if iid stride is no more than sid_stride)
    if out_val.stride_of(nd::Axis(0)) <= out_val.stride_of(nd::Axis(1)) {
        nd::par_azip!((mut out_col in out_val.axis_iter_mut(nd::Axis(1)),
                    sid_i_in_ptr in sid_index) {
            let in_col = in_val.index_axis(nd::Axis(1), *sid_i_in_ptr);
            for did_i in 0..did_count
            {
                for (iid_i_out, iid_i_in_ptr) in iid_index.iter().enumerate() {
                    out_col[(iid_i_out,did_i)] = in_col[(*iid_i_in_ptr,did_i)].into();
                }
            }
        });
    } else {
        //If output is C-order, transpose input and output and recurse
        let val_in_t = in_val.view().permuted_axes([1, 0, 2]);
        let mut val_out_t = out_val.view_mut().permuted_axes([1, 0, 2]); // !!!cmk out_val or val_out -- be consistent
        return matrix_subset_no_alloc(&val_in_t, &sid_index, &iid_index, &mut val_out_t);
    }

    return Ok(());
}

pub fn impute_and_zero_mean_snps<
    T: Default + Copy + Debug + Sync + Send + Float + ToPrimitive + FromPrimitive,
>(
    val: &mut nd::ArrayViewMut2<'_, T>,
    beta_not_unit_variance: bool,
    beta_a: f64,
    beta_b: f64,
    apply_in_place: bool,
    use_stats: bool,
    stats: &mut nd::ArrayViewMut2<'_, T>,
) -> Result<(), BedError> {
    let two = T::one() + T::one();

    // If output is F-order (or in general if iid stride is no more than sid_stride)
    if val.stride_of(nd::Axis(0)) <= val.stride_of(nd::Axis(1)) {
        nd::par_azip!((mut col in val.axis_iter_mut(nd::Axis(1)),
                   mut stats_row in stats.axis_iter_mut(nd::Axis(0))){
        _process_sid(
            &mut col,
            apply_in_place,
            use_stats,
            &mut stats_row,
            beta_a,
            beta_b,
            beta_not_unit_variance,
            two).unwrap();
        });
        return Ok(());
    } else {
        //If C-order
        return _process_all_iids(
            val,
            apply_in_place,
            use_stats,
            stats,
            beta_not_unit_variance,
            beta_a,
            beta_b,
            two,
        );
    }
}

fn find_factor<T: Default + Copy + Debug + Sync + Send + Float + ToPrimitive + FromPrimitive>(
    beta_not_unit_variance: bool,
    beta_a: f64,
    beta_b: f64,
    mean_s: T,
    std: T,
) -> T {
    if beta_not_unit_variance {
        let beta_dist = Beta::new(beta_a, beta_b).unwrap();
        let mut maf = mean_s.to_f64().unwrap() / 2.0;
        if maf > 0.5 {
            maf = 1.0 - maf;
        }
        return T::from_f64(beta_dist.pdf(maf)).unwrap();
    } else {
        return T::one() / std;
    }
}

fn _process_sid<T: Default + Copy + Debug + Sync + Send + Float + ToPrimitive + FromPrimitive>(
    col: &mut nd::ArrayViewMut1<'_, T>,
    apply_in_place: bool,
    use_stats: bool,
    stats_row: &mut nd::ArrayViewMut1<'_, T>,
    beta_a: f64,
    beta_b: f64,
    beta_not_unit_variance: bool,
    two: T,
) -> Result<(), BedError> {
    if !use_stats {
        let mut n_observed = T::zero();
        let mut sum_s = T::zero(); //the sum of a SNP over all observed individuals
        let mut sum2_s = T::zero(); //the sum of the squares of the SNP over all observed individuals

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
            return Err(BedError::NoIndividuals);
        }
        let mean_s = sum_s / n_observed; //compute the mean over observed individuals for the current SNP
        let mean2_s: T = sum2_s / n_observed; //compute the mean of the squared SNP

        if mean_s.is_nan() || (beta_not_unit_variance && ((mean_s > two) || (mean_s < T::zero()))) {
            return Err(BedError::IllegalSnpMean);
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

            let factor = find_factor(beta_not_unit_variance, beta_a, beta_b, mean_s, std);

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
    beta_not_unit_variance: bool,
    beta_a: f64,
    beta_b: f64,
    two: T,
) -> Result<(), BedError> {
    let sid_count = val.shape()[1];

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
                {
                if !v.is_nan() {
                    *n_observed_ptr = *n_observed_ptr + T::one();
                    *sum_s_ptr = *sum_s_ptr + v;
                    *sum2_s_ptr = *sum2_s_ptr + v * v;
                }
            }
            );
        }

        // O(sid_count)
        nd::par_azip!((mut stats_row in stats.axis_iter_mut(nd::Axis(0)),
                &n_observed in &n_observed_array,
                &sum_s in &sum_s_array,
                &sum2_s in &sum2_s_array)
        {
            if n_observed < T::one() {
                //LATER make it work (in some form) for n of 0
                panic!("no individuals");
                // !!! cmk return Err(BedError::NoIndividuals);
            }
            let mean_s = sum_s / n_observed; //compute the mean over observed individuals for the current SNP
            let mean2_s: T = sum2_s / n_observed; //compute the mean of the squared SNP

            if mean_s.is_nan()
                || (beta_not_unit_variance && ((mean_s > two) || (mean_s < T::zero())))
            {
                panic!("IllegalSnpMean")
                // !!!cmk return Err(BedError::IllegalSnpMean);
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
    }

    if apply_in_place {
        // O(sid_count)
        let mut factor_array = nd::Array1::<T>::zeros(stats.shape()[0]);
        nd::par_azip!((factor_ptr in &mut factor_array, stats_row in stats.axis_iter_mut(nd::Axis(0)))
            {
                *factor_ptr = find_factor(
                    beta_not_unit_variance,
                    beta_a,
                    beta_b,
                    stats_row[0],
                    stats_row[1],
                );
            }
        );

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

pub fn create_pool(num_threads: usize) -> rayon::ThreadPool {
    return rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
}

mod python_module;
mod tests;

// !!!cmk add default values
// !!!cmk add thread control
// !!!cmk Don't forget the Bed writer and log gamma, etc
