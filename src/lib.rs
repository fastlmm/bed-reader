#![warn(missing_docs)]
// Inspired by C++ version by Chris Widmer and Carl Kadie

// See: https://towardsdatascience.com/nine-rules-for-writing-python-extensions-in-rust-d35ea3a4ec29?sk=f8d808d5f414154fdb811e4137011437
// for an article on how this project uses Rust to create a Python extension.

// For Rust API tips see https://rust-lang.github.io/api-guidelines/necessities.html
#![doc = include_str!("../README-rust.md")]
//! ## Main Functions
//!
//! | Function | Description |
//! | -------- | ----------- |
//! | [`Bed::new`](struct.Bed.html#method.new) or [`Bed::builder`](struct.Bed.html#method.builder) | Open a PLINK .bed file for reading genotype data and metadata. |
//! | [`ReadOptions::builder`](struct.ReadOptions.html#method.builder) | Read genotype data. Supports indexing and options. |
//! | [`WriteOptions::builder`](struct.WriteOptions.html#method.builder) | Write values to a file in PLINK .bed format. Supports metadata and options. |
//!
//! ### `Bed` Metadata Methods
//!
//! After using [`Bed::new`](struct.Bed.html#method.new) or [`Bed::builder`](struct.Bed.html#method.builder) to open a PLINK .bed file for reading, use
//! these methods to see metadata.
//!
//! | Method | Description |
//! | -------- | ----------- |
//! | [`iid_count`](struct.Bed.html#method.iid_count) | Number of individuals (samples) |
//! | [`sid_count`](struct.Bed.html#method.sid_count) | Number of SNPs (variants) |
//! | [`dim`](struct.Bed.html#method.dim) | Number of individuals and SNPs |
//! | [`fid`](struct.Bed.html#method.fid) | Family id of each of individual (sample) |
//! | [`iid`](struct.Bed.html#method.iid) | Individual id of each of individual (sample) |
//! | [`father`](struct.Bed.html#method.father) | Father id of each of individual (sample) |
//! | [`mother`](struct.Bed.html#method.mother) | Mother id of each of individual (sample) |
//! | [`sex`](struct.Bed.html#method.sex) | Sex of each individual (sample) |
//! | [`pheno`](struct.Bed.html#method.pheno) | A phenotype for each individual (seldom used) |
//! | [`chromosome`](struct.Bed.html#method.chromosome) | Chromosome of each SNP (variant) |
//! | [`sid`](struct.Bed.html#method.sid) | SNP Id of each SNP (variant) |
//! | [`cm_position`](struct.Bed.html#method.cm_position) | Centimorgan position of each SNP (variant) |
//! | [`bp_position`](struct.Bed.html#method.bp_position) | Base-pair position of each SNP (variant) |
//! | [`allele_1`](struct.Bed.html#method.allele_1) | First allele of each SNP (variant) |
//! | [`allele_2`](struct.Bed.html#method.allele_2) | Second allele of each SNP (variant) |
//! | [`metadata`](struct.Bed.html#method.metadata) | All the metadata returned as a [`struct.Metadata`](struct.Metadata.html) |
//!
//! ### `ReadOptions`
//!
//! When using [`ReadOptions::builder`](struct.ReadOptions.html#method.builder) to read genotype data, use these options to
//! specify a desired numeric type,
//! which individuals (samples) to read, which SNPs (variants) to read, etc.
//!
//! | Option | Description |
//! | -------- | ----------- |
//! | [`i8`](struct.ReadOptionsBuilder.html#method.i8) | Read values as i8 |
//! | [`f32`](struct.ReadOptionsBuilder.html#method.f32) | Read values as f32 |
//! | [`f64`](struct.ReadOptionsBuilder.html#method.f64) | Read values as f64 |
//! | [`iid_index`](struct.ReadOptionsBuilder.html#method.iid_index) | Index of individuals (samples) to read (defaults to all)|
//! | [`sid_index`](struct.ReadOptionsBuilder.html#method.sid_index) | Index of SNPs (variants) to read (defaults to all) |
//! | [`f`](struct.ReadOptionsBuilder.html#method.f) | Order of the output array, Fortran-style (default) |
//! | [`c`](struct.ReadOptionsBuilder.html#method.c) | Order of the output array, C-style |
//! | [`is_f`](struct.ReadOptionsBuilder.html#method.is_f) | Is order of the output array Fortran-style? (defaults to true)|
//! | [`missing_value`](struct.ReadOptionsBuilder.html#method.missing_value) | Value to use for missing values (defaults to -127 or NaN) |
//! | [`count_a1`](struct.ReadOptionsBuilder.html#method.count_a1) | Count the number allele 1 (default) |
//! | [`count_a2`](struct.ReadOptionsBuilder.html#method.count_a2) | Count the number allele 2 |
//! | [`is_a1_counted`](struct.ReadOptionsBuilder.html#method.is_a1_counted) | Is allele 1 counted? (defaults to true) |
//! | [`num_threads`](struct.ReadOptionsBuilder.html#method.num_threads) | Number of threads to use (defaults to all processors) |
//!
//! ### [`Index`](enum.Index.html) Expressions
//!
//! Select which individuals (samples) and SNPs (variants) to read by using these
//! [`iid_index`](struct.ReadOptionsBuilder.html#method.iid_index) and/or
//! [`sid_index`](struct.ReadOptionsBuilder.html#method.sid_index) expressions.
//!
//! | Example | Type | Description |
//! | -------- | --- | ----------- |
//! | nothing | `()` | All |
//! | `2` | `isize` | Index position 2 |
//! | `-1` | `isize` | Last index position |
//! | `vec![0, 10, -2]` | `Vec<isize>` | Index positions 0, 10, and 2nd from last |
//! | `[0, 10, -2]` | `[isize]` and `[isize;n]` | Index positions 0, 10, and 2nd from last |
//! | `ndarray::array![0, 10, -2]` | `ndarray::Array1<isize>` | Index positions 0, 10, and 2nd from last |
//! | `10..20` | `Range<usize>` | Index positions 10 (inclusive) to 20 (exclusive). *Note: Rust ranges don't support negatives* |
//! | `..=19` | `RangeInclusive<usize>` | Index positions 0 (inclusive) to 19 (inclusive). *Note: Rust ranges don't support negatives* |
//! | *any Rust ranges* | `Range*<usize>` | *Note: Rust ranges don't support negatives* |
//! | `s![10..20;2]` | `ndarray::SliceInfo1` | Index positions 10 (inclusive) to 20 (exclusive) in steps of 2 |
//! | `s![-20..-10;-2]` | `ndarray::SliceInfo1` | 10th from last (exclusive) to 20th from last (inclusive), in steps of -2 |
//! | `vec![true, false, true]` | `Vec<bool>`| Index positions 0 and 2. |
//! | `[true, false, true]` | `[bool]` and `[bool;n]`| Index positions 0 and 2.|
//! | `ndarray::array![true, false, true]` | `ndarray::Array1<bool>`| Index positions 0 and 2.|
//!
//! ### Environment Variables
//!
//! * `BED_READER_NUM_THREADS`
//! * `NUM_THREADS`
//!
//! If [`ReadOptionsBuilder::num_threads`](struct.ReadOptionsBuilder.html#method.num_threads)
//! or [`WriteOptionsBuilder::num_threads`](struct.WriteOptionsBuilder.html#method.num_threads) is not specified,
//! the number of threads to use is determined by these environment variable (in order of priority):
//! If neither of these environment variables are set, all processors are used.
//!
//! * `BED_READER_DATA_DIR`
//!
//! Any requested sample file will be downloaded to this directory. If the environment variable is not set,
//! a cache folder, appropriate to the OS, will be used.

mod python_module;
mod tests;
use anyinput::anyinput;
use bed_cloud::BedCloud;
use core::fmt::Debug;
use derive_builder::{Builder, UninitializedFieldError};
use fetch_data::{FetchData, FetchDataError};
use futures_util::pin_mut;
use futures_util::StreamExt;
use nd::ShapeBuilder;
use ndarray as nd;
use object_store::delimited::newline_delimited_stream;
use object_store::path::Path as StorePath;
use object_store::ObjectStore;
use std::cmp::Ordering;
use std::collections::HashSet;
use std::convert::TryFrom;
use std::fs::{self};
use std::io::Write;
use std::ops::{
    Bound, Deref, Range, RangeBounds, RangeFrom, RangeInclusive, RangeTo, RangeToInclusive,
};
use std::rc::Rc;
use std::str::Utf8Error;
use std::{
    env,
    fs::File,
    io::{BufRead, BufReader, BufWriter},
    ops::RangeFull,
    path::{Path, PathBuf},
};

use byteorder::{LittleEndian, ReadBytesExt};
use dpc_pariter::{scope, IteratorExt};
use num_traits::{abs, Float, FromPrimitive, Signed, ToPrimitive};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rayon::{iter::ParallelBridge, ThreadPoolBuildError};
use statrs::distribution::{Beta, Continuous};
use std::io::Read;
use std::io::Seek;
use std::io::SeekFrom;
use std::num::{ParseFloatError, ParseIntError};
use std::ops::AddAssign;
use std::ops::{Div, Sub};
use thiserror::Error;
use tokio::task::JoinError;
/// cmk docks
pub mod bed_cloud;

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

/// All possible errors returned by this library and the libraries it depends on.
// Based on `<https://nick.groenen.me/posts/rust-error-handling/#the-library-error-type>`
#[derive(Error, Debug)]
pub enum BedErrorPlus {
    #[allow(missing_docs)]
    #[error(transparent)]
    BedError(#[from] BedError),

    #[allow(missing_docs)]
    #[error(transparent)]
    IOError(#[from] std::io::Error),

    #[allow(missing_docs)]
    #[error(transparent)]
    ThreadPoolError(#[from] ThreadPoolBuildError),

    #[allow(missing_docs)]
    #[error(transparent)]
    ParseIntError(#[from] ParseIntError),

    #[allow(missing_docs)]
    #[error(transparent)]
    UninitializedFieldError(#[from] ::derive_builder::UninitializedFieldError),

    #[allow(missing_docs)]
    #[error(transparent)]
    ParseFloatError(#[from] ParseFloatError),

    #[allow(missing_docs)]
    #[error(transparent)]
    FetchData(#[from] FetchDataError),

    #[allow(missing_docs)]
    #[error(transparent)]
    ObjectStoreError(#[from] object_store::Error),

    #[allow(missing_docs)]
    #[error(transparent)]
    ObjectStorePathError(#[from] object_store::path::Error),

    #[allow(missing_docs)]
    #[error(transparent)]
    JoinError(#[from] JoinError),

    #[allow(missing_docs)]
    #[error(transparent)]
    Utf8Error(#[from] Utf8Error),
}
// https://docs.rs/thiserror/1.0.23/thiserror/

/// All errors specific to this library.
#[derive(Error, Debug, Clone)]
pub enum BedError {
    #[allow(missing_docs)]
    #[error("Ill-formed BED file. BED file header is incorrect or length is wrong. '{0}'")]
    IllFormed(String),

    #[allow(missing_docs)]
    #[error(
        "Ill-formed BED file. BED file header is incorrect. Expected mode to be 0 or 1. '{0}'"
    )]
    BadMode(String),

    #[allow(missing_docs)]
    #[error("Attempt to write illegal value to BED file. Only 0,1,2,missing allowed. '{0}'")]
    BadValue(String),

    #[allow(missing_docs)]
    #[error("Multithreading resulted in panic(s)")]
    PanickedThread(),

    #[allow(missing_docs)]
    #[error("No individual observed for the SNP.")]
    NoIndividuals,

    #[allow(missing_docs)]
    #[error("Illegal SNP mean.")]
    IllegalSnpMean,

    #[allow(missing_docs)]
    #[error("Index to individual larger than the number of individuals. (Index value {0})")]
    IidIndexTooBig(isize),

    #[allow(missing_docs)]
    #[error("Index to SNP larger than the number of SNPs. (Index value {0})")]
    SidIndexTooBig(isize),

    #[allow(missing_docs)]
    #[error("Length of iid_index ({0}) and sid_index ({1}) must match dimensions of output array ({2},{3}).")]
    IndexMismatch(usize, usize, usize, usize),

    #[allow(missing_docs)]
    #[error("Indexes ({0},{1}) too big for files")]
    IndexesTooBigForFiles(usize, usize),

    #[allow(missing_docs)]
    #[error("Subset: length of iid_index ({0}) and sid_index ({1}) must match dimensions of output array ({2},{3}).")]
    SubsetMismatch(usize, usize, usize, usize),

    #[allow(missing_docs)]
    #[error("Cannot convert beta values to/from float 64")]
    CannotConvertBetaToFromF64,

    #[allow(missing_docs)]
    #[error("Cannot create Beta Dist with given parameters ({0},{1})")]
    CannotCreateBetaDist(f64, f64),

    #[allow(missing_docs)]
    #[error("Cannot use skipped metadata '{0}'")]
    CannotUseSkippedMetadata(String),

    #[allow(missing_docs)]
    #[error("Index starts at {0} but ends at {1}")]
    StartGreaterThanEnd(usize, usize),

    #[allow(missing_docs)]
    #[error("Step of zero not allowed")]
    StepZero,

    #[allow(missing_docs)]
    #[error("Index starts at {0} but count is {1}")]
    StartGreaterThanCount(usize, usize),

    #[allow(missing_docs)]
    #[error("Index ends at {0} but count is {1}")]
    EndGreaterThanCount(usize, usize),

    #[allow(missing_docs)]
    #[error("Adding new axis not allowed")]
    NewAxis,

    #[allow(missing_docs)]
    #[error("Expect 1-D NDArray SliceInfo")]
    NdSliceInfoNot1D,

    #[allow(missing_docs)]
    #[error("Expect {0} fields but find only {1} in '{2}'")]
    MetadataFieldCount(usize, usize, String),

    #[allow(missing_docs)]
    #[error("{0}_count values of {1} and {2} are inconsistent")]
    InconsistentCount(String, usize, usize),

    #[allow(missing_docs)]
    #[error("Expect bool arrays and vectors to be length {0}, not {1}")]
    BoolArrayVectorWrongLength(usize, usize),

    #[allow(missing_docs)]
    #[error("Expect ndarray of shape ({0}, {1}), but found shape ({2}, {3})")]
    InvalidShape(usize, usize, usize, usize),

    #[allow(missing_docs)]
    #[error("Can't write '{0}' metadata if some fields are None")]
    MetadataMissingForWrite(String),

    #[allow(missing_docs)]
    #[error("Unknown or bad sample file '{0}'")]
    UnknownOrBadSampleFile(String),

    #[allow(missing_docs)]
    #[error("The registry of sample files is invalid")]
    SampleRegistryProblem(),

    #[allow(missing_docs)]
    #[error("Samples construction failed with error: {0}")]
    SamplesConstructionFailed(String),

    #[allow(missing_docs)]
    #[error("Downloaded sample file not seen: {0}")]
    DownloadedSampleFileNotSeen(String),

    #[allow(missing_docs)]
    #[error("Downloaded sample file has wrong hash: {0},expected: {1}, actual: {2}")]
    DownloadedSampleFileWrongHash(String, String, String),

    #[allow(missing_docs)]
    #[error("Cannot create cache directory")]
    CannotCreateCacheDir(),
}

// Trait alias

/// A trait alias, used internally, for the values of a .bed file, namely i8, f32, f64.
pub trait BedVal:
    Copy + Default + From<i8> + Debug + Sync + Send + Sync + Missing + PartialEq
{
}
impl<T> BedVal for T where
    T: Copy + Default + From<i8> + Debug + Sync + Send + Sync + Missing + PartialEq
{
}

fn create_pool(num_threads: usize) -> Result<rayon::ThreadPool, Box<BedErrorPlus>> {
    match rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
    {
        Err(e) => Err(Box::new(e.into())),
        Ok(pool) => Ok(pool),
    }
}

#[allow(clippy::too_many_arguments)]
#[anyinput]
fn read_no_alloc<TVal: BedVal>(
    path: AnyPath,
    iid_count: usize,
    sid_count: usize,
    is_a1_counted: bool,
    iid_index: &[isize],
    sid_index: &[isize],
    missing_value: TVal,
    num_threads: usize,
    val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.
) -> Result<(), Box<BedErrorPlus>> {
    create_pool(num_threads)?.install(|| {
        let (buf_reader, bytes_vector) = open_and_check(path)?;

        match bytes_vector[2] {
            0 => {
                // We swap 'iid' and 'sid' and then reverse the axes.
                let mut val_t = val.view_mut().reversed_axes();
                internal_read_no_alloc(
                    buf_reader,
                    path,
                    sid_count,
                    iid_count,
                    is_a1_counted,
                    sid_index,
                    iid_index,
                    missing_value,
                    &mut val_t,
                )
            }
            1 => internal_read_no_alloc(
                buf_reader,
                path,
                iid_count,
                sid_count,
                is_a1_counted,
                iid_index,
                sid_index,
                missing_value,
                val,
            ),
            _ => Err(Box::new(BedError::BadMode(path_ref_to_string(path)).into())),
        }
    })?;
    Ok(())
}

#[anyinput]
fn path_ref_to_string(path: AnyPath) -> String {
    PathBuf::from(path).display().to_string()
}

impl From<BedError> for Box<BedErrorPlus> {
    fn from(err: BedError) -> Self {
        Box::new(BedErrorPlus::BedError(err))
    }
}
impl From<std::io::Error> for Box<BedErrorPlus> {
    fn from(err: std::io::Error) -> Self {
        Box::new(BedErrorPlus::IOError(err))
    }
}
impl From<ThreadPoolBuildError> for Box<BedErrorPlus> {
    fn from(err: ThreadPoolBuildError) -> Self {
        Box::new(BedErrorPlus::ThreadPoolError(err))
    }
}
impl From<ParseIntError> for Box<BedErrorPlus> {
    fn from(err: ParseIntError) -> Self {
        Box::new(BedErrorPlus::ParseIntError(err))
    }
}

impl From<ParseFloatError> for Box<BedErrorPlus> {
    fn from(err: ParseFloatError) -> Self {
        Box::new(BedErrorPlus::ParseFloatError(err))
    }
}

impl From<::derive_builder::UninitializedFieldError> for Box<BedErrorPlus> {
    fn from(err: ::derive_builder::UninitializedFieldError) -> Self {
        Box::new(BedErrorPlus::UninitializedFieldError(err))
    }
}
impl From<FetchDataError> for Box<BedErrorPlus> {
    fn from(err: FetchDataError) -> Self {
        Box::new(BedErrorPlus::FetchData(err))
    }
}

#[anyinput]
fn open_and_check(path: AnyPath) -> Result<(BufReader<File>, Vec<u8>), Box<BedErrorPlus>> {
    let mut buf_reader = BufReader::new(File::open(path)?);
    let mut bytes_vector: Vec<u8> = vec![0; CB_HEADER_USIZE];
    buf_reader.read_exact(&mut bytes_vector)?;
    if (BED_FILE_MAGIC1 != bytes_vector[0]) || (BED_FILE_MAGIC2 != bytes_vector[1]) {
        return Err(Box::new(
            BedError::IllFormed(path_ref_to_string(path)).into(),
        ));
    }
    Ok((buf_reader, bytes_vector))
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

/// A trait alias, used internally, to provide default missing values for i8, f32, f64.
pub trait Missing {
    /// The default missing value for a type such as i8, f32, and f64.
    fn missing() -> Self;
}

impl Missing for f64 {
    fn missing() -> Self {
        f64::NAN
    }
}

impl Missing for f32 {
    fn missing() -> Self {
        f32::NAN
    }
}

impl Missing for i8 {
    fn missing() -> Self {
        -127i8
    }
}

// We make this generic instead of u64, so that we can test it via u8
fn try_div_4<T: Max + TryFrom<usize> + Sub<Output = T> + Div<Output = T> + Ord>(
    in_iid_count: usize,
    in_sid_count: usize,
    cb_header: T,
) -> Result<(usize, T), Box<BedErrorPlus>> {
    // 4 genotypes per byte so round up without overflow
    let in_iid_count_div4 = if in_iid_count > 0 {
        (in_iid_count - 1) / 4 + 1
    } else {
        0
    };
    let in_iid_count_div4_t = match T::try_from(in_iid_count_div4) {
        Ok(v) => v,
        Err(_) => {
            return Err(Box::new(
                BedError::IndexesTooBigForFiles(in_iid_count, in_sid_count).into(),
            ))
        }
    };
    let in_sid_count_t = match T::try_from(in_sid_count) {
        Ok(v) => v,
        Err(_) => {
            return Err(Box::new(
                BedError::IndexesTooBigForFiles(in_iid_count, in_sid_count).into(),
            ))
        }
    };

    let m: T = Max::max(); // Don't know how to move this into the next line.
    if in_sid_count > 0 && (m - cb_header) / in_sid_count_t < in_iid_count_div4_t {
        return Err(Box::new(
            BedError::IndexesTooBigForFiles(in_iid_count, in_sid_count).into(),
        ));
    }

    Ok((in_iid_count_div4, in_iid_count_div4_t))
}

#[allow(clippy::too_many_arguments)]
#[anyinput]
fn internal_read_no_alloc<TVal: BedVal>(
    mut buf_reader: BufReader<File>,
    path: AnyPath,
    in_iid_count: usize,
    in_sid_count: usize,
    is_a1_counted: bool,
    iid_index: &[isize],
    sid_index: &[isize],
    missing_value: TVal,
    out_val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.
) -> Result<(), Box<BedErrorPlus>> {
    // Check the file length

    let (in_iid_count_div4, in_iid_count_div4_u64) =
        try_div_4(in_iid_count, in_sid_count, CB_HEADER_U64)?;
    // "as" and math is safe because of early checks
    let file_len = buf_reader.seek(SeekFrom::End(0))?;
    let file_len2 = in_iid_count_div4_u64 * (in_sid_count as u64) + CB_HEADER_U64;
    if file_len != file_len2 {
        return Err(Box::new(
            BedError::IllFormed(path_ref_to_string(path)).into(),
        ));
    }

    // Check and precompute for each iid_index

    let (i_div_4_array, i_mod_4_times_2_array) =
        check_and_precompute_iid_index(in_iid_count, iid_index)?;

    // Check and compute work for each sid_index

    let from_two_bits_to_value = set_up_two_bits_to_value(is_a1_counted, missing_value);
    let lower_sid_count = -(in_sid_count as isize);
    let upper_sid_count: isize = (in_sid_count as isize) - 1;
    // See https://morestina.net/blog/1432/parallel-stream-processing-with-rayon
    // Possible optimization: We could try to read only the iid info needed
    // Possible optimization: We could read snp in their input order instead of their output order
    sid_index
        .iter()
        .map(|in_sid_i_signed| {
            // Turn signed sid_index into unsigned sid_index (or error)
            let in_sid_i = if (0..=upper_sid_count).contains(in_sid_i_signed) {
                *in_sid_i_signed as u64
            } else if (lower_sid_count..=-1).contains(in_sid_i_signed) {
                (in_sid_count - ((-in_sid_i_signed) as usize)) as u64
            } else {
                return Err(Box::new(BedErrorPlus::BedError(BedError::SidIndexTooBig(
                    *in_sid_i_signed,
                ))));
            };

            // Read the iid info for one snp from the disk
            let mut bytes_vector: Vec<u8> = vec![0; in_iid_count_div4];
            let pos: u64 = in_sid_i * in_iid_count_div4_u64 + CB_HEADER_U64; // "as" and math is safe because of early checks
            buf_reader.seek(SeekFrom::Start(pos))?;
            buf_reader.read_exact(&mut bytes_vector)?;
            Ok(bytes_vector)
        })
        // Zip in the column of the output array
        .zip(out_val.axis_iter_mut(nd::Axis(1)))
        // In parallel, decompress the iid info and put it in its column
        .par_bridge() // This seems faster that parallel zip
        .try_for_each(|(bytes_vector_result, mut col)| match bytes_vector_result {
            Err(e) => Err(e),
            Ok(bytes_vector) => {
                for out_iid_i in 0..iid_index.len() {
                    let i_div_4 = i_div_4_array[out_iid_i];
                    let i_mod_4_times_2 = i_mod_4_times_2_array[out_iid_i];
                    let genotype_byte: u8 = (bytes_vector[i_div_4] >> i_mod_4_times_2) & 0x03;
                    col[out_iid_i] = from_two_bits_to_value[genotype_byte as usize];
                }
                Ok(())
            }
        })?;

    Ok(())
}

type Array1Usize = nd::ArrayBase<nd::OwnedRepr<usize>, nd::Dim<[usize; 1]>>;
type Array1U8 = nd::ArrayBase<nd::OwnedRepr<u8>, nd::Dim<[usize; 1]>>;

fn check_and_precompute_iid_index(
    in_iid_count: usize,
    iid_index: &[isize],
) -> Result<(Array1Usize, Array1U8), Box<BedErrorPlus>> {
    let lower_iid_count = -(in_iid_count as isize);
    let upper_iid_count: isize = (in_iid_count as isize) - 1;
    let mut i_div_4_array = nd::Array1::<usize>::zeros(iid_index.len());
    let mut i_mod_4_times_2_array = nd::Array1::<u8>::zeros(iid_index.len());
    let mut result_list: Vec<Result<(), BedError>> = vec![Ok(()); iid_index.len()];
    nd::par_azip!((in_iid_i_signed in iid_index,
        i_div_4 in &mut i_div_4_array,
        i_mod_4_times_2 in &mut i_mod_4_times_2_array,
        result in &mut result_list
    )
    {
        let in_iid_i = if (0..=upper_iid_count).contains(in_iid_i_signed) {
            *result = Ok(());
            *in_iid_i_signed as usize
        } else if (lower_iid_count..=-1).contains(in_iid_i_signed) {
            *result = Ok(());
            in_iid_count - ((-in_iid_i_signed) as usize)
        } else {
            *result = Err(BedError::IidIndexTooBig(
                *in_iid_i_signed,
            ));
            0
        };

        *i_div_4 = in_iid_i / 4;
        *i_mod_4_times_2 = (in_iid_i % 4 * 2) as u8;
    });
    result_list
        .iter()
        .par_bridge()
        .try_for_each(|x| (*x).clone())?;
    Ok((i_div_4_array, i_mod_4_times_2_array))
}

fn set_up_two_bits_to_value<TVal: From<i8>>(count_a1: bool, missing_value: TVal) -> [TVal; 4] {
    let homozygous_primary_allele = TVal::from(0); // Major Allele
    let heterozygous_allele = TVal::from(1);
    let homozygous_secondary_allele = TVal::from(2); // Minor Allele

    if count_a1 {
        [
            homozygous_secondary_allele, // look-up 0
            missing_value,               // look-up 1
            heterozygous_allele,         // look-up 2
            homozygous_primary_allele,   // look-up 3
        ]
    } else {
        [
            homozygous_primary_allele,   // look-up 0
            missing_value,               // look-up 1
            heterozygous_allele,         // look-up 2
            homozygous_secondary_allele, // look-up 3
        ]
    }
}

// Thanks to Dawid for his dpc-pariter library that makes this function scale.
// https://dpc.pw/adding-parallelism-to-your-rust-iterators
#[anyinput]
fn write_val<S, TVal>(
    path: AnyPath,
    val: &nd::ArrayBase<S, nd::Ix2>,
    is_a1_counted: bool,
    missing: TVal,
    num_threads: usize,
) -> Result<(), Box<BedErrorPlus>>
where
    S: nd::Data<Elem = TVal>,
    TVal: BedVal,
{
    let (iid_count, sid_count) = val.dim();

    // 4 genotypes per byte so round up
    let (iid_count_div4, _) = try_div_4(iid_count, sid_count, CB_HEADER_U64)?;

    // We create and write to a file.
    // If there is an error, we will delete it.
    if let Err(e) = write_internal(
        path,
        iid_count_div4,
        val,
        is_a1_counted,
        missing,
        num_threads,
    ) {
        // Clean up the file
        let _ = fs::remove_file(path);
        Err(e)
    } else {
        Ok(())
    }
}

// https://www.reddit.com/r/rust/comments/mo4s8e/difference_between_reference_and_view_in_ndarray/
#[anyinput]
fn write_internal<S, TVal>(
    path: AnyPath,
    iid_count_div4: usize,
    //val: &nd::ArrayView2<'_, TVal>,
    val: &nd::ArrayBase<S, nd::Ix2>,
    is_a1_counted: bool,
    missing: TVal,
    num_threads: usize,
) -> Result<(), Box<BedErrorPlus>>
where
    S: nd::Data<Elem = TVal>,
    TVal: BedVal,
{
    let mut writer = BufWriter::new(File::create(path)?);
    writer.write_all(&[BED_FILE_MAGIC1, BED_FILE_MAGIC2, 0x01])?;

    #[allow(clippy::eq_op)]
    let use_nan = missing != missing; // generic NAN test
    let zero_code = if is_a1_counted { 3u8 } else { 0u8 };
    let two_code = if is_a1_counted { 0u8 } else { 3u8 };

    let homozygous_primary_allele = TVal::from(0); // Major Allele
    let heterozygous_allele = TVal::from(1);
    let homozygous_secondary_allele = TVal::from(2); // Minor Allele

    scope(|scope| {
        val.axis_iter(nd::Axis(1))
            .parallel_map_scoped(scope, {
                move |column| {
                    // Convert each column into a bytes_vector
                    let mut bytes_vector: Vec<u8> = vec![0; iid_count_div4]; // inits to 0
                    for (iid_i, &v0) in column.iter().enumerate() {
                        #[allow(clippy::eq_op)]
                        let genotype_byte = if v0 == homozygous_primary_allele {
                            zero_code
                        } else if v0 == heterozygous_allele {
                            2
                        } else if v0 == homozygous_secondary_allele {
                            two_code
                        //                    v0 !=v0 is generic NAN check
                        } else if (use_nan && v0 != v0) || (!use_nan && v0 == missing) {
                            1
                        } else {
                            return Err(BedError::BadValue(path_ref_to_string(path)));
                        };
                        // Possible optimization: We could pre-compute the conversion, the division, the mod, and the multiply*2
                        let i_div_4 = iid_i / 4;
                        let i_mod_4 = iid_i % 4;
                        bytes_vector[i_div_4] |= genotype_byte << (i_mod_4 * 2);
                    }
                    Ok(bytes_vector)
                }
            })
            .threads(num_threads)
            .try_for_each(|bytes_vector| {
                // Write the bytes vector, they must be in order.
                writer.write_all(&bytes_vector?)?;
                Ok(())
            })
    })
    .map_err(|_e| BedError::PanickedThread())?
}

#[anyinput]
fn count_lines(path: AnyPath) -> Result<usize, Box<BedErrorPlus>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let count = reader.lines().count();
    Ok(count)
}

#[allow(dead_code)]
enum Dist {
    Unit,
    Beta { a: f64, b: f64 },
}

#[allow(dead_code)]
fn impute_and_zero_mean_snps<
    T: Default + Copy + Debug + Sync + Send + Sync + Float + ToPrimitive + FromPrimitive,
>(
    val: &mut nd::ArrayViewMut2<'_, T>,
    dist: Dist,
    apply_in_place: bool,
    use_stats: bool,
    stats: &mut nd::ArrayViewMut2<'_, T>,
) -> Result<(), Box<BedErrorPlus>> {
    let two = T::one() + T::one();

    // If output is F-order (or in general if iid stride is no more than sid_stride)
    if val.stride_of(nd::Axis(0)) <= val.stride_of(nd::Axis(1)) {
        let result_list = nd::Zip::from(val.axis_iter_mut(nd::Axis(1)))
            .and(stats.axis_iter_mut(nd::Axis(0)))
            .par_map_collect(|mut col, mut stats_row| {
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

        Ok(())
    } else {
        //If C-order
        _process_all_iids(val, apply_in_place, use_stats, stats, dist, two)
    }
}

// Later move the other fast-lmm functions into their own package
#[allow(dead_code)]
fn find_factor<
    T: Default + Copy + Debug + Sync + Send + Sync + Float + ToPrimitive + FromPrimitive,
>(
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
        if let Some(b) = T::from_f64(beta_dist.pdf(maf)) {
            Ok(b)
        } else {
            Err(BedError::CannotConvertBetaToFromF64)
        }
    } else {
        Ok(T::one() / std)
    }
}

#[allow(dead_code)]
fn _process_sid<
    T: Default + Copy + Debug + Sync + Send + Sync + Float + ToPrimitive + FromPrimitive,
>(
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
            return Err(BedError::NoIndividuals);
        }
        let mean_s = sum_s / n_observed; //compute the mean over observed individuals for the current SNP
        let mean2_s: T = sum2_s / n_observed; //compute the mean of the squared SNP

        if mean_s.is_nan()
            || (matches!(dist, Dist::Beta { a: _, b: _ })
                && ((mean_s > two) || (mean_s < T::zero())))
        {
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

            let factor = find_factor(dist, mean_s, std)?;

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
    Ok(())
}

#[allow(dead_code)]
fn _process_all_iids<
    T: Default + Copy + Debug + Sync + Send + Sync + Float + ToPrimitive + FromPrimitive,
>(
    val: &mut nd::ArrayViewMut2<'_, T>,
    apply_in_place: bool,
    use_stats: bool,
    stats: &mut nd::ArrayViewMut2<'_, T>,
    dist: Dist,
    two: T,
) -> Result<(), Box<BedErrorPlus>> {
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
    Ok(())
}

#[allow(dead_code)]
#[anyinput]
fn file_b_less_aatbx(
    a_filename: AnyPath,
    offset: u64,
    iid_count: usize,
    b1: &mut nd::ArrayViewMut2<'_, f64>,
    aatb: &mut nd::ArrayViewMut2<'_, f64>,
    atb: &mut nd::ArrayViewMut2<'_, f64>,
    log_frequency: usize,
) -> Result<(), Box<BedErrorPlus>> {
    //speed idea from C++:
    //Are copies really needed?
    //is F, vc C order the best?
    //would bigger snp blocks be better

    let (a_sid_count, b_sid_count) = atb.dim();
    if log_frequency > 0 {
        println!(
            "file_b_less_aatbx: iid_count={}, {}x{} output",
            iid_count, a_sid_count, b_sid_count
        );
    };

    // Open the file and move to the starting sid
    let mut buf_reader = BufReader::new(File::open(a_filename)?);
    buf_reader.seek(SeekFrom::Start(offset))?;

    let mut sid_reuse = vec![f64::NAN; iid_count];
    for (a_sid_index, mut atb_row) in atb.axis_iter_mut(nd::Axis(0)).enumerate() {
        if log_frequency > 0 && a_sid_index % log_frequency == 0 {
            println!(
                "   working on train_sid_index={} of {} (iid_count={}, b_sid_count={})",
                a_sid_index, a_sid_count, iid_count, b_sid_count
            );
        }

        buf_reader.read_f64_into::<LittleEndian>(&mut sid_reuse)?;

        nd::par_azip!(
            (mut atb_element in atb_row.axis_iter_mut(nd::Axis(0)),
            b1_col in b1.axis_iter(nd::Axis(1)),
            mut aatb_col in aatb.axis_iter_mut(nd::Axis(1)))
        {
            let mut atbi = 0.0;
            for iid_index in 0..iid_count {
                atbi += sid_reuse[iid_index] * b1_col[iid_index];
            }
            atb_element[()] = atbi;
            for iid_index in 0..iid_count {
                aatb_col[iid_index] -= sid_reuse[iid_index] * atbi;
            }
        });
    }
    Ok(())
}

#[allow(dead_code)]
fn read_into_f64(src: &mut BufReader<File>, dst: &mut [f64]) -> std::io::Result<()> {
    src.read_f64_into::<LittleEndian>(dst)
}

#[allow(dead_code)]
fn read_into_f32(src: &mut BufReader<File>, dst: &mut [f32]) -> std::io::Result<()> {
    src.read_f32_into::<LittleEndian>(dst)
}

/* Here are Python algorithms that shows how to do a low-memory multiply A (or A.T) x B (or B.T)
   They are used by file_ata_piece and file_aat_piece with some optimizations for A and B being the same.

output_list = [np.zeros((4,4)) for i in range(4)]

# a.T.dot(b)
for a_col2 in range(0,4,2): # 1 pass through A, returning output chunk about the same size writing in one pass
    buffer_a2 = a[:,a_col2:a_col2+2]
    for b_col in range(4): # A1/a1 passes through B
        buffer_b = b[:,b_col]
        for i in range(4):
            b_val = buffer_b[i]
            a_slice = buffer_a2[i,:]
            for k in range(2): # A1/a1 * A0 passes through the output
                output_list[0][a_col2+k,b_col] += a_slice[k]*b_val

# a.dot(b.T)
for out_col2 in range(0,4,2): # 1 pass through output, returning chunk on each pass
    for col in range(4): # O1/o1 passes through A and B
        buffer_a = a[:,col]
        buffer_b = b[:,col]
        for k in range(2):
            for i in range(4):
                output_list[1][i,out_col2+k] += buffer_a[i]*buffer_b[out_col2+k]

# a.T.dot(b.T)
for a_col2 in range(0,4,2): # 1 pass through A, returning an output chunk on each pass
    buffer_a2 = a[:,a_col2:a_col2+2]
    for b_col in range(4):
        buffer_b = b[:,b_col]
        for i in range(4):
            b_val = buffer_b[i]
            for k in range(2):
                output_list[2][a_col2+k,i] += buffer_a2[b_col,k]*b_val

# a.dot(b)  - but should instead do  (b.T.dot(a.T)).T
for b_col2 in range(0,4,2): #Transpose of preceding one
    buffer_b2 = b[:,b_col2:b_col2+2]
    for a_col in range(4):
        buffer_a = a[:,a_col]
        for i in range(4):
            a_val = buffer_a[i]
            for k in range(2):
                output_list[3][i,b_col2+k] += buffer_b2[a_col,k]*a_val


for output in output_list:
    print(output)
 */

// Given A, a matrix in Fortran order in a file
// with row_count rows and col_count columns,
// and given a starting column,
// returns part of A.T x A, the column vs column product.
// The piece piece returned has dimensions
// (col_count-col_start) x ncols
// where ncols <= (col_count-col_start)
// Makes only one pass through the file.
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
#[anyinput]
fn file_ata_piece<T: Float + Send + Sync + Sync + AddAssign>(
    path: AnyPath,
    offset: u64,
    row_count: usize,
    col_count: usize,
    col_start: usize,
    ata_piece: &mut nd::ArrayViewMut2<'_, T>,
    log_frequency: usize,
    read_into: fn(&mut BufReader<File>, &mut [T]) -> std::io::Result<()>,
) -> Result<(), Box<BedErrorPlus>> {
    let (nrows, ncols) = ata_piece.dim();
    if (col_start >= col_count)
        || (col_start + nrows != col_count)
        || (col_start + ncols > col_count)
    {
        return Err(Box::new(BedErrorPlus::BedError(
            BedError::CannotConvertBetaToFromF64,
        )));
    }

    _file_ata_piece_internal(
        path,
        offset,
        row_count,
        col_start,
        ata_piece,
        log_frequency,
        read_into,
    )
}

#[allow(dead_code)]
#[anyinput]
fn _file_ata_piece_internal<T: Float + Send + Sync + Sync + AddAssign>(
    path: AnyPath,
    offset: u64,
    row_count: usize,
    col_start: usize,
    ata_piece: &mut nd::ArrayViewMut2<'_, T>,
    log_frequency: usize,
    read_into: fn(&mut BufReader<File>, &mut [T]) -> std::io::Result<()>,
) -> Result<(), Box<BedErrorPlus>> {
    let (nrows, ncols) = ata_piece.dim();
    if log_frequency > 0 {
        println!(
            "file_ata_piece: col_start={}, {}x{} output",
            col_start, nrows, ncols
        );
    };

    // Open the file and move to the starting col
    let mut buf_reader = BufReader::new(File::open(path)?);
    buf_reader.seek(SeekFrom::Start(
        offset + col_start as u64 * row_count as u64 * std::mem::size_of::<T>() as u64,
    ))?;

    let mut col_save_list: Vec<Vec<T>> = vec![];
    let mut col_reuse = vec![T::nan(); row_count];

    for (col_rel_index, mut ata_row) in ata_piece.axis_iter_mut(nd::Axis(0)).enumerate() {
        if log_frequency > 0 && col_rel_index % log_frequency == 0 {
            println!("   working on {} of {}", col_rel_index, nrows);
        }

        // Read next col and save if in range
        let col = if col_save_list.len() < ncols {
            let mut col_save = vec![T::nan(); row_count];
            read_into(&mut buf_reader, &mut col_save)?;
            col_save_list.push(col_save);
            col_save_list.last().unwrap() // unwrap is OK here
        } else {
            read_into(&mut buf_reader, &mut col_reuse)?;
            &col_reuse
        };

        // Multiple saved sids with new sid
        let mut ata_row_trimmed = ata_row.slice_mut(nd::s![..col_save_list.len()]);
        nd::par_azip!((
            col_in_range in &col_save_list,
            mut ata_val in ata_row_trimmed.axis_iter_mut(nd::Axis(0))
        )
        {
            ata_val[()] = col_product(col_in_range, col);
        });
    }

    // Reflect the new product values
    for row_index in 0usize..ncols - 1 {
        for col_index in row_index..ncols {
            ata_piece[(row_index, col_index)] = ata_piece[(col_index, row_index)];
        }
    }
    Ok(())
}

#[allow(dead_code)]
fn col_product<T: Float + AddAssign>(col_i: &[T], col_j: &[T]) -> T {
    assert!(col_i.len() == col_j.len()); // real assert
    let mut product = T::zero();
    for row_index in 0..col_i.len() {
        product += col_i[row_index] * col_j[row_index];
    }
    product
}

// Given A, a matrix in Fortran order in a file
// with row_count rows and col_count columns,
// and given a starting column,
// returns part of A x A.T, the row vs row product.
// The piece piece returned has dimensions
// (row_count-row_start) x ncols
// where ncols <= (row_count-row_start)
// Makes only one pass through the file.
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
#[anyinput]
fn file_aat_piece<T: Float + Sync + Send + Sync + AddAssign>(
    path: AnyPath,
    offset: u64,
    row_count: usize,
    col_count: usize,
    row_start: usize,
    aat_piece: &mut nd::ArrayViewMut2<'_, T>,
    log_frequency: usize,
    read_into: fn(&mut BufReader<File>, &mut [T]) -> std::io::Result<()>,
) -> Result<(), Box<BedErrorPlus>> {
    let (nrows, ncols) = aat_piece.dim();

    if log_frequency > 0 {
        println!(
            "file_aat_piece: row_start={}, {}x{} output",
            row_start, nrows, ncols
        );
    };

    if (row_start >= row_count)
        || (row_start + nrows != row_count)
        || (row_start + ncols > row_count)
    {
        return Err(Box::new(BedErrorPlus::BedError(
            BedError::CannotConvertBetaToFromF64,
        )));
    }

    aat_piece.fill(T::zero());

    // Open the file and move to the starting col
    let mut buf_reader = BufReader::new(File::open(path)?);

    let mut col = vec![T::nan(); row_count - row_start];

    for col_index in 0..col_count {
        if log_frequency > 0 && col_index % log_frequency == 0 {
            println!("   working on {} of {}", col_index, col_count);
        }

        // Read next col
        buf_reader.seek(SeekFrom::Start(
            offset + (col_index * row_count + row_start) as u64 * std::mem::size_of::<T>() as u64,
        ))?;
        read_into(&mut buf_reader, &mut col)?;

        nd::par_azip!(
            (index row_index1,
            mut aat_col in aat_piece.axis_iter_mut(nd::Axis(1))
        )
        {
            let val1 = col[row_index1];
            for row_index0 in row_index1..nrows {
                aat_col[row_index0] += val1 * col[row_index0];
            }
        });
    }

    // Notice that ata reflects and aat doesn't. They don't need
    // to be the same, but they could be.
    Ok(())
}

// References: https://www.youtube.com/watch?v=0zOg8_B71gE&t=22s
// https://deterministic.space/elegant-apis-in-rust.html
// https://rust-lang.github.io/api-guidelines/
// https://ricardomartins.cc/2016/08/03/convenient_and_idiomatic_conversions_in_rust

/// Represents the metadata from PLINK .fam and .bim files.
///
/// Construct with [`Metadata::builder`](struct.Metadata.html#method.builder) or [`Metadata::new`](struct.Metadata.html#method.new).
///
/// # Example
///
/// Extract metadata from a file.
/// Create a random file with the same metadata.
/// ```
/// use ndarray as nd;
/// use bed_reader::{Bed, WriteOptions, sample_bed_file};
/// use ndarray_rand::{rand::prelude::StdRng, rand::SeedableRng, rand_distr::Uniform, RandomExt};
///
/// let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
/// let metadata = bed.metadata()?;
/// let shape = bed.dim()?;
///
/// let mut rng = StdRng::seed_from_u64(0);
/// let val = nd::Array::random_using(shape, Uniform::from(-1..3), &mut rng);
///
/// let temp_out = temp_testdir::TempDir::default();
/// let output_file = temp_out.join("random.bed");
/// WriteOptions::builder(output_file)
///     .metadata(&metadata)
///     .missing_value(-1)
///     .write(&val)?;
/// # use bed_reader::BedErrorPlus;
/// # Ok::<(), Box<BedErrorPlus>>(())
/// ```
#[derive(Clone, Debug, Builder, PartialEq)]
#[builder(build_fn(private, name = "build_no_file_check", error = "BedErrorPlus"))]
pub struct Metadata {
    #[builder(setter(custom))]
    #[builder(default = "None")]
    fid: Option<Rc<nd::Array1<String>>>,
    #[builder(setter(custom))]
    #[builder(default = "None")]
    iid: Option<Rc<nd::Array1<String>>>,
    #[builder(setter(custom))]
    #[builder(default = "None")]
    father: Option<Rc<nd::Array1<String>>>,
    #[builder(setter(custom))]
    #[builder(default = "None")]
    mother: Option<Rc<nd::Array1<String>>>,

    // i32 based on https://www.cog-genomics.org/plink2/formats#bim
    #[builder(setter(custom))]
    #[builder(default = "None")]
    sex: Option<Rc<nd::Array1<i32>>>,
    #[builder(setter(custom))]
    #[builder(default = "None")]
    pheno: Option<Rc<nd::Array1<String>>>,

    #[builder(setter(custom))]
    #[builder(default = "None")]
    chromosome: Option<Rc<nd::Array1<String>>>,
    #[builder(setter(custom))]
    #[builder(default = "None")]
    sid: Option<Rc<nd::Array1<String>>>,
    #[builder(setter(custom))]
    #[builder(default = "None")]
    cm_position: Option<Rc<nd::Array1<f32>>>,
    #[builder(setter(custom))]
    #[builder(default = "None")]
    bp_position: Option<Rc<nd::Array1<i32>>>,
    #[builder(setter(custom))]
    #[builder(default = "None")]
    allele_1: Option<Rc<nd::Array1<String>>>,
    #[builder(setter(custom))]
    #[builder(default = "None")]
    allele_2: Option<Rc<nd::Array1<String>>>,
}

fn lazy_or_skip_count<T>(array: &Option<Rc<nd::Array1<T>>>) -> Option<usize> {
    array.as_ref().map(|array| array.len())
}

/// Represents a PLINK .bed file that is open for reading genotype data and metadata.
///
/// Construct with [`Bed::new`](struct.Bed.html#method.new) or [`Bed::builder`](struct.Bed.html#method.builder).
///
/// # Example
///
/// Open a file for reading. Then, read the individual (sample) ids
/// and all the genotype data.
/// ```
/// use ndarray as nd;
/// use bed_reader::{Bed, ReadOptions, sample_bed_file};
/// use bed_reader::assert_eq_nan;
///
/// let file_name = sample_bed_file("small.bed")?;
/// let mut bed = Bed::new(file_name)?;
/// println!("{:?}", bed.iid()?); // Outputs ndarray ["iid1", "iid2", "iid3"]
/// let val = ReadOptions::builder().f64().read(&mut bed)?;
///
/// assert_eq_nan(
///     &val,
///     &nd::array![
///         [1.0, 0.0, f64::NAN, 0.0],
///         [2.0, 0.0, f64::NAN, 2.0],
///         [0.0, 1.0, 2.0, 0.0]
///     ],
/// );
/// # use bed_reader::BedErrorPlus;
/// # Ok::<(), Box<BedErrorPlus>>(())
/// ```
#[derive(Clone, Debug, Builder)]
#[builder(build_fn(private, name = "build_no_file_check", error = "BedErrorPlus"))]
pub struct Bed {
    // https://stackoverflow.com/questions/32730714/what-is-the-right-way-to-store-an-immutable-path-in-a-struct
    // don't emit a setter, but keep the field declaration on the builder
    /// The file name or path of the .bed file.
    #[builder(setter(custom))]
    path: PathBuf,

    #[builder(setter(custom))]
    #[builder(default = "None")]
    fam_path: Option<PathBuf>,

    #[builder(setter(custom))]
    #[builder(default = "None")]
    bim_path: Option<PathBuf>,

    #[builder(setter(custom))]
    #[builder(default = "true")]
    is_checked_early: bool,

    #[builder(setter(custom))]
    #[builder(default = "None")]
    iid_count: Option<usize>,

    #[builder(setter(custom))]
    #[builder(default = "None")]
    sid_count: Option<usize>,

    #[builder(setter(custom))]
    metadata: Metadata,

    #[builder(setter(custom))]
    skip_set: HashSet<MetadataFields>,
}

/// All Metadata fields.
///
/// Used by [`Metadata::read_fam`](struct.Metadata.html#method.read_fam) and
/// [`Metadata::read_bim`](struct.Metadata.html#method.read_bim) to skip reading
/// specified metadata fields.
#[derive(Debug, PartialEq, Eq, Copy, Clone, Ord, PartialOrd, Hash)]
pub enum MetadataFields {
    #[allow(missing_docs)]
    Fid,
    #[allow(missing_docs)]
    Iid,
    #[allow(missing_docs)]
    Father,
    #[allow(missing_docs)]
    Mother,
    #[allow(missing_docs)]
    Sex,
    #[allow(missing_docs)]
    Pheno,
    #[allow(missing_docs)]
    Chromosome,
    #[allow(missing_docs)]
    Sid,
    #[allow(missing_docs)]
    CmPosition,
    #[allow(missing_docs)]
    BpPosition,
    #[allow(missing_docs)]
    Allele1,
    #[allow(missing_docs)]
    Allele2,
}

impl BedBuilder {
    #[anyinput]
    fn new(path: AnyPath) -> Self {
        Self {
            path: Some(path.to_owned()),
            fam_path: None,
            bim_path: None,

            is_checked_early: None,
            iid_count: None,
            sid_count: None,

            metadata: Some(Metadata::new()),
            skip_set: Some(HashSet::new()),
        }
    }

    /// Create [`Bed`](struct.Bed.html) from the builder.
    ///
    /// > See [`Bed::builder`](struct.Bed.html#method.builder) for more details and examples.
    pub fn build(&self) -> Result<Bed, Box<BedErrorPlus>> {
        let mut bed = self.build_no_file_check()?;

        if bed.is_checked_early {
            open_and_check(&bed.path)?;
        }

        (bed.iid_count, bed.sid_count) = bed.metadata.check_counts(bed.iid_count, bed.sid_count)?;

        Ok(bed)
    }

    // https://stackoverflow.com/questions/38183551/concisely-initializing-a-vector-of-strings
    // https://stackoverflow.com/questions/65250496/how-to-convert-intoiteratoritem-asrefstr-to-iteratoritem-str-in-rust

    /// Override the family id (fid) values found in the .fam file.
    ///
    /// By default, if fid values are needed and haven't already been found,
    /// they will be read from the .fam file.
    /// Providing them here avoids that file read and provides a way to give different values.
    #[anyinput]
    pub fn fid(mut self, fid: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because BedBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_fid(fid);
        self
    }

    /// Override the individual id (iid) values found in the .fam file.
    ///
    /// By default, if iid values are needed and haven't already been found,
    /// they will be read from the .fam file.
    /// Providing them here avoids that file read and provides a way to give different values.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, assert_eq_nan, sample_bed_file};
    /// let file_name = sample_bed_file("small.bed")?;
    /// use bed_reader::ReadOptions;
    ///
    /// let mut bed = Bed::builder(file_name)
    ///    .iid(["sample1", "sample2", "sample3"])
    ///    .build()?;
    /// println!("{:?}", bed.iid()?); // Outputs ndarray ["sample1", "sample2", "sample3"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn iid(mut self, iid: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because BedBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_iid(iid);
        self
    }

    /// Override the father values found in the .fam file.
    ///nd
    /// By default, if father values are needed and haven't already been found,
    /// they will be read from the .fam file.
    /// Providing them here avoids that file read and provides a way to gi&ve different values.
    #[anyinput]
    pub fn father(mut self, father: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because BedBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_father(father);
        self
    }

    /// Override the mother values found in the .fam file.
    ///
    /// By default, if mother values are needed and haven't already been found,
    /// they will be read from the .fam file.
    /// Providing them here avoids that file read and provides a way to give different values.
    #[anyinput]
    pub fn mother(mut self, mother: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because BedBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_mother(mother);
        self
    }

    /// Override the sex values found in the .fam file.
    ///
    /// By default, if sex values are needed and haven't already been found,
    /// they will be read from the .fam file.
    /// Providing them here avoids that file read and provides a way to give different values.
    #[anyinput]
    pub fn sex(mut self, sex: AnyIter<i32>) -> Self {
        // Unwrap will always work because BedBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_sex(sex);
        self
    }

    /// Override the phenotype values found in the .fam file.
    ///
    /// Note that the phenotype values in the .fam file are seldom used.
    /// By default, if phenotype values are needed and haven't already been found,
    /// they will be read from the .fam file.
    /// Providing them here avoids that file read and provides a way to give different values.
    #[anyinput]
    pub fn pheno(mut self, pheno: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because BedBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_pheno(pheno);
        self
    }

    /// Override the chromosome values found in the .bim file.
    ///
    /// By default, if chromosome values are needed and haven't already been found,
    /// they will be read from the .bim file.
    /// Providing them here avoids that file read and provides a way to give different values.
    #[anyinput]
    pub fn chromosome(mut self, chromosome: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because BedBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_chromosome(chromosome);
        self
    }

    /// Override the SNP id (sid) values found in the .fam file.
    ///
    /// By default, if sid values are needed and haven't already been found,
    /// they will be read from the .bim file.
    /// Providing them here avoids that file read and provides a way to give different values.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, assert_eq_nan, sample_bed_file};
    /// let file_name = sample_bed_file("small.bed")?;
    ///
    /// let mut bed = Bed::builder(file_name)
    ///    .sid(["SNP1", "SNP2", "SNP3", "SNP4"])
    ///    .build()?;
    /// println!("{:?}", bed.sid()?); // Outputs ndarray ["SNP1", "SNP2", "SNP3", "SNP4"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn sid(mut self, sid: AnyIter<AnyString>) -> Self {
        self.metadata.as_mut().unwrap().set_sid(sid);
        self
    }

    /// Override the centimorgan position values found in the .bim file.
    ///
    /// By default, if centimorgan position values are needed and haven't already been found,
    /// they will be read from the .bim file.
    /// Providing them here avoids that file read and provides a way to give different values.
    #[anyinput]
    pub fn cm_position(mut self, cm_position: AnyIter<f32>) -> Self {
        // Unwrap will always work because BedBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_cm_position(cm_position);
        self
    }

    /// Override the base-pair position values found in the .bim file.
    ///
    /// By default, if base-pair position values are needed and haven't already been found,
    /// they will be read from the .bim file.
    /// Providing them here avoids that file read and provides a way to give different values.
    #[anyinput]
    pub fn bp_position(mut self, bp_position: AnyIter<i32>) -> Self {
        // Unwrap will always work because BedBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_bp_position(bp_position);
        self
    }

    /// Override the allele 1 values found in the .bim file.
    ///
    /// By default, if allele 1 values are needed and haven't already been found,
    /// they will be read from the .bim file.
    /// Providing them here avoids that file read and provides a way to give different values.
    #[anyinput]
    pub fn allele_1(mut self, allele_1: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because BedBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_allele_1(allele_1);
        self
    }

    /// Override the allele 2 values found in the .bim file.
    ///
    /// By default, if allele 2 values are needed and haven't already been found,
    /// they will be read from the .bim file.
    /// Providing them here avoids that file read and provides a way to give different values.
    #[anyinput]
    pub fn allele_2(mut self, allele_2: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because BedBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_allele_2(allele_2);
        self
    }

    /// Set the number of individuals (samples) in the data.
    ///
    /// By default, if this number is needed, it will be found
    /// and remembered
    /// by opening the .fam file and quickly counting the number
    /// of lines. Providing the number thus avoids a file read.
    pub fn iid_count(mut self, count: usize) -> Self {
        self.iid_count = Some(Some(count));
        self
    }

    /// Set the number of SNPs in the data.
    ///
    /// By default, if this number is needed, it will be found
    /// and remembered
    /// by opening the .bim file and quickly counting the number
    /// of lines. Providing the number thus avoids a file read.
    pub fn sid_count(mut self, count: usize) -> Self {
        self.sid_count = Some(Some(count));
        self
    }

    /// Don't check the header of the .bed file until and unless the file is actually read.
    ///
    /// By default, when a [`Bed`](struct.Bed.html) struct is created, the .bed
    /// file header is checked. This stops that early check.
    pub fn skip_early_check(mut self) -> Self {
        self.is_checked_early = Some(false);
        self
    }

    /// Set the path to the .fam file.
    ///
    /// If not set, the .fam file will be assumed
    /// to have the same name as the .bed file, but with the extension .fam.
    ///
    /// # Example:
    /// Read .bed, .fam, and .bim files with non-standard names.
    /// ```
    /// use bed_reader::{Bed, ReadOptions, sample_files};
    /// let deb_maf_mib = sample_files(["small.deb", "small.maf", "small.mib"])?;
    /// let mut bed = Bed::builder(&deb_maf_mib[0])
    ///    .fam_path(&deb_maf_mib[1])
    ///    .bim_path(&deb_maf_mib[2])
    ///    .build()?;
    /// println!("{:?}", bed.iid()?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// println!("{:?}", bed.sid()?); // Outputs ndarray ["sid1", "sid2", "sid3", "sid4"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn fam_path(mut self, path: AnyPath) -> Self {
        self.fam_path = Some(Some(path.to_owned()));
        self
    }

    /// Set the path to the .bim file.
    ///
    /// If not set, the .bim file will be assumed
    /// to have the same name as the .bed file, but with the extension .bim.
    ///
    /// # Example:
    /// Read .bed, .fam, and .bim files with non-standard names.
    /// ```
    /// use bed_reader::{Bed, ReadOptions, sample_files};
    /// let deb_maf_mib = sample_files(["small.deb", "small.maf", "small.mib"])?;
    /// let mut bed = Bed::builder(&deb_maf_mib[0])
    ///    .fam_path(&deb_maf_mib[1])
    ///    .bim_path(&deb_maf_mib[2])
    ///    .build()?;
    /// println!("{:?}", bed.iid()?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// println!("{:?}", bed.sid()?); // Outputs ndarray ["sid1", "sid2", "sid3", "sid4"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn bim_path(mut self, path: AnyPath) -> Self {
        self.bim_path = Some(Some(path.to_owned()));
        self
    }

    /// Don't read the fid information from the .fam file.
    ///
    /// By default, when the .fam is read, the fid (the family id) is recorded.
    /// This stops that recording. This is useful if the fid is not needed.
    /// Asking for the fid after skipping it results in an error.    
    pub fn skip_fid(mut self) -> Self {
        // Unwrap will always work because BedBuilder starting with some skip_set
        self.skip_set.as_mut().unwrap().insert(MetadataFields::Fid);
        self
    }

    /// Don't read the iid information from the .fam file.
    ///
    /// By default, when the .fam is read, the iid (the individual id) is recorded.
    /// This stops that recording. This is useful if the iid is not needed.
    /// Asking for the iid after skipping it results in an error.
    pub fn skip_iid(mut self) -> Self {
        // Unwrap will always work because BedBuilder starting with some skip_set
        self.skip_set.as_mut().unwrap().insert(MetadataFields::Iid);
        self
    }

    /// Don't read the father information from the .fam file.
    ///
    /// By default, when the .fam is read, the father id is recorded.
    /// This stops that recording. This is useful if the father id is not needed.
    /// Asking for the father id after skipping it results in an error.    
    pub fn skip_father(mut self) -> Self {
        // Unwrap will always work because BedBuilder starting with some skip_set
        self.skip_set
            .as_mut()
            .unwrap()
            .insert(MetadataFields::Father);
        self
    }

    /// Don't read the mother information from the .fam file.
    ///
    /// By default, when the .fam is read, the mother id is recorded.
    /// This stops that recording. This is useful if the mother id is not needed.
    /// Asking for the mother id after skipping it results in an error.    
    pub fn skip_mother(mut self) -> Self {
        // Unwrap will always work because BedBuilder starting with some skip_set
        self.skip_set
            .as_mut()
            .unwrap()
            .insert(MetadataFields::Mother);
        self
    }

    /// Don't read the sex information from the .fam file.
    ///
    /// By default, when the .fam is read, the sex is recorded.
    /// This stops that recording. This is useful if sex is not needed.
    /// Asking for sex after skipping it results in an error.    
    pub fn skip_sex(mut self) -> Self {
        // Unwrap will always work because BedBuilder starting with some skip_set
        self.skip_set.as_mut().unwrap().insert(MetadataFields::Sex);
        self
    }

    /// Don't read the phenotype information from the .fam file.
    ///
    /// Note that the phenotype information in the .fam file is
    /// seldom used.
    ///
    /// By default, when the .fam is read, the phenotype is recorded.
    /// This stops that recording. This is useful if this phenotype
    /// information is not needed.
    /// Asking for the phenotype after skipping it results in an error.    
    pub fn skip_pheno(mut self) -> Self {
        // Unwrap will always work because BedBuilder starting with some skip_set
        self.skip_set
            .as_mut()
            .unwrap()
            .insert(MetadataFields::Pheno);
        self
    }

    /// Don't read the chromosome information from the .bim file.
    ///
    /// By default, when the .bim is read, the chromosome is recorded.
    /// This stops that recording. This is useful if the chromosome is not needed.
    /// Asking for the chromosome after skipping it results in an error.    
    pub fn skip_chromosome(mut self) -> Self {
        // Unwrap will always work because BedBuilder starting with some skip_set
        self.skip_set
            .as_mut()
            .unwrap()
            .insert(MetadataFields::Chromosome);
        self
    }

    /// Don't read the SNP id information from the .bim file.
    ///
    /// By default, when the .bim is read, the sid (SNP id) is recorded.
    /// This stops that recording. This is useful if the sid is not needed.
    /// Asking for the sid after skipping it results in an error.    
    pub fn skip_sid(mut self) -> Self {
        // Unwrap will always work because BedBuilder starting with some skip_set
        self.skip_set.as_mut().unwrap().insert(MetadataFields::Sid);
        self
    }

    /// Don't read the centimorgan position information from the .bim file.
    ///
    /// By default, when the .bim is read, the cm position is recorded.
    /// This stops that recording. This is useful if the cm position is not needed.
    /// Asking for the cm position after skipping it results in an error.    
    pub fn skip_cm_position(mut self) -> Self {
        // Unwrap will always work because BedBuilder starting with some skip_set
        self.skip_set
            .as_mut()
            .unwrap()
            .insert(MetadataFields::CmPosition);
        self
    }

    /// Don't read the base-pair position information from the .bim file.
    ///
    /// By default, when the .bim is read, the bp position is recorded.
    /// This stops that recording. This is useful if the bp position is not needed.
    /// Asking for the cp position after skipping it results in an error.    
    pub fn skip_bp_position(mut self) -> Self {
        // Unwrap will always work because BedBuilder starting with some skip_set
        self.skip_set
            .as_mut()
            .unwrap()
            .insert(MetadataFields::BpPosition);
        self
    }

    /// Don't read the allele 1 information from the .bim file.
    ///
    /// By default, when the .bim is read, allele 1 is recorded.
    /// This stops that recording. This is useful if allele 1 is not needed.
    /// Asking for allele 1 after skipping it results in an error.    
    pub fn skip_allele_1(mut self) -> Self {
        // Unwrap will always work because BedBuilder starting with some skip_set
        self.skip_set
            .as_mut()
            .unwrap()
            .insert(MetadataFields::Allele1);
        self
    }

    /// Don't read the allele 2 information from the .bim file.
    ///
    /// By default, when the .bim is read, allele 2 is recorded.
    /// This stops that recording. This is useful if allele 2 is not needed.
    /// Asking for allele 2 after skipping it results in an error.    
    pub fn skip_allele_2(mut self) -> Self {
        // Unwrap will always work because BedBuilder starting with some skip_set
        self.skip_set
            .as_mut()
            .unwrap()
            .insert(MetadataFields::Allele2);
        self
    }

    /// Override the metadata in the .fam and .bim files with info merged in from a [`Metadata`](struct.Metadata.html).
    ///
    /// # Example
    ///
    /// In the example, we create a [`Metadata`](struct.Metadata.html) with iid
    /// and sid arrays. Next, we use [`BedBuilder`](struct.BedBuilder.html) to override the fid array
    /// and an iid array. Then, we add the metadata to the [`BedBuilder`](struct.BedBuilder.html),
    /// overwriting iid (again) and overriding sid. Finally, we print these
    /// three arrays and chromosome. Chromosome was never overridden so
    /// it is read from the *.bim file.
    ///```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, Metadata, sample_bed_file};
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let metadata = Metadata::builder()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build()?;
    /// let mut bed = Bed::builder(file_name)
    ///     .fid(["f1", "f2", "f3"])
    ///     .iid(["x1", "x2", "x3"])
    ///     .metadata(&metadata)
    ///     .build()?;
    /// println!("{0:?}", bed.fid()?);  // Outputs ndarray ["f1", "f2", "f3"]
    /// println!("{0:?}", bed.iid()?);  // Outputs ndarray ["i1", "i2", "i3"]
    /// println!("{0:?}", bed.sid()?);  // Outputs ndarray ["s1", "s2", "s3", "s4"]
    /// println!("{0:?}", bed.chromosome()?);  // Outputs ndarray ["1", "1", "5", "Y"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn metadata(mut self, metadata: &Metadata) -> Self {
        self.metadata = Some(
            Metadata::builder()
                .metadata(&self.metadata.unwrap()) // unwrap is ok because we know we have metadata
                .metadata(metadata) // consistent counts will be check later by the BedBuilder
                .build_no_file_check()
                .unwrap(), // unwrap is ok because nothing can go wrong
        );

        self
    }
}

#[anyinput]
fn to_metadata_path(
    bed_path: AnyPath,
    metadata_path: &Option<PathBuf>,
    extension: AnyString,
) -> PathBuf {
    if let Some(metadata_path) = metadata_path {
        metadata_path.to_owned()
    } else {
        bed_path.with_extension(extension)
    }
}

impl Bed {
    /// Attempts to open a PLINK .bed file for reading. Supports options.
    ///
    /// > Also see [`Bed::new`](struct.Bed.html#method.new), which does not support options.
    ///
    /// The options, [listed here](struct.BedBuilder.html#implementations), can:
    ///  * set the path of the .fam and/or .bim file
    ///  * override some metadata, for example, replace the individual ids.
    ///  * set the number of individuals (samples) or SNPs (variants)
    ///  * control checking the validity of the .bed file's header
    ///  * skip reading selected metadata
    ///
    /// Note that this method is a lazy about holding files, so unlike `std::fs::File::open(&path)`, it
    /// will not necessarily lock the file(s).
    ///
    /// # Errors
    /// By default, this method will return an error if the file is missing or its header
    /// is ill-formed. It will also return an error if the options contradict each other.
    /// See [`BedError`](enum.BedError.html) and [`BedErrorPlus`](enum.BedErrorPlus.html)
    /// for all possible errors.
    ///
    /// # Examples
    /// List individual (sample) [`iid`](struct.Bed.html#method.iid) and
    /// SNP (variant) [`sid`](struct.Bed.html#method.sid),
    /// then [`read`](struct.Bed.html#method.read) the whole file.
    ///
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, assert_eq_nan, sample_bed_file};
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::builder(file_name).build()?;
    /// println!("{:?}", bed.iid()?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// println!("{:?}", bed.sid()?); // Outputs ndarray ["snp1", "snp2", "snp3", "snp4"]
    /// let val = bed.read::<f64>()?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1.0, 0.0, f64::NAN, 0.0],
    ///         [2.0, 0.0, f64::NAN, 2.0],
    ///         [0.0, 1.0, 2.0, 0.0]
    ///     ],
    /// );
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    ///
    /// Replace [`iid`](struct.Bed.html#method.iid).
    /// ```
    /// # use ndarray as nd;
    /// # use bed_reader::{Bed, ReadOptions, assert_eq_nan, sample_bed_file};
    /// # let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::builder(file_name)
    ///    .iid(["sample1", "sample2", "sample3"])
    ///    .build()?;
    /// println!("{:?}", bed.iid()?); // Outputs ndarray ["sample1", "sample2", "sample3"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    /// Give the number of individuals (samples) and SNPs (variants) so that the .fam and
    /// .bim files need never be opened.
    /// ```
    /// # use ndarray as nd;
    /// # use bed_reader::{Bed, ReadOptions, assert_eq_nan, sample_bed_file};
    /// # let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::builder(file_name).iid_count(3).sid_count(4).build()?;
    /// let val = bed.read::<f64>()?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1.0, 0.0, f64::NAN, 0.0],
    ///         [2.0, 0.0, f64::NAN, 2.0],
    ///         [0.0, 1.0, 2.0, 0.0]
    ///     ],
    /// );
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    /// Mark some properties as "don’t read or offer".
    /// ```
    /// # use ndarray as nd;
    /// # use bed_reader::{Bed, ReadOptions, assert_eq_nan, sample_bed_file};
    /// # let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::builder(file_name)
    ///     .skip_father()
    ///     .skip_mother()
    ///     .skip_sex()
    ///     .skip_pheno()
    ///     .skip_allele_1()
    ///     .skip_allele_2()
    ///     .build()?;
    /// println!("{:?}", bed.iid()?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// bed.allele_2().expect_err("Can't be read");
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    ///
    #[anyinput]
    pub fn builder(path: AnyPath) -> BedBuilder {
        BedBuilder::new(path)
    }

    /// Attempts to open a PLINK .bed file for reading. Does not support options.
    ///
    /// > Also see [`Bed::builder`](struct.Bed.html#method.builder), which does support options.
    ///
    /// Note that this method is a lazy about holding files, so unlike `std::fs::File::open(&path)`, it
    /// will not necessarily lock the file(s).
    ///
    /// # Errors
    /// By default, this method will return an error if the file is missing or its header
    /// is ill-formed. See [`BedError`](enum.BedError.html) and [`BedErrorPlus`](enum.BedErrorPlus.html)
    /// for all possible errors.
    ///
    /// # Examples
    /// List individual (sample) [`iid`](struct.Bed.html#method.iid) and
    /// SNP (variant) [`sid`](struct.Bed.html#method.sid),
    /// then [`read`](struct.Bed.html#method.read) the whole file.
    ///
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, assert_eq_nan, sample_bed_file};
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// println!("{:?}", bed.iid()?); // Outputs ndarray: ["iid1", "iid2", "iid3"]
    /// println!("{:?}", bed.sid()?); // Outputs ndarray: ["sid1", "sid2", "sid3", "sid4"]
    /// let val = bed.read::<f64>()?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1.0, 0.0, f64::NAN, 0.0],
    ///         [2.0, 0.0, f64::NAN, 2.0],
    ///         [0.0, 1.0, 2.0, 0.0]
    ///     ],
    /// );
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    ///
    /// Open the file and read data for one SNP (variant)
    /// at index position 2.
    /// ```
    /// # use ndarray as nd;
    /// # use bed_reader::{Bed, ReadOptions, assert_eq_nan, sample_bed_file};
    /// # let file_name = sample_bed_file("small.bed")?;
    ///
    /// let mut bed = Bed::new(file_name)?;
    /// let val = ReadOptions::builder().sid_index(2).f64().read(&mut bed)?;
    ///
    /// assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn new(path: AnyPath) -> Result<Self, Box<BedErrorPlus>> {
        Bed::builder(path).build()
    }

    /// Number of individuals (samples)
    ///
    /// If this number is needed, it will be found
    /// by opening the .fam file and quickly counting the number
    /// of lines. Once found, the number will be remembered.
    /// The file read can be avoided by setting the
    /// number with [`BedBuilder::iid_count`](struct.BedBuilder.html#method.iid_count)
    /// or, for example, [`BedBuilder::iid`](struct.BedBuilder.html#method.iid).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, assert_eq_nan, sample_bed_file};
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let iid_count = bed.iid_count()?;
    ///
    /// assert!(iid_count == 3);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn iid_count(&mut self) -> Result<usize, Box<BedErrorPlus>> {
        if let Some(iid_count) = self.iid_count {
            Ok(iid_count)
        } else {
            let fam_path = self.fam_path();
            let iid_count = count_lines(fam_path)?;
            self.iid_count = Some(iid_count);
            Ok(iid_count)
        }
    }

    /// Number of SNPs (variants)
    ///
    /// If this number is needed, it will be found
    /// by opening the .bim file and quickly counting the number
    /// of lines. Once found, the number will be remembered.
    /// The file read can be avoided by setting the
    /// number with [`BedBuilder::sid_count`](struct.BedBuilder.html#method.sid_count)
    /// or, for example, [`BedBuilder::sid`](struct.BedBuilder.html#method.sid).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, assert_eq_nan, sample_bed_file};
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let sid_count = bed.sid_count()?;
    ///
    /// assert!(sid_count == 4);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn sid_count(&mut self) -> Result<usize, Box<BedErrorPlus>> {
        if let Some(sid_count) = self.sid_count {
            Ok(sid_count)
        } else {
            let bim_path = self.bim_path();
            let sid_count = count_lines(bim_path)?;
            self.sid_count = Some(sid_count);
            Ok(sid_count)
        }
    }

    /// Number of individuals (samples) and SNPs (variants)
    ///
    /// If these numbers aren't known, they will be found
    /// by opening the .fam and .bim files and quickly counting the number
    /// of lines. Once found, the numbers will be remembered.
    /// The file read can be avoided by setting the
    /// number with [`BedBuilder::iid_count`](struct.BedBuilder.html#method.iid_count)
    /// and [`BedBuilder::sid_count`](struct.BedBuilder.html#method.sid_count).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let dim = bed.dim()?;
    ///
    /// assert!(dim == (3,4));
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn dim(&mut self) -> Result<(usize, usize), Box<BedErrorPlus>> {
        Ok((self.iid_count()?, self.sid_count()?))
    }

    /// Family id of each of individual (sample)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .fam file. Once found, this ndarray
    /// and other information in the .fam file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedBuilder::fid`](struct.BedBuilder.html#method.fid).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let fid = bed.fid()?;
    /// println!("{fid:?}"); // Outputs ndarray ["fid1", "fid1", "fid2"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn fid(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_fam::<String>(self.metadata.fid.is_none(), MetadataFields::Fid, "fid")?;
        Ok(self.metadata.fid.as_ref().unwrap()) //unwrap always works because of lazy_fam
    }

    /// Individual id of each of individual (sample)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .fam file. Once found, this ndarray
    /// and other information in the .fam file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedBuilder::iid`](struct.BedBuilder.html#method.iid).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let iid = bed.iid()?;    ///
    /// println!("{iid:?}"); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn iid(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_fam::<String>(self.metadata.iid.is_none(), MetadataFields::Iid, "iid")?;
        Ok(self.metadata.iid.as_ref().unwrap()) //unwrap always works because of lazy_fam
    }

    /// Father id of each of individual (sample)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .fam file. Once found, this ndarray
    /// and other information in the .fam file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedBuilder::father`](struct.BedBuilder.html#method.father).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let father = bed.father()?;
    /// println!("{father:?}"); // Outputs ndarray ["iid23", "iid23", "iid22"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())    
    pub fn father(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_fam::<String>(
            self.metadata.father.is_none(),
            MetadataFields::Father,
            "father",
        )?;
        Ok(self.metadata.father.as_ref().unwrap()) //unwrap always works because of lazy_fam
    }

    /// Mother id of each of individual (sample)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .fam file. Once found, this ndarray
    /// and other information in the .fam file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedBuilder::mother`](struct.BedBuilder.html#method.mother).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let mother = bed.mother()?;
    /// println!("{mother:?}"); // Outputs ndarray ["iid34", "iid34", "iid33"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn mother(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_fam::<String>(
            self.metadata.mother.is_none(),
            MetadataFields::Mother,
            "mother",
        )?;
        Ok(self.metadata.mother.as_ref().unwrap()) //unwrap always works because of lazy_fam
    }

    /// Sex each of individual (sample)
    ///
    /// 0 is unknown, 1 is male, 2 is female
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .fam file. Once found, this ndarray
    /// and other information in the .fam file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedBuilder::sex`](struct.BedBuilder.html#method.sex).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let sex = bed.sex()?;
    /// println!("{sex:?}"); // Outputs ndarray [1, 2, 0]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn sex(&mut self) -> Result<&nd::Array1<i32>, Box<BedErrorPlus>> {
        self.unlazy_fam::<String>(self.metadata.sex.is_none(), MetadataFields::Sex, "sex")?;
        Ok(self.metadata.sex.as_ref().unwrap()) //unwrap always works because of lazy_fam
    }

    /// A phenotype for each individual (seldom used)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .fam file. Once found, this ndarray
    /// and other information in the .fam file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedBuilder::pheno`](struct.BedBuilder.html#method.pheno).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let pheno = bed.pheno()?;
    /// println!("{pheno:?}"); // Outputs ndarray ["red", "red", "blue"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn pheno(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_fam::<String>(
            self.metadata.pheno.is_none(),
            MetadataFields::Pheno,
            "pheno",
        )?;
        Ok(self.metadata.pheno.as_ref().unwrap()) //unwrap always works because of lazy_fam
    }

    /// Chromosome of each SNP (variant)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .bim file. Once found, this ndarray
    /// and other information in the .bim file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedBuilder::chromosome`](struct.BedBuilder.html#method.chromosome).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let chromosome = bed.chromosome()?;
    /// println!("{chromosome:?}"); // Outputs ndarray ["1", "1", "5", "Y"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn chromosome(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_bim::<String>(
            self.metadata.chromosome.is_none(),
            MetadataFields::Chromosome,
            "chromosome",
        )?;
        Ok(self.metadata.chromosome.as_ref().unwrap()) //unwrap always works because of lazy_bim
    }

    /// SNP id of each SNP (variant)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .bim file. Once found, this ndarray
    /// and other information in the .bim file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedBuilder::sid`](struct.BedBuilder.html#method.sid).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let sid = bed.sid()?;
    /// println!("{sid:?}"); // Outputs ndarray "sid1", "sid2", "sid3", "sid4"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn sid(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_bim::<String>(self.metadata.sid.is_none(), MetadataFields::Sid, "sid")?;
        Ok(self.metadata.sid.as_ref().unwrap()) //unwrap always works because of lazy_bim
    }

    /// Centimorgan position of each SNP (variant)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .bim file. Once found, this ndarray
    /// and other information in the .bim file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedBuilder::cm_position`](struct.BedBuilder.html#method.cm_position).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let cm_position = bed.cm_position()?;
    /// println!("{cm_position:?}"); // Outputs ndarray [100.4, 2000.5, 4000.7, 7000.9]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn cm_position(&mut self) -> Result<&nd::Array1<f32>, Box<BedErrorPlus>> {
        self.unlazy_bim::<String>(
            self.metadata.cm_position.is_none(),
            MetadataFields::CmPosition,
            "cm_position",
        )?;
        Ok(self.metadata.cm_position.as_ref().unwrap()) //unwrap always works because of lazy_bim
    }

    /// Base-pair position of each SNP (variant)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .bim file. Once found, this ndarray
    /// and other information in the .bim file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedBuilder::bp_position`](struct.BedBuilder.html#method.bp_position).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let bp_position = bed.bp_position()?;
    /// println!("{bp_position:?}"); // Outputs ndarray [1, 100, 1000, 1004]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn bp_position(&mut self) -> Result<&nd::Array1<i32>, Box<BedErrorPlus>> {
        self.unlazy_bim::<String>(
            self.metadata.bp_position.is_none(),
            MetadataFields::BpPosition,
            "bp_position",
        )?;
        Ok(self.metadata.bp_position.as_ref().unwrap()) //unwrap always works because of lazy_bim
    }

    /// First allele of each SNP (variant)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .bim file. Once found, this ndarray
    /// and other information in the .bim file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedBuilder::allele_1`](struct.BedBuilder.html#method.allele_1).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let allele_1 = bed.allele_1()?;
    /// println!("{allele_1:?}"); // Outputs ndarray ["A", "T", "A", "T"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn allele_1(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_bim::<String>(
            self.metadata.allele_1.is_none(),
            MetadataFields::Allele1,
            "allele_1",
        )?;
        Ok(self.metadata.allele_1.as_ref().unwrap()) //unwrap always works because of lazy_bim
    }

    /// Second allele of each SNP (variant)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .bim file. Once found, this ndarray
    /// and other information in the .bim file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedBuilder::allele_2`](struct.BedBuilder.html#method.allele_2).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let allele_2 = bed.allele_2()?;
    /// println!("{allele_2:?}"); // Outputs ndarray ["A", "C", "C", "G"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn allele_2(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_bim::<String>(
            self.metadata.allele_2.is_none(),
            MetadataFields::Allele2,
            "allele_2",
        )?;
        Ok(self.metadata.allele_2.as_ref().unwrap()) //unwrap always works because of lazy_bim
    }

    /// [`Metadata`](struct.Metadata.html) for this dataset, for example, the individual (sample) Ids.
    ///
    /// This returns a struct with 12 fields. Each field is a ndarray.
    /// The struct will always be new, but the 12 ndarrays will be
    /// shared with this [`Bed`](struct.Bed.html).
    ///
    /// If the needed, the metadata will be read from the .fam and/or .bim files.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, sample_bed_file};
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let metadata = bed.metadata()?;
    /// println!("{0:?}", metadata.iid()); // Outputs Some(["iid1", "iid2", "iid3"] ...)
    /// println!("{0:?}", metadata.sid()); // Outputs Some(["sid1", "sid2", "sid3", "sid4"] ...)
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    pub fn metadata(&mut self) -> Result<Metadata, Box<BedErrorPlus>> {
        self.fam()?;
        self.bim()?;
        Ok(self.metadata.clone())
    }

    /// Return the path of the .bed file.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Return the path of the .fam file.
    pub fn fam_path(&mut self) -> PathBuf {
        // We need to clone the path because self might mutate later
        if let Some(path) = &self.fam_path {
            path.clone()
        } else {
            let path = to_metadata_path(&self.path, &self.fam_path, "fam");
            self.fam_path = Some(path.clone());
            path
        }
    }

    /// Return the path of the .bim file.
    pub fn bim_path(&mut self) -> PathBuf {
        // We need to clone the path because self might mutate later
        if let Some(path) = &self.bim_path {
            path.clone()
        } else {
            let path = to_metadata_path(&self.path, &self.bim_path, "bim");
            self.bim_path = Some(path.clone());
            path
        }
    }

    /// Read genotype data.
    ///
    /// > Also see [`ReadOptions::builder`](struct.ReadOptions.html#method.builder) which supports selection and options.
    ///
    /// # Errors
    /// See [`BedError`](enum.BedError.html) and [`BedErrorPlus`](enum.BedErrorPlus.html)
    /// for all possible errors.
    ///
    /// # Examples
    /// Read all data in a .bed file.
    ///
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = bed.read::<f64>()?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1.0, 0.0, f64::NAN, 0.0],
    ///         [2.0, 0.0, f64::NAN, 2.0],
    ///         [0.0, 1.0, 2.0, 0.0]
    ///     ],
    /// );
    ///
    /// // Your output array can be f32, f64, or i8
    /// let val = bed.read::<i8>()?;
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1, 0, -127, 0],
    ///         [2, 0, -127, 2],
    ///         [0, 1, 2, 0]
    ///     ],
    /// );
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```    
    pub fn read<TVal: BedVal>(&mut self) -> Result<nd::Array2<TVal>, Box<BedErrorPlus>> {
        let read_options = ReadOptions::<TVal>::builder().build()?;
        self.read_with_options(&read_options)
    }

    /// Read genotype data with options, into a preallocated array.
    ///
    /// > Also see [`ReadOptionsBuilder::read_and_fill`](struct.ReadOptionsBuilder.html#method.read_and_fill).
    ///
    /// Note that options [`ReadOptions::f`](struct.ReadOptions.html#method.f),
    /// [`ReadOptions::c`](struct.ReadOptions.html#method.c), and [`ReadOptions::is_f`](struct.ReadOptionsBuilder.html#method.is_f)
    /// are ignored. Instead, the order of the preallocated array is used.
    ///
    /// # Errors
    /// See [`BedError`](enum.BedError.html) and [`BedErrorPlus`](enum.BedErrorPlus.html)
    /// for all possible errors.
    ///
    /// # Example
    ///
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// // Read the SNPs indexed by 2.
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let read_options = ReadOptions::builder().sid_index(2).build()?;
    /// let mut val = nd::Array2::<f64>::default((3, 1));
    /// bed.read_and_fill_with_options(&mut val.view_mut(), &read_options)?;
    ///
    /// assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```  
    pub fn read_and_fill_with_options<TVal: BedVal>(
        &mut self,
        val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.,
        read_options: &ReadOptions<TVal>,
    ) -> Result<(), Box<BedErrorPlus>> {
        let iid_count = self.iid_count()?;
        let sid_count = self.sid_count()?;

        let num_threads = compute_num_threads(read_options.num_threads)?;

        // If we already have a Vec<isize>, reference it. If we don't, create one and reference it.
        let iid_hold = Hold::new(&read_options.iid_index, iid_count)?;
        let iid_index = iid_hold.as_ref();
        let sid_hold = Hold::new(&read_options.sid_index, sid_count)?;
        let sid_index = sid_hold.as_ref();

        let dim = val.dim();
        if dim != (iid_index.len(), sid_index.len()) {
            return Err(Box::new(
                BedError::InvalidShape(iid_index.len(), sid_index.len(), dim.0, dim.1).into(),
            ));
        }

        read_no_alloc(
            &self.path,
            iid_count,
            sid_count,
            read_options.is_a1_counted,
            iid_index,
            sid_index,
            read_options.missing_value,
            num_threads,
            &mut val.view_mut(),
        )?;

        Ok(())
    }

    /// Read all genotype data into a preallocated array.
    ///
    /// > Also see [`ReadOptions::builder`](struct.ReadOptions.html#method.builder).
    ///
    /// # Errors
    /// See [`BedError`](enum.BedError.html) and [`BedErrorPlus`](enum.BedErrorPlus.html)
    /// for all possible errors.
    ///
    /// # Example
    ///
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let mut val = nd::Array2::<i8>::default(bed.dim()?);
    /// bed.read_and_fill(&mut val.view_mut())?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1, 0, -127, 0],
    ///         [2, 0, -127, 2],
    ///         [0, 1, 2, 0]
    ///     ],
    /// );
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn read_and_fill<TVal: BedVal>(
        &mut self,
        val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.,
    ) -> Result<(), Box<BedErrorPlus>> {
        let read_options = ReadOptions::<TVal>::builder().build()?;
        self.read_and_fill_with_options(val, &read_options)
    }

    /// Read genotype data with options.
    ///
    /// > Also see [`ReadOptions::builder`](struct.ReadOptions.html#method.builder).
    ///
    /// # Errors
    /// See [`BedError`](enum.BedError.html) and [`BedErrorPlus`](enum.BedErrorPlus.html)
    /// for all possible errors.
    ///
    /// # Example
    ///
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// // Read the SNPs indexed by 2.
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let read_options = ReadOptions::builder().sid_index(2).f64().build()?;
    /// let val = bed.read_with_options(&read_options)?;
    ///
    /// assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```  
    pub fn read_with_options<TVal: BedVal>(
        &mut self,
        read_options: &ReadOptions<TVal>,
    ) -> Result<nd::Array2<TVal>, Box<BedErrorPlus>> {
        let iid_count_in = self.iid_count()?;
        let sid_count_in = self.sid_count()?;
        let iid_count_out = read_options.iid_index.len(iid_count_in)?;
        let sid_count_out = read_options.sid_index.len(sid_count_in)?;
        let shape = ShapeBuilder::set_f((iid_count_out, sid_count_out), read_options.is_f);
        let mut val = nd::Array2::<TVal>::default(shape);

        self.read_and_fill_with_options(&mut val.view_mut(), read_options)?;

        Ok(val)
    }
    /// Write genotype data with default metadata.
    ///
    /// > Also see [`WriteOptions::builder`](struct.WriteOptions.html#method.builder), which supports metadata and options.
    ///
    /// # Errors
    /// See [`BedError`](enum.BedError.html) and [`BedErrorPlus`](enum.BedErrorPlus.html)
    /// for all possible errors.
    ///
    /// # Example
    /// In this example, write genotype data using default metadata.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    ///
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    ///
    /// let val = nd::array![[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]];
    /// Bed::write(&val, &output_file)?;
    ///
    /// // If we then read the new file and list the chromosome property,
    /// // it is an array of zeros, the default chromosome value.
    /// let mut bed2 = Bed::new(&output_file)?;
    /// println!("{:?}", bed2.chromosome()?); // Outputs ndarray ["0", "0", "0", "0"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn write<S: nd::Data<Elem = TVal>, TVal: BedVal>(
        val: &nd::ArrayBase<S, nd::Ix2>,
        path: &Path,
    ) -> Result<(), Box<BedErrorPlus>> {
        WriteOptions::builder(path).write(val)
    }

    /// Given an 2D array of genotype data and a [`WriteOptions`](struct.WriteOptionsBuilder.html), write to a .bed file.
    ///
    /// > Also see [`WriteOptionsBuilder::write`](struct.WriteOptionsBuilder.html#method.write), which creates
    /// > a [`WriteOptions`](struct.WriteOptionsBuilder.html) and writes to file in one step.
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    ///
    /// let val = nd::array![
    ///     [1.0, 0.0, f64::NAN, 0.0],
    ///     [2.0, 0.0, f64::NAN, 2.0],
    ///     [0.0, 1.0, 2.0, 0.0]
    /// ];
    ///
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .iid(["iid1", "iid2", "iid3"])
    ///     .sid(["sid1", "sid2", "sid3", "sid4"])
    ///     .build(3,4)?;
    ///
    /// Bed::write_with_options(&val, &write_options)?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn write_with_options<S, TVal>(
        val: &nd::ArrayBase<S, nd::Ix2>,
        write_options: &WriteOptions<TVal>,
    ) -> Result<(), Box<BedErrorPlus>>
    where
        S: nd::Data<Elem = TVal>,
        TVal: BedVal,
    {
        let (iid_count, sid_count) = val.dim();
        if iid_count != write_options.iid_count() {
            return Err(BedError::InconsistentCount(
                "iid".to_string(),
                write_options.iid_count(),
                iid_count,
            )
            .into());
        }
        if sid_count != write_options.sid_count() {
            return Err(BedError::InconsistentCount(
                "sid".to_string(),
                write_options.sid_count(),
                sid_count,
            )
            .into());
        }

        let num_threads = compute_num_threads(write_options.num_threads)?;
        write_val(
            &write_options.path,
            val,
            write_options.is_a1_counted,
            write_options.missing_value,
            num_threads,
        )?;

        if !write_options.skip_fam() {
            if let Err(e) = write_options.metadata.write_fam(write_options.fam_path()) {
                // Clean up the file
                let _ = fs::remove_file(&write_options.fam_path);
                return Err(e);
            }
        }

        if !write_options.skip_bim() {
            if let Err(e) = write_options.metadata.write_bim(write_options.bim_path()) {
                // Clean up the file
                let _ = fs::remove_file(&write_options.bim_path);
                return Err(e);
            }
        }

        Ok(())
    }

    fn unlazy_fam<T: FromStringArray<T>>(
        &mut self,
        is_none: bool,
        field_index: MetadataFields,
        name: &str,
    ) -> Result<(), Box<BedErrorPlus>> {
        if self.skip_set.contains(&field_index) {
            return Err(BedError::CannotUseSkippedMetadata(name.to_string()).into());
        }
        if is_none {
            self.fam()?
        }
        Ok(())
    }

    fn unlazy_bim<T: FromStringArray<T>>(
        &mut self,
        is_none: bool,
        field_index: MetadataFields,
        name: &str,
    ) -> Result<(), Box<BedErrorPlus>> {
        if self.skip_set.contains(&field_index) {
            return Err(BedError::CannotUseSkippedMetadata(name.to_string()).into());
        }
        if is_none {
            self.bim()?
        }
        Ok(())
    }

    fn fam(&mut self) -> Result<(), Box<BedErrorPlus>> {
        let fam_path = self.fam_path();

        let (metadata, count) = self.metadata.read_fam(fam_path, &self.skip_set)?;
        self.metadata = metadata;

        match self.iid_count {
            Some(iid_count) => {
                if iid_count != count {
                    return Err(
                        BedError::InconsistentCount("iid".to_string(), iid_count, count).into(),
                    );
                }
            }
            None => {
                self.iid_count = Some(count);
            }
        }
        Ok(())
    }

    fn bim(&mut self) -> Result<(), Box<BedErrorPlus>> {
        let bim_path = self.bim_path();

        let (metadata, count) = self.metadata.read_bim(bim_path, &self.skip_set)?;
        self.metadata = metadata;

        match self.sid_count {
            Some(sid_count) => {
                if sid_count != count {
                    return Err(
                        BedError::InconsistentCount("sid".to_string(), sid_count, count).into(),
                    );
                }
            }
            None => {
                self.sid_count = Some(count);
            }
        }
        Ok(())
    }
}

/// If we already have a Vec<isize> remember a reference to it.
/// If we don't, then create one.
enum Hold<'a> {
    Copy(Vec<isize>),
    Ref(&'a Vec<isize>),
}

impl Hold<'_> {
    fn new(index: &Index, count: usize) -> Result<Hold, Box<BedErrorPlus>> {
        let hold = if let Index::Vec(vec) = index {
            Hold::Ref(vec)
        } else {
            Hold::Copy(index.to_vec(count)?)
        };
        Ok(hold)
    }

    fn as_ref(&self) -> &Vec<isize> {
        match self {
            Hold::Ref(vec) => vec,
            Hold::Copy(ref vec) => vec,
        }
    }
}

fn compute_num_threads(option_num_threads: Option<usize>) -> Result<usize, Box<BedErrorPlus>> {
    let num_threads = if let Some(num_threads) = option_num_threads {
        num_threads
    } else if let Ok(num_threads) = env::var("BED_READER_NUM_THREADS") {
        num_threads.parse::<usize>()?
    } else if let Ok(num_threads) = env::var("NUM_THREADS") {
        num_threads.parse::<usize>()?
    } else {
        0
    };
    Ok(num_threads)
}

fn compute_max_concurrent_requests(
    option_max_concurrent_requests: Option<usize>,
) -> Result<usize, Box<BedErrorPlus>> {
    let max_concurrent_requests =
        if let Some(max_concurrent_requests) = option_max_concurrent_requests {
            max_concurrent_requests
        // } else if let Ok(num_threads) = env::var("BED_READER_NUM_THREADS") {
        //     num_threads.parse::<usize>()?
        // } else if let Ok(num_threads) = env::var("NUM_THREADS") {
        //     num_threads.parse::<usize>()?
        } else {
            10
        };
    Ok(max_concurrent_requests)
}

fn compute_max_chunk_size(
    option_max_chunk_size: Option<usize>,
) -> Result<usize, Box<BedErrorPlus>> {
    let max_chunk_size = if let Some(max_chunk_size) = option_max_chunk_size {
        max_chunk_size
    // } else if let Ok(num_threads) = env::var("BED_READER_NUM_THREADS") {
    //     num_threads.parse::<usize>()?
    // } else if let Ok(num_threads) = env::var("NUM_THREADS") {
    //     num_threads.parse::<usize>()?
    } else {
        10
    };
    Ok(max_chunk_size)
}

impl Index {
    // We can't define a 'From' because we want to add count at the last moment.
    // Later Would be nice to not always allocate a new vec, maybe with Rc<[T]>?
    // Even better would be to support an iterator from Index (an enum with fields).

    /// Turns an [`Index`](enum.Index.html) into a vector of usize indexes. Negative means count from end.
    pub fn to_vec(&self, count: usize) -> Result<Vec<isize>, Box<BedErrorPlus>> {
        let count_signed = count as isize;
        match self {
            Index::All => Ok((0..count_signed).collect()),
            Index::Vec(vec) => Ok(vec.to_vec()),
            Index::NDArrayBool(nd_array_bool) => {
                if nd_array_bool.len() != count {
                    return Err(
                        BedError::BoolArrayVectorWrongLength(count, nd_array_bool.len()).into(),
                    );
                }
                Ok(nd_array_bool
                    .iter()
                    .enumerate()
                    .filter(|(_, b)| **b)
                    .map(|(i, _)| i as isize)
                    .collect())
            }
            Index::NDSliceInfo(nd_slice_info) => {
                Ok(RangeNdSlice::new(nd_slice_info, count)?.to_vec())
            }
            Index::RangeAny(range_any) => {
                let range = range_any.to_range(count)?;
                Ok(range.map(|i| i as isize).collect::<Vec<isize>>())
            }
            Index::NDArray(nd_array) => Ok(nd_array.to_vec()),
            Index::One(one) => Ok(vec![*one]),
            Index::VecBool(vec_bool) => {
                if vec_bool.len() != count {
                    return Err(BedError::BoolArrayVectorWrongLength(count, vec_bool.len()).into());
                }
                Ok(vec_bool
                    .iter()
                    .enumerate()
                    .filter(|(_, b)| **b)
                    .map(|(i, _)| i as isize)
                    .collect())
            }
        }
    }
}

/// Type alias for 1-D slices of NDArrays.
pub type SliceInfo1 =
    nd::SliceInfo<[nd::SliceInfoElem; 1], nd::Dim<[usize; 1]>, nd::Dim<[usize; 1]>>;

/// A specification of which individuals (samples) or SNPs (variants) to read.
///
/// See the [Table of Index Expressions](index.html#index-expressions)
/// for a list of expressions for selecting individuals (sample)
/// and SNPs (variants).
///
/// By default, all individuals or SNPs are read.
/// The indices can be specified as:
///   * an index (negative numbers count from the end)
///   * a vector or ndarray of indices
///   * a Rust range (negatives not allowed)
///   * a vector or ndarray of booleans
///   * an ndarray slice (negative indexing and steps allowed)
///
/// # Examples
/// ```
/// use ndarray as nd;
/// use bed_reader::{Bed, ReadOptions, sample_bed_file};
/// use bed_reader::assert_eq_nan;
/// use ndarray::s;
///
/// let file_name = sample_bed_file("some_missing.bed")?;
/// let mut bed = Bed::new(file_name)?;
/// println!("{:?}", bed.dim()?); // prints (100, 100)
///
/// // Read all individuals and all SNPs
/// let val = ReadOptions::builder().f64().read(&mut bed)?;
/// assert!(val.dim() == (100, 100));
///
/// // Read the individual at index position 10 and all SNPs
/// let val = ReadOptions::builder().iid_index(10).f64().read(&mut bed)?;
/// assert!(val.dim() == (1, 100));
///
/// // Read the individuals at index positions 0,5, 1st-from-the-end and
/// // the SNP at index position 3
/// let val = ReadOptions::builder()
///     .iid_index(vec![0, 5, -1])
///     .sid_index(3)
///     .f64()
///     .read(&mut bed)?;
/// assert!(val.dim() == (3, 1));
/// // Repeat, but with an ndarray
/// let val = ReadOptions::builder()
///     .iid_index(nd::array![0, 5, -1])
///     .sid_index(3)
///     .f64()
///     .read(&mut bed)?;
/// assert!(val.dim() == (3, 1));
/// // Repeat, but with an Rust array
/// let val = ReadOptions::builder()
///     .iid_index([0, 5, -1])
///     .sid_index(3)
///     .f64()
///     .read(&mut bed)?;
/// assert!(val.dim() == (3, 1));

/// // Create a boolean ndarray identifying SNPs in chromosome 5,
/// // then select those SNPs.
/// let chrom_5 = bed.chromosome()?.map(|elem| elem == "5");
/// let val = ReadOptions::builder()
///     .sid_index(chrom_5)
///     .f64()
///     .read(&mut bed)?;
/// assert!(val.dim() == (100, 6));

/// // Use ndarray's slice macro, [`s!`](https://docs.rs/ndarray/latest/ndarray/macro.s.html),
/// // to select every 2nd individual and every 3rd SNP.
/// let val = ReadOptions::builder()
///     .iid_index(s![..;2])
///     .sid_index(s![..;3])
///     .f64()
///     .read(&mut bed)?;
/// assert!(val.dim() == (50, 34));
/// // Use ndarray's slice macro, [`s!`](https://docs.rs/ndarray/latest/ndarray/macro.s.html),
/// // to select the 10th-from-last individual to the last, in reverse order,
/// // and every 3rd SNP in reverse order.)
/// let val = ReadOptions::builder()
///     .iid_index(s![-10..;-1])
///     .sid_index(s![..;-3])
///     .f64()
///     .read(&mut bed)?;
/// assert!(val.dim() == (10, 34));
/// # use bed_reader::BedErrorPlus;
/// # Ok::<(), Box<BedErrorPlus>>(())
/// ```

#[derive(Debug, Clone)]
pub enum Index {
    // Could implement an enumerator, but it is complex and requires a 'match' on each next()
    //     https://stackoverflow.com/questions/65272613/how-to-implement-intoiterator-for-an-enum-of-iterable-variants
    #[allow(missing_docs)]
    All,
    #[allow(missing_docs)]
    One(isize),
    #[allow(missing_docs)]
    Vec(Vec<isize>),
    #[allow(missing_docs)]
    NDArray(nd::Array1<isize>),
    #[allow(missing_docs)]
    VecBool(Vec<bool>),
    #[allow(missing_docs)]
    NDArrayBool(nd::Array1<bool>),
    #[allow(missing_docs)]
    NDSliceInfo(SliceInfo1),
    #[allow(missing_docs)]
    RangeAny(RangeAny),
}

#[doc(hidden)]
/// Used internally to represent Rust ranges such as `0..10`, `..10`, etc.
#[derive(Debug, Clone)]
pub struct RangeAny {
    start: Option<usize>,
    end: Option<usize>,
}

impl RangeAny {
    fn new<T: RangeBounds<usize>>(range_thing: &T) -> RangeAny {
        let start_bound = range_thing.start_bound();
        let start = match start_bound {
            Bound::Included(&start) => Some(start),
            Bound::Excluded(&start) => Some(start + 1),
            Bound::Unbounded => None,
        };

        let end_bound = range_thing.end_bound();
        let end = match end_bound {
            Bound::Included(&end) => Some(end + 1),
            Bound::Excluded(&end) => Some(end),
            Bound::Unbounded => None,
        };
        RangeAny { start, end }
    }

    // https://stackoverflow.com/questions/55925523/array-cannot-be-indexed-by-rangefull
    fn to_range(&self, count: usize) -> Result<Range<usize>, Box<BedErrorPlus>> {
        let start = if let Some(start) = self.start {
            start
        } else {
            0
        };
        let end = if let Some(end) = self.end { end } else { count };
        if start > end {
            Err(BedError::StartGreaterThanEnd(start, end).into())
        } else {
            Ok(Range { start, end })
        }
    }

    fn len(&self, count: usize) -> Result<usize, Box<BedErrorPlus>> {
        let range = self.to_range(count)?;
        Ok(range.end - range.start)
    }

    fn is_empty(&self, count: usize) -> Result<bool, Box<BedErrorPlus>> {
        Ok(self.len(count)? == 0)
    }
}

#[doc(hidden)]
#[derive(Debug, Clone)]
/// Used internally to represent NDArray Slices such as s![..], s![0..;2], s![0..10;-1]
pub struct RangeNdSlice {
    start: usize,
    end: usize,
    step: usize,
    is_reversed: bool,
}

// https://www.geeksforgeeks.org/find-ceil-ab-without-using-ceil-function/
fn div_ceil(a: usize, b: usize) -> usize {
    (a + b - 1) / b
}

impl RangeNdSlice {
    fn len(&self) -> usize {
        if self.start > self.end {
            0
        } else {
            div_ceil(self.end - self.start, self.step)
        }
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    // https://docs.rs/ndarray/0.15.4/ndarray/struct.ArrayBase.html#slicing
    fn to_vec(&self) -> Vec<isize> {
        if self.start >= self.end {
            Vec::new()
        } else if !self.is_reversed {
            (self.start..self.end)
                .step_by(self.step)
                .map(|i| i as isize)
                .collect()
        } else {
            // https://docs.rs/ndarray/latest/ndarray/macro.s.html
            let size = self.len();
            let mut vec: Vec<isize> = Vec::<isize>::with_capacity(size);
            let mut i = self.end - 1;
            while i >= self.start {
                vec.push(i as isize);
                if i < self.step {
                    break;
                }
                i -= self.step;
            }
            vec
        }
    }

    fn new(nd_slice_info: &SliceInfo1, count: usize) -> Result<Self, Box<BedErrorPlus>> {
        //  self.to_vec(count).len(),
        // https://docs.rs/ndarray/0.15.4/ndarray/struct.ArrayBase.html#method.slice_collapse
        // Error in the following cases
        // * SliceInfo is not a 1-dimensional or is a NewAxis
        // * Step is 0
        // * Start is greater than count
        // * End is greater than count
        // As with ndarray, Start can be greater than End is allowed
        // and means the slice is empty.
        if nd_slice_info.in_ndim() != 1 || nd_slice_info.out_ndim() != 1 {
            return Err(BedError::NdSliceInfoNot1D.into());
        }

        let slice_info_elem = nd_slice_info[0];
        match slice_info_elem {
            nd::SliceInfoElem::Slice { start, end, step } => {
                // https://docs.rs/ndarray/0.15.4/ndarray/enum.SliceInfoElem.html
                // s![..], 0,None,1
                // s![a..b;2] a,b,2
                // s![a..;-1], from a to end in reverse order
                // start index; negative are counted from the back of the axis
                // end index; negative are counted from the back of the axis; when not present the default is the full length of the axis.
                // step size in elements; the default is 1, for every element.
                // A range with step size. end is an exclusive index. Negative start or end indexes are counted from the back of the axis. If end is None, the slice extends to the end of the axis.
                let step2: usize;
                let is_reverse2: bool;
                match step.cmp(&0) {
                    Ordering::Greater => {
                        step2 = step as usize;
                        is_reverse2 = false;
                    }
                    Ordering::Less => {
                        step2 = (-step) as usize;
                        is_reverse2 = true;
                    }
                    Ordering::Equal => {
                        return Err(BedError::StepZero.into());
                    }
                }

                let start2 = if start >= 0 {
                    let start3 = start as usize;
                    if start3 > count {
                        return Err(BedError::StartGreaterThanCount(start3, count).into());
                    } else {
                        start3
                    }
                } else {
                    let start3 = (-start) as usize;
                    if start3 > count {
                        return Err(BedError::StartGreaterThanCount(start3, count).into());
                    }
                    count - start3
                };

                let end2 = if let Some(end) = end {
                    if end >= 0 {
                        let end3 = end as usize;
                        if end3 > count {
                            return Err(BedError::EndGreaterThanCount(end3, count).into());
                        } else {
                            end3
                        }
                    } else {
                        let end3 = (-end) as usize;
                        if end3 > count {
                            return Err(BedError::EndGreaterThanCount(end3, count).into());
                        }
                        count - end3
                    }
                } else {
                    count
                };

                Ok(RangeNdSlice {
                    start: start2,
                    end: end2,
                    step: step2,
                    is_reversed: is_reverse2,
                })
            }
            nd::SliceInfoElem::Index(index) => Ok(RangeNdSlice {
                start: index as usize,
                end: index as usize + 1,
                step: 1,
                is_reversed: false,
            }),
            nd::SliceInfoElem::NewAxis => Err(BedError::NewAxis.into()),
        }
    }
}

impl Index {
    /// Returns the number of elements in an [`Index`](enum.Index.html).
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self, count: usize) -> Result<usize, Box<BedErrorPlus>> {
        match self {
            Index::All => Ok(count),
            Index::One(_) => Ok(1),
            Index::Vec(vec) => Ok(vec.len()),
            Index::NDArray(nd_array) => Ok(nd_array.len()),
            Index::VecBool(vec_bool) => Ok(vec_bool.iter().filter(|&b| *b).count()),
            Index::NDArrayBool(nd_array_bool) => Ok(nd_array_bool.iter().filter(|&b| *b).count()),
            Index::NDSliceInfo(nd_slice_info) => Ok(RangeNdSlice::new(nd_slice_info, count)?.len()),
            Index::RangeAny(range_any) => range_any.len(count),
        }
    }

    /// Returns true if the [`Index`](enum.Index.html) is empty.
    pub fn is_empty(&self, count: usize) -> Result<bool, Box<BedErrorPlus>> {
        match self {
            Index::All => Ok(count == 0),
            Index::One(_) => Ok(false),
            Index::Vec(vec) => Ok(vec.is_empty()),
            Index::NDArray(nd_array) => Ok(nd_array.is_empty()),
            Index::VecBool(vec_bool) => Ok(!vec_bool.iter().any(|&b| b)),
            Index::NDArrayBool(nd_array_bool) => Ok(!nd_array_bool.iter().any(|&b| b)),
            Index::NDSliceInfo(nd_slice_info) => {
                Ok(RangeNdSlice::new(nd_slice_info, count)?.is_empty())
            }
            Index::RangeAny(range_any) => range_any.is_empty(count),
        }
    }
}

impl From<SliceInfo1> for Index {
    fn from(slice_info: SliceInfo1) -> Index {
        Index::NDSliceInfo(slice_info)
    }
}
impl From<&SliceInfo1> for Index {
    fn from(slice_info: &SliceInfo1) -> Index {
        Index::NDSliceInfo(slice_info.to_owned())
    }
}

impl From<RangeFull> for Index {
    fn from(range_thing: RangeFull) -> Index {
        Index::RangeAny(RangeAny::new(&range_thing))
    }
}

impl From<&RangeFull> for Index {
    fn from(range_thing: &RangeFull) -> Index {
        Index::RangeAny(RangeAny::new(range_thing))
    }
}

impl From<Range<usize>> for Index {
    fn from(range_thing: Range<usize>) -> Index {
        Index::RangeAny(RangeAny::new(&range_thing))
    }
}

impl From<&Range<usize>> for Index {
    fn from(range_thing: &Range<usize>) -> Index {
        Index::RangeAny(RangeAny::new(range_thing))
    }
}

impl From<RangeFrom<usize>> for Index {
    fn from(range_thing: RangeFrom<usize>) -> Index {
        Index::RangeAny(RangeAny::new(&range_thing))
    }
}

impl From<&RangeFrom<usize>> for Index {
    fn from(range_thing: &RangeFrom<usize>) -> Index {
        Index::RangeAny(RangeAny::new(range_thing))
    }
}

impl From<RangeInclusive<usize>> for Index {
    fn from(range_thing: RangeInclusive<usize>) -> Index {
        Index::RangeAny(RangeAny::new(&range_thing))
    }
}

impl From<&RangeInclusive<usize>> for Index {
    fn from(range_thing: &RangeInclusive<usize>) -> Index {
        Index::RangeAny(RangeAny::new(range_thing))
    }
}

impl From<RangeTo<usize>> for Index {
    fn from(range_thing: RangeTo<usize>) -> Index {
        Index::RangeAny(RangeAny::new(&range_thing))
    }
}

impl From<&RangeTo<usize>> for Index {
    fn from(range_thing: &RangeTo<usize>) -> Index {
        Index::RangeAny(RangeAny::new(range_thing))
    }
}

impl From<RangeToInclusive<usize>> for Index {
    fn from(range_thing: RangeToInclusive<usize>) -> Index {
        Index::RangeAny(RangeAny::new(&range_thing))
    }
}

impl From<&RangeToInclusive<usize>> for Index {
    fn from(range_thing: &RangeToInclusive<usize>) -> Index {
        Index::RangeAny(RangeAny::new(range_thing))
    }
}

impl From<&[isize]> for Index {
    fn from(array: &[isize]) -> Index {
        Index::Vec(array.to_vec())
    }
}

impl<const N: usize> From<[isize; N]> for Index {
    fn from(array: [isize; N]) -> Index {
        Index::Vec(array.to_vec())
    }
}

impl<const N: usize> From<&[isize; N]> for Index {
    fn from(array: &[isize; N]) -> Index {
        Index::Vec(array.to_vec())
    }
}

impl From<&nd::ArrayView1<'_, isize>> for Index {
    fn from(view: &nd::ArrayView1<isize>) -> Index {
        Index::NDArray(view.to_owned())
    }
}

impl From<nd::ArrayView1<'_, isize>> for Index {
    fn from(view: nd::ArrayView1<isize>) -> Index {
        Index::NDArray(view.to_owned())
    }
}

impl From<Vec<isize>> for Index {
    fn from(vec: Vec<isize>) -> Index {
        Index::Vec(vec)
    }
}
impl From<&Vec<isize>> for Index {
    fn from(vec_ref: &Vec<isize>) -> Index {
        Index::Vec(vec_ref.clone())
    }
}

impl From<nd::ArrayView1<'_, bool>> for Index {
    fn from(view: nd::ArrayView1<bool>) -> Index {
        Index::NDArrayBool(view.to_owned())
    }
}

impl From<&nd::ArrayView1<'_, bool>> for Index {
    fn from(view: &nd::ArrayView1<bool>) -> Index {
        Index::NDArrayBool(view.to_owned())
    }
}

impl From<&Vec<bool>> for Index {
    fn from(vec_ref: &Vec<bool>) -> Index {
        Index::VecBool(vec_ref.clone())
    }
}

impl From<&[bool]> for Index {
    fn from(array: &[bool]) -> Index {
        Index::VecBool(array.to_vec())
    }
}

impl<const N: usize> From<[bool; N]> for Index {
    fn from(array: [bool; N]) -> Index {
        Index::VecBool(array.to_vec())
    }
}

impl<const N: usize> From<&[bool; N]> for Index {
    fn from(array: &[bool; N]) -> Index {
        Index::VecBool(array.to_vec())
    }
}

impl From<isize> for Index {
    fn from(one: isize) -> Index {
        Index::One(one)
    }
}
impl From<&isize> for Index {
    fn from(one: &isize) -> Index {
        Index::One(one.to_owned())
    }
}

impl From<nd::Array1<isize>> for Index {
    fn from(nd_array: nd::Array1<isize>) -> Index {
        Index::NDArray(nd_array)
    }
}

impl From<&nd::Array1<isize>> for Index {
    fn from(nd_array: &nd::Array1<isize>) -> Index {
        Index::NDArray(nd_array.to_owned())
    }
}

impl From<nd::Array1<bool>> for Index {
    fn from(nd_array_bool: nd::Array1<bool>) -> Index {
        Index::NDArrayBool(nd_array_bool)
    }
}

impl From<&nd::Array1<bool>> for Index {
    fn from(nd_array_bool: &nd::Array1<bool>) -> Index {
        Index::NDArrayBool(nd_array_bool.clone())
    }
}

impl From<Vec<bool>> for Index {
    fn from(vec_bool: Vec<bool>) -> Index {
        Index::VecBool(vec_bool)
    }
}

impl From<()> for Index {
    fn from(_: ()) -> Index {
        Index::All
    }
}

// See https://nullderef.com/blog/rust-parameters/

/// Represents options for reading genotype data from a PLINK .bed file.
///
/// Construct with [`ReadOptions::builder`](struct.ReadOptions.html#method.builder).
///
/// See the [Table of ReadOptions](index.html#readoptions)
/// for a list of the supported options.
/// See the [Table of Index Expressions](index.html#index-expressions)
/// for a list of expressions for selecting individuals (sample)
/// and SNPs (variants).
#[derive(Debug, Clone, Builder)]
#[builder(build_fn(error = "BedErrorPlus"))]
pub struct ReadOptions<TVal: BedVal> {
    /// Value to use for missing values (defaults to -127 or NaN)
    ///
    /// -127 is the default for i8 and NaN is the default for f32 and f64.
    ///
    /// In this example, the missing value is set to -1:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = ReadOptions::builder().missing_value(-1).i8().read(&mut bed)?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1, 0, -1, 0],
    ///         [2, 0, -1, 2],
    ///         [0, 1, 2, 0]
    ///     ],
    /// );
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[builder(default = "TVal::missing()")]
    missing_value: TVal,

    /// Select which individual (sample) values to read -- Defaults to all.
    ///
    /// Can select with a signed number, various lists of signed numbers,
    /// ranges, and various lists of booleans.
    ///
    /// See the [Table of Index Expressions](index.html#index-expressions)
    /// for a list of the supported index expressions.
    ///
    /// # Examples:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, assert_eq_nan, sample_bed_file};
    /// use ndarray::s;
    ///
    /// let file_name = sample_bed_file("some_missing.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    ///
    /// // Read the individual at index position 3
    ///
    /// let val = ReadOptions::builder()
    ///     .iid_index(3)
    ///     .f64()
    ///     .read(&mut bed)?;
    /// assert!(val.dim() == (1, 100));
    ///
    /// // Read the individuals at index positions 0, 5, and 1st-from-last.
    ///
    /// let val = ReadOptions::builder()
    ///     .iid_index([0, 5, -1])
    ///     .f64()
    ///     .read(&mut bed)?;
    ///
    /// assert!(val.dim() == (3, 100));
    ///
    /// // Read the individuals at index positions 20 (inclusive) to 30 (exclusive).
    ///
    /// let val = ReadOptions::builder()
    ///     .iid_index(20..30)
    ///     .f64()
    ///     .read(&mut bed)?;
    ///
    /// assert!(val.dim() == (10, 100));
    ///
    /// // Read the individuals at every 2nd index position.
    ///
    /// let val = ReadOptions::builder()
    ///     .iid_index(s![..;2])
    ///     .f64()
    ///     .read(&mut bed)?;
    ///
    /// assert!(val.dim() == (50, 100));
    ///
    /// // Read chromosome 5 of the female individuals.
    ///
    /// let female = bed.sex()?.map(|elem| *elem == 2);
    /// let chrom_5 = bed.chromosome()?.map(|elem| elem == "5");
    /// let val = ReadOptions::builder()
    ///     .iid_index(female)
    ///     .sid_index(chrom_5)
    ///     .f64()
    ///     .read(&mut bed)?;
    ///
    /// assert!(val.dim() == (50, 6));
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[builder(default = "Index::All")]
    #[builder(setter(into))]
    iid_index: Index,

    /// Select which SNPs (variant) values to read -- Defaults to all.
    ///
    /// Can select with a signed number, various lists of signed numbers,
    /// ranges, and various lists of booleans.
    ///
    /// See the [Table of Index Expressions](index.html#index-expressions)
    /// for a list of the supported index expressions.
    ///
    /// # Examples:
    /// ```
    /// use ndarray as nd;
    /// use ndarray::s;
    /// use bed_reader::{Bed, ReadOptions, assert_eq_nan, sample_bed_file};
    ///
    /// let file_name = sample_bed_file("some_missing.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    ///
    /// // Read the SNP at index position 3
    ///
    /// let val = ReadOptions::builder()
    ///     .sid_index(3)
    ///     .f64()
    ///     .read(&mut bed)?;
    /// assert!(val.dim() == (100, 1));
    ///
    /// // Read the SNPs at index positions 0, 5, and 1st-from-last.
    ///
    /// let val = ReadOptions::builder()
    ///     .sid_index([0, 5, -1])
    ///     .f64()
    ///     .read(&mut bed)?;
    ///
    /// assert!(val.dim() == (100, 3));
    ///
    /// // Read the SNPs at index positions 20 (inclusive) to 30 (exclusive).
    ///
    /// let val = ReadOptions::builder()
    ///     .sid_index(20..30)
    ///     .f64()
    ///     .read(&mut bed)?;
    ///
    /// assert!(val.dim() == (100, 10));
    ///
    /// // Read the SNPs at every 2nd index position.
    ///
    /// let val = ReadOptions::builder()
    ///     .sid_index(s![..;2])
    ///     .f64()
    ///     .read(&mut bed)?;
    ///
    /// assert!(val.dim() == (100, 50));
    ///
    /// // Read chromosome 5 of the female individuals.
    ///
    /// let female = bed.sex()?.map(|elem| *elem == 2);
    /// let chrom_5 = bed.chromosome()?.map(|elem| elem == "5");
    /// let val = ReadOptions::builder()
    ///     .iid_index(female)
    ///     .sid_index(chrom_5)
    ///     .f64()
    ///     .read(&mut bed)?;
    ///
    /// assert!(val.dim() == (50, 6));
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[builder(default = "Index::All")]
    #[builder(setter(into))]
    sid_index: Index,

    /// Sets if the order of the output array is Fortran-style -- Default is true.
    ///
    /// "Fortran order" is also called "column-major order" [Wikipedia](https://en.wikipedia.org/wiki/Row-_and_column-major_order).
    ///
    /// Also see [`f`](struct.ReadOptionsBuilder.html#method.f) and [`c`](struct.ReadOptionsBuilder.html#method.c).
    #[builder(default = "true")]
    is_f: bool,

    /// Sets if allele 1 is counted. Default is true.
    ///
    /// Also see [`count_a1`](struct.ReadOptionsBuilder.html#method.count_a1) and [`count_a2`](struct.ReadOptionsBuilder.html#method.count_a2).
    #[builder(default = "true")]
    is_a1_counted: bool,

    /// Number of threads to use (defaults to all processors)
    ///
    /// Can also be set with an environment variable.
    /// See [Environment Variables](index.html#environment-variables).
    ///
    /// In this example, we read using only one thread.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = ReadOptions::builder().num_threads(1).i8().read(&mut bed)?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1, 0, -127, 0],
    ///         [2, 0, -127, 2],
    ///         [0, 1, 2, 0]
    ///     ],
    /// );
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[builder(default, setter(strip_option))]
    num_threads: Option<usize>,

    /// cmk docs
    #[builder(default, setter(strip_option))]
    max_concurrent_requests: Option<usize>,

    /// cmk docs
    #[builder(default, setter(strip_option))]
    max_chunk_size: Option<usize>,
}

impl<TVal: BedVal> ReadOptions<TVal> {
    /// Read genotype data. Supports selection and options.
    ///
    /// > Also see [`Bed::read`](struct.Bed.html#method.read) (read without options).
    /// > To fill a preallocated ndarray, see [`ReadOptionsBuilder::read_and_fill`](struct.ReadOptionsBuilder.html#method.read_and_fill).
    ///
    /// See the [Table of ReadOptions](index.html#readoptions)
    /// for a list of the supported options.
    /// See the [Table of Index Expressions](index.html#index-expressions)
    /// for a list of expressions for selecting individuals (sample)
    /// and SNPs (variants).
    ///
    /// # Errors
    /// See [`BedError`](enum.BedError.html) and [`BedErrorPlus`](enum.BedErrorPlus.html)
    /// for all possible errors.
    ///
    /// # Examples
    ///
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// // Read all data from a .bed file into an ndarray of f64.
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = ReadOptions::builder().f64().read(&mut bed)?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1.0, 0.0, f64::NAN, 0.0],
    ///         [2.0, 0.0, f64::NAN, 2.0],
    ///         [0.0, 1.0, 2.0, 0.0]
    ///     ],
    /// );
    ///
    /// // Read the SNPs indexed by 2.
    /// let val = ReadOptions::builder().sid_index(2).f64().read(&mut bed)?;
    ///
    /// assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
    ///
    /// // Read the SNPs indexed by 2, 3, and 4th from last.
    /// let val = ReadOptions::builder()
    ///     .sid_index([2, 3, -4])
    ///     .f64()
    ///     .read(&mut bed)?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![[f64::NAN, 0.0, 1.0], [f64::NAN, 2.0, 2.0], [2.0, 0.0, 0.0]],
    /// );
    ///
    /// //  Read SNPs from 1 (inclusive) to 4 (exclusive).
    /// let val = ReadOptions::builder()
    ///     .sid_index(1..4)
    ///     .f64()
    ///     .read(&mut bed)?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![[0.0, f64::NAN, 0.0], [0.0, f64::NAN, 2.0], [1.0, 2.0, 0.0]],
    /// );
    ///
    /// // Print unique chrom values. Then, read all SNPs in chrom 5.
    /// use std::collections::HashSet;
    ///
    /// println!("{:?}", bed.chromosome()?.iter().collect::<HashSet<_>>());
    /// // This outputs: {"1", "5", "Y"}.
    /// let val = ReadOptions::builder()
    ///     .sid_index(bed.chromosome()?.map(|elem| elem == "5"))
    ///     .f64()
    ///     .read(&mut bed)?;
    ///
    /// assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
    ///
    /// // Read 1st individual (across all SNPs).
    /// let val = ReadOptions::builder().iid_index(0).f64().read(&mut bed)?;
    /// assert_eq_nan(&val, &nd::array![[1.0, 0.0, f64::NAN, 0.0]]);
    ///
    /// // Read every 2nd individual.
    /// use ndarray::s;
    ///
    /// let val = ReadOptions::builder()
    ///     .iid_index(s![..;2])
    ///     .f64()
    ///     .read(&mut bed)?;
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![[1.0, 0.0, f64::NAN, 0.0], [0.0, 1.0, 2.0, 0.0]],
    /// );
    ///
    /// // Read last and 2nd-to-last individuals and the last SNP
    /// let val = ReadOptions::builder()
    ///     .iid_index([-1,-2])
    ///     .sid_index(-1)
    ///     .f64()
    ///     .read(&mut bed)?;
    ///
    /// assert_eq_nan(&val, &nd::array![[0.0],[2.0]]);
    ///
    /// // The output array can be f32, f64, or i8
    /// let val = ReadOptions::builder().i8().read(&mut bed)?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1, 0, -127, 0],
    ///         [2, 0, -127, 2],
    ///         [0, 1, 2, 0]
    ///     ],
    /// );
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn builder() -> ReadOptionsBuilder<TVal> {
        ReadOptionsBuilder::default()
    }

    /// Value to be used for missing values (defaults to -127 or NaN).
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let read_options = ReadOptions::builder().sid_index([2, 3, 0]).i8().build()?;
    /// assert_eq!(read_options.missing_value(), -127);
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = bed.read_with_options(&read_options)?;

    /// assert_eq_nan(&val, &nd::array![[-127, 0, 1], [-127, 2, 2], [2, 0, 0]]);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn missing_value(&self) -> TVal {
        self.missing_value
    }

    /// Index of individuals (samples) to read (defaults to all).
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let read_options = ReadOptions::builder().sid_index([2, 3, 0]).i8().build()?;
    /// println!("{0:?}", read_options.iid_index()); // Outputs 'All'
    /// println!("{0:?}", read_options.sid_index()); // Outputs 'Vec([2, 3, 0])'
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = bed.read_with_options(&read_options)?;

    /// assert_eq_nan(&val, &nd::array![[-127, 0, 1], [-127, 2, 2], [2, 0, 0]]);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn iid_index(&self) -> &Index {
        &self.iid_index
    }

    /// Index of SNPs (variants) to read (defaults to all).
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let read_options = ReadOptions::builder().sid_index([2, 3, 0]).i8().build()?;
    /// println!("{0:?}", read_options.iid_index()); // Outputs 'All'
    /// println!("{0:?}", read_options.sid_index()); // Outputs 'Vec([2, 3, 0])'
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = bed.read_with_options(&read_options)?;

    /// assert_eq_nan(&val, &nd::array![[-127, 0, 1], [-127, 2, 2], [2, 0, 0]]);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn sid_index(&self) -> &Index {
        &self.sid_index
    }

    /// Is the order of the output array Fortran-style (defaults to true).
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let read_options = ReadOptions::builder().sid_index([2, 3, 0]).i8().build()?;
    /// assert_eq!(read_options.is_f(), true);
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = bed.read_with_options(&read_options)?;

    /// assert_eq_nan(&val, &nd::array![[-127, 0, 1], [-127, 2, 2], [2, 0, 0]]);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn is_f(&self) -> bool {
        self.is_f
    }

    /// If allele 1 will be counted (defaults to true).
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let read_options = ReadOptions::builder().sid_index([2, 3, 0]).i8().build()?;
    /// assert_eq!(read_options.is_a1_counted(), true);
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = bed.read_with_options(&read_options)?;

    /// assert_eq_nan(&val, &nd::array![[-127, 0, 1], [-127, 2, 2], [2, 0, 0]]);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn is_a1_counted(&self) -> bool {
        self.is_a1_counted
    }

    /// Number of threads to be used (`None` means set with
    /// [Environment Variables](index.html#environment-variables) or use all processors).
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let read_options = ReadOptions::builder().sid_index([2, 3, 0]).i8().build()?;
    /// assert_eq!(read_options.num_threads(), None);
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = bed.read_with_options(&read_options)?;

    /// assert_eq_nan(&val, &nd::array![[-127, 0, 1], [-127, 2, 2], [2, 0, 0]]);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn num_threads(&self) -> Option<usize> {
        self.num_threads
    }
}

impl<TVal: BedVal> ReadOptionsBuilder<TVal> {
    /// > See [`ReadOptions::builder`](struct.ReadOptions.html#method.builder) for details and examples.
    pub fn read(&self, bed: &mut Bed) -> Result<nd::Array2<TVal>, Box<BedErrorPlus>> {
        let read_options = self.build()?;
        bed.read_with_options(&read_options)
    }

    /// cmk
    pub async fn read_cloud<TArc>(
        &self,
        bed_cloud: &mut BedCloud<TArc>,
    ) -> Result<nd::Array2<TVal>, Box<BedErrorPlus>>
    where
        TArc: Clone + Deref + Send + Sync + 'static,
        TArc::Target: ObjectStore + Send + Sync,
    {
        let read_options = self.build()?;
        bed_cloud.read_with_options(&read_options).await
    }

    /// Read genotype data with options, into a preallocated array.
    ///
    /// > Also see [`Bed::read_and_fill`](struct.Bed.html#method.read_and_fill) and
    /// > [`Bed::read_and_fill_with_options`](struct.Bed.html#method.read_and_fill_with_options).
    ///
    /// Note that options [`ReadOptions::f`](struct.ReadOptions.html#method.f),
    /// [`ReadOptions::c`](struct.ReadOptions.html#method.c), and [`ReadOptions::is_f`](struct.ReadOptionsBuilder.html#method.is_f)
    /// are ignored. Instead, the order of the preallocated array is used.
    ///
    /// # Errors
    /// See [`BedError`](enum.BedError.html) and [`BedErrorPlus`](enum.BedErrorPlus.html)
    /// for all possible errors.
    ///
    /// # Example
    ///
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// // Read the SNPs indexed by 2.
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let mut val = nd::Array2::<f64>::default((3, 1));
    /// ReadOptions::builder()
    ///     .sid_index(2)
    ///     .read_and_fill(&mut bed, &mut val.view_mut())?;
    ///
    /// assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn read_and_fill(
        &self,
        bed: &mut Bed,
        val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.
    ) -> Result<(), Box<BedErrorPlus>> {
        let read_options = self.build()?;
        bed.read_and_fill_with_options(val, &read_options)
    }

    // cmk should there be read_and_fill_cloud?

    /// Order of the output array, Fortran-style (default)
    ///
    /// Also called "column-major order" [Wikipedia](https://en.wikipedia.org/wiki/Row-_and_column-major_order).
    ///
    /// Also see [`is_f`](struct.ReadOptionsBuilder.html#method.is_f) and [`c`](struct.ReadOptionsBuilder.html#method.c).
    pub fn f(&mut self) -> &mut Self {
        self.is_f(true);
        self
    }

    /// Order of the output array, C (default)
    ///
    /// Also called "row-major order" [Wikipedia](https://en.wikipedia.org/wiki/Row-_and_column-major_order).
    ///
    /// Also see [`is_f`](struct.ReadOptionsBuilder.html#method.is_f) and [`f`](struct.ReadOptionsBuilder.html#method.f).
    pub fn c(&mut self) -> &mut Self {
        self.is_f(false);
        self
    }

    /// Count the number allele 1 (default and PLINK standard).
    ///
    /// Also see [`is_a1_counted`](struct.ReadOptionsBuilder.html#method.is_a1_counted) and [`count_a2`](struct.ReadOptionsBuilder.html#method.count_a2).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = ReadOptions::builder().count_a1().i8().read(&mut bed)?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1, 0, -127, 0],
    ///         [2, 0, -127, 2],
    ///         [0, 1, 2, 0]
    ///     ],
    /// );
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn count_a1(&mut self) -> &mut Self {
        self.is_a1_counted = Some(true);
        self
    }

    /// Count the number allele 2.
    ///
    /// Also see [`is_a1_counted`](struct.ReadOptionsBuilder.html#method.is_a1_counted) and [`count_a1`](struct.ReadOptionsBuilder.html#method.count_a1).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = ReadOptions::builder().count_a2().i8().read(&mut bed)?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1, 2, -127, 2],
    ///         [0, 2, -127, 0],
    ///         [2, 1, 0, 2]
    ///     ],
    /// );
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn count_a2(&mut self) -> &mut Self {
        self.is_a1_counted = Some(false);
        self
    }
}

impl ReadOptionsBuilder<i8> {
    /// Output an ndarray of i8.
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = ReadOptions::builder().i8().read(&mut bed)?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1, 0, -127, 0],
    ///         [2, 0, -127, 2],
    ///         [0, 1, 2, 0]
    ///     ],
    /// );
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn i8(&mut self) -> &mut Self {
        self
    }
}

impl ReadOptionsBuilder<f32> {
    /// Output an ndarray of f32.
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = ReadOptions::builder().f32().read(&mut bed)?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1.0, 0.0, f32::NAN, 0.0],
    ///         [2.0, 0.0, f32::NAN, 2.0],
    ///         [0.0, 1.0, 2.0, 0.0]
    ///     ],
    /// );
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```    
    pub fn f32(&mut self) -> &mut Self {
        self
    }
}

impl ReadOptionsBuilder<f64> {
    /// Output an ndarray of f64.
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, ReadOptions, sample_bed_file};
    /// use bed_reader::assert_eq_nan;
    ///
    /// let file_name = sample_bed_file("small.bed")?;
    /// let mut bed = Bed::new(file_name)?;
    /// let val = ReadOptions::builder().f64().read(&mut bed)?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1.0, 0.0, f64::NAN, 0.0],
    ///         [2.0, 0.0, f64::NAN, 2.0],
    ///         [0.0, 1.0, 2.0, 0.0]
    ///     ],
    /// );
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```    
    pub fn f64(&mut self) -> &mut Self {
        self
    }
}

/// Represents options for writing genotype data and metadata to a PLINK .bed file.
///
/// Construct with [`WriteOptions::builder`](struct.WriteOptions.html#method.builder).
#[derive(Clone, Debug, Builder)]
#[builder(build_fn(skip))]
pub struct WriteOptions<TVal>
where
    TVal: BedVal,
{
    #[builder(setter(custom))]
    path: PathBuf,

    #[builder(setter(custom))]
    fam_path: PathBuf,

    #[builder(setter(custom))]
    bim_path: PathBuf,

    #[builder(setter(custom))]
    metadata: Metadata,

    #[builder(setter(custom), default = "true")]
    is_a1_counted: bool,

    #[builder(default, setter(custom))]
    num_threads: Option<usize>,

    #[builder(default = "TVal::missing()", setter(custom))]
    missing_value: TVal,

    #[builder(setter(custom), default = "false")]
    skip_fam: bool,

    #[builder(setter(custom), default = "false")]
    skip_bim: bool,
}

impl<TVal> WriteOptions<TVal>
where
    TVal: BedVal,
{
    /// Write values to a file in PLINK .bed format. Supports metadata and options.
    ///
    /// > Also see [`Bed::write`](struct.Bed.html#method.write), which does not support metadata or options.
    ///
    /// The options, [listed here](struct.WriteOptionsBuilder.html#implementations), can specify the:
    ///  * items of metadata, for example the individual ids or the SNP ids
    ///  * a non-default path for the .fam and/or .bim files
    ///  * a non-default value that represents missing data
    ///  * whether the first allele is counted (default) or the second
    ///  * number of threads to use for writing
    ///  * a [`Metadata`](struct.Metadata.html)
    ///
    /// # Examples
    /// In this example, all metadata is given one item at a time.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    ///
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let val = nd::array![
    ///     [1.0, 0.0, f64::NAN, 0.0],
    ///     [2.0, 0.0, f64::NAN, 2.0],
    ///     [0.0, 1.0, 2.0, 0.0]
    /// ];
    /// WriteOptions::builder(output_file)
    ///     .fid(["fid1", "fid1", "fid2"])
    ///     .iid(["iid1", "iid2", "iid3"])
    ///     .father(["iid23", "iid23", "iid22"])
    ///     .mother(["iid34", "iid34", "iid33"])
    ///     .sex([1, 2, 0])
    ///     .pheno(["red", "red", "blue"])
    ///     .chromosome(["1", "1", "5", "Y"])
    ///     .sid(["sid1", "sid2", "sid3", "sid4"])
    ///     .cm_position([100.4, 2000.5, 4000.7, 7000.9])
    ///     .bp_position([1, 100, 1000, 1004])
    ///     .allele_1(["A", "T", "A", "T"])
    ///     .allele_2(["A", "C", "C", "G"])
    ///     .write(&val)?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    /// Here, no metadata is given, so default values are assigned.
    /// If we then read the new file and list the chromosome property,
    /// it is an array of zeros, the default chromosome value.
    /// ```
    /// # use ndarray as nd;
    /// # use bed_reader::{Bed, WriteOptions};
    /// # let output_folder = temp_testdir::TempDir::default();
    /// let output_file2 = output_folder.join("small2.bed");
    /// let val = nd::array![[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]];
    ///
    /// WriteOptions::builder(&output_file2).write(&val)?;
    ///
    /// let mut bed2 = Bed::new(&output_file2)?;
    /// println!("{:?}", bed2.chromosome()?); // Outputs ndarray ["0", "0", "0", "0"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn builder(path: AnyPath) -> WriteOptionsBuilder<TVal> {
        WriteOptionsBuilder::new(path)
    }

    /// Family id of each of individual (sample). Defaults to "0"'s
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.fid()); // Outputs ndarray ["0", "0", "0"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn fid(&self) -> &nd::Array1<String> {
        // unwrap always works because the WriteOptions constructor fills all metadata.
        self.metadata.fid.as_ref().unwrap()
    }

    /// Individual id of each of individual (sample). Defaults to "iid1", "iid2" ...
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.iid()); // Outputs ndarray ["i1", "i2", "i3"]
    ///
    /// let val = nd::array![
    ///     [1.0, 0.0, f64::NAN, 0.0],
    ///     [2.0, 0.0, f64::NAN, 2.0],
    ///     [0.0, 1.0, 2.0, 0.0]
    /// ];
    /// Bed::write_with_options(&val, &write_options)?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn iid(&self) -> &nd::Array1<String> {
        // unwrap always works because the WriteOptions constructor fills all metadata.
        self.metadata.iid.as_ref().unwrap()
    }

    ///  Father id of each of individual (sample). Defaults to "0"'s
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::WriteOptions;
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.father()); // Outputs ndarray ["0", "0", "0"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn father(&self) -> &nd::Array1<String> {
        // unwrap always works because the WriteOptions constructor fills all metadata.
        self.metadata.father.as_ref().unwrap()
    }

    ///  Mother id of each of individual (sample). Defaults to "0"'s
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::WriteOptions;
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.mother()); // Outputs ndarray ["0", "0", "0"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn mother(&self) -> &nd::Array1<String> {
        // unwrap always works because the WriteOptions constructor fills all metadata.
        self.metadata.mother.as_ref().unwrap()
    }

    ///  Sex of each of individual (sample). Defaults to 0's
    ///
    /// 0 is unknown, 1 is male, 2 is female
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::WriteOptions;
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.sex()); // Outputs ndarray [0, 0, 0]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn sex(&self) -> &nd::Array1<i32> {
        // unwrap always works because the WriteOptions constructor fills all metadata.
        self.metadata.sex.as_ref().unwrap()
    }

    ///  Phenotype of each of individual (sample). Seldom used. Defaults to 0's
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::WriteOptions;
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.pheno()); // Outputs ndarray ["0", "0", "0"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn pheno(&self) -> &nd::Array1<String> {
        // unwrap always works because the WriteOptions constructor fills all metadata.
        self.metadata.pheno.as_ref().unwrap()
    }

    ///  Chromosome of each of SNP (variant). Defaults to "0"'s
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::WriteOptions;
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.chromosome()); // Outputs ndarray ["0", "0", "0", "0"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn chromosome(&self) -> &nd::Array1<String> {
        // unwrap always works because the WriteOptions constructor fills all metadata.
        self.metadata.chromosome.as_ref().unwrap()
    }

    ///  SNP id of each of SNP (variant). Defaults to "sid1", "sid2", ...
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.sid()); // Outputs ndarray ["s1", "s2", "s3", "s4"]
    ///
    /// let val = nd::array![
    ///     [1.0, 0.0, f64::NAN, 0.0],
    ///     [2.0, 0.0, f64::NAN, 2.0],
    ///     [0.0, 1.0, 2.0, 0.0]
    /// ];
    /// Bed::write_with_options(&val, &write_options)?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn sid(&self) -> &nd::Array1<String> {
        // unwrap always works because the WriteOptions constructor fills all metadata.
        self.metadata.sid.as_ref().unwrap()
    }

    /// Centimorgan position of each SNP (variant). Defaults to 0.0's.
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::WriteOptions;
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.cm_position()); // Outputs ndarray [0.0, 0.0, 0.0, 0.0]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn cm_position(&self) -> &nd::Array1<f32> {
        // unwrap always works because the WriteOptions constructor fills all metadata.
        self.metadata.cm_position.as_ref().unwrap()
    }

    /// Base-pair position of each SNP (variant). Defaults to 0's.
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.bp_position()); // Outputs ndarray [0, 0, 0, 0]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn bp_position(&self) -> &nd::Array1<i32> {
        // unwrap always works because the WriteOptions constructor fills all metadata.
        self.metadata.bp_position.as_ref().unwrap()
    }

    /// First allele of each SNP (variant). Defaults to "A1"
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.allele_1()); // Outputs ndarray ["A1", "A1", "A1", "A1"]
    /// println!("{0:?}", write_options.allele_2()); // Outputs ndarray ["A2", "A2", "A2", "A2"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn allele_1(&self) -> &nd::Array1<String> {
        // unwrap always works because the WriteOptions constructor fills all metadata.
        self.metadata.allele_1.as_ref().unwrap()
    }

    /// Second allele of each SNP (variant). Defaults to "A2"
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.allele_1()); // Outputs ndarray ["A1", "A1", "A1", "A1"]
    /// println!("{0:?}", write_options.allele_2()); // Outputs ndarray ["A2", "A2", "A2", "A2"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn allele_2(&self) -> &nd::Array1<String> {
        // unwrap always works because the WriteOptions constructor fills all metadata.
        self.metadata.allele_2.as_ref().unwrap()
    }

    /// [`Metadata`](struct.Metadata.html) for this [`WriteOptions`](struct.WriteOptions.html), for example, the individual (sample) Ids.
    ///
    /// This returns a struct with 12 fields. Each field is a ndarray.
    /// The struct will always be new, but the 12 ndarrays will be
    /// shared with this [`WriteOptions`](struct.WriteOptions.html).
    ///
    /// If the needed, default values will be used.
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// let metadata = write_options.metadata();
    /// println!("{0:?}", metadata.iid()); // Outputs optional ndarray Some(["i1", "i2", "i3"])
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn metadata(&self) -> Metadata {
        self.metadata.clone()
    }

    /// The number of individuals (samples)
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// assert_eq!(write_options.iid_count(), 3);
    /// assert_eq!(write_options.sid_count(), 4);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn iid_count(&self) -> usize {
        self.iid().len()
    }

    /// The number of SNPs (variants)
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// assert_eq!(write_options.iid_count(), 3);
    /// assert_eq!(write_options.sid_count(), 4);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn sid_count(&self) -> usize {
        self.sid().len()
    }

    /// Number of individuals (samples) and SNPs (variants)
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// assert_eq!(write_options.dim(), (3, 4));
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn dim(&self) -> (usize, usize) {
        (self.iid_count(), self.sid_count())
    }

    /// Path to .bed file.
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.path()); // Outputs "...small.bed"
    /// println!("{0:?}", write_options.fam_path()); // Outputs "...small.fam"
    /// println!("{0:?}", write_options.bim_path()); // Outputs "...small.bim"
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn path(&self) -> &PathBuf {
        &self.path
    }

    /// Path to .fam file.
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.path()); // Outputs "...small.bed"
    /// println!("{0:?}", write_options.fam_path()); // Outputs "...small.fam"
    /// println!("{0:?}", write_options.bim_path()); // Outputs "...small.bim"
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn fam_path(&self) -> &PathBuf {
        &self.fam_path
    }

    /// Path to .bim file.
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// println!("{0:?}", write_options.path()); // Outputs "...small.bed"
    /// println!("{0:?}", write_options.fam_path()); // Outputs "...small.fam"
    /// println!("{0:?}", write_options.bim_path()); // Outputs "...small.bim"
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn bim_path(&self) -> &PathBuf {
        &self.bim_path
    }

    /// If allele 1 will be counted (defaults to true).
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .i8()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// assert!(write_options.is_a1_counted());
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn is_a1_counted(&self) -> bool {
        self.is_a1_counted
    }

    /// Number of threads to be used (`None` means set with
    /// [Environment Variables](index.html#environment-variables) or use all processors).
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .i8()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// assert!(write_options.num_threads().is_none());
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn num_threads(&self) -> Option<usize> {
        self.num_threads
    }

    /// Value to be used for missing values (defaults to -127 or NaN).
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .i8()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    ///
    /// assert!(write_options.missing_value() == -127);
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn missing_value(&self) -> TVal {
        self.missing_value
    }

    /// If skipping writing .fam file.
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .i8()
    ///     .skip_fam()
    ///     .skip_bim()
    ///     .build(3, 4)?;
    /// assert!(write_options.skip_fam());
    /// assert!(write_options.skip_bim());
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn skip_fam(&self) -> bool {
        self.skip_fam
    }

    /// If skipping writing .bim file.
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .i8()
    ///     .skip_fam()
    ///     .skip_bim()
    ///     .build(3, 4)?;
    /// assert!(write_options.skip_fam());
    /// assert!(write_options.skip_bim());
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn skip_bim(&self) -> bool {
        self.skip_bim
    }
}

impl<TVal> WriteOptionsBuilder<TVal>
where
    TVal: BedVal,
{
    /// Creates a new [`WriteOptions`](struct.WriteOptions.html) with the options given and then writes a .bed (and .fam and .bim) file.
    ///
    /// See [`WriteOptions`](struct.WriteOptions.html) for details and examples.
    pub fn write<S: nd::Data<Elem = TVal>>(
        &mut self,
        val: &nd::ArrayBase<S, nd::Ix2>,
    ) -> Result<(), Box<BedErrorPlus>> {
        let (iid_count, sid_count) = val.dim();
        let write_options = self.build(iid_count, sid_count)?;
        Bed::write_with_options(val, &write_options)?;

        Ok(())
    }

    /// Set the family id (fid) values for each individual (sample).
    ///
    /// Defaults to zeros.
    ///
    /// > See [`WriteOptions`](struct.WriteOptions.html) for examples.
    ///
    #[anyinput]
    pub fn fid(mut self, fid: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because WriteOptionsBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_fid(fid);
        self
    }

    /// Set the individual id (iid) values for each individual (sample).
    ///
    /// Defaults to "iid1", "iid2", ...
    ///
    /// > See [`WriteOptions`](struct.WriteOptions.html) for examples.
    ///
    #[anyinput]
    pub fn iid(mut self, iid: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because WriteOptionsBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_iid(iid);
        self
    }

    /// Set the father id values for each individual (sample).
    ///
    /// Defaults to zeros.
    ///
    /// > See [`WriteOptions`](struct.WriteOptions.html) for examples.
    ///
    #[anyinput]
    pub fn father(mut self, father: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because WriteOptionsBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_father(father);
        self
    }

    /// Set the mother id values for each individual (sample).
    ///
    /// Defaults to zeros.
    ///
    /// > See [`WriteOptions`](struct.WriteOptions.html) for examples.
    ///
    #[anyinput]
    pub fn mother(mut self, mother: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because WriteOptionsBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_mother(mother);
        self
    }

    /// Set the sex for each individual (sample).
    ///
    /// 0 is unknown (default), 1 is male, 2 is female
    #[anyinput]
    pub fn sex(mut self, sex: AnyIter<i32>) -> Self {
        // Unwrap will always work because WriteOptionsBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_sex(sex);
        self
    }

    /// Set a phenotype for each individual (sample). Seldom used.
    ///
    /// Defaults to zeros.
    ///
    /// > See [`WriteOptions`](struct.WriteOptions.html) for examples.
    ///
    #[anyinput]
    pub fn pheno(mut self, pheno: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because WriteOptionsBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_pheno(pheno);
        self
    }

    /// Set the chromosome for each SNP (variant).
    ///
    /// Defaults to zeros.
    #[anyinput]
    pub fn chromosome(mut self, chromosome: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because WriteOptionsBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_chromosome(chromosome);
        self
    }

    /// Set the SNP id (sid) for each SNP (variant).
    ///
    /// Defaults to "sid1", "sid2", ...
    ///
    /// > See [`WriteOptions`](struct.WriteOptions.html) for examples.
    ///
    #[anyinput]
    pub fn sid(mut self, sid: AnyIter<AnyString>) -> Self {
        self.metadata.as_mut().unwrap().set_sid(sid);
        self
    }

    /// Set the centimorgan position for each SNP (variant).
    ///
    /// Defaults to zeros.
    #[anyinput]
    pub fn cm_position(mut self, cm_position: AnyIter<f32>) -> Self {
        // Unwrap will always work because WriteOptionsBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_cm_position(cm_position);
        self
    }

    /// Set the base-pair position for each SNP (variant).
    ///
    /// Defaults to zeros.
    ///
    /// > See [`WriteOptions`](struct.WriteOptions.html) for examples.
    ///
    #[anyinput]
    pub fn bp_position(mut self, bp_position: AnyIter<i32>) -> Self {
        // Unwrap will always work because WriteOptionsBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_bp_position(bp_position);
        self
    }

    /// Set the first allele for each SNP (variant).
    ///
    /// Defaults to "A1", A1" ...
    ///
    /// > See [`WriteOptions`](struct.WriteOptions.html) for examples.
    ///
    #[anyinput]
    pub fn allele_1(mut self, allele_1: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because WriteOptionsBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_allele_1(allele_1);
        self
    }

    /// Set the second allele for each SNP (variant).
    ///
    /// Defaults to "A2", A2" ...
    ///
    /// > See [`WriteOptions`](struct.WriteOptions.html) for examples.
    ///
    #[anyinput]
    pub fn allele_2(mut self, allele_2: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because WriteOptionsBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_allele_2(allele_2);
        self
    }

    /// Merge metadata from a [`Metadata`](struct.Metadata.html).
    ///
    /// If a field is set in both [`Metadata`](struct.Metadata.html)'s,
    /// it will be overridden.
    ///
    /// # Example
    ///
    /// Extract metadata from a file.
    /// Create a random file with the same metadata.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions, sample_bed_file};
    /// use ndarray_rand::{rand::prelude::StdRng, rand::SeedableRng, rand_distr::Uniform, RandomExt};
    ///
    /// let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
    /// let metadata = bed.metadata()?;
    /// let shape = bed.dim()?;
    ///
    /// let mut rng = StdRng::seed_from_u64(0);
    /// let val = nd::Array::random_using(shape, Uniform::from(-1..3), &mut rng);
    ///
    /// let temp_out = temp_testdir::TempDir::default();
    /// let output_file = temp_out.join("random.bed");
    /// WriteOptions::builder(output_file)
    ///     .metadata(&metadata)
    ///     .missing_value(-1)
    ///     .write(&val)?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn metadata(mut self, metadata: &Metadata) -> Self {
        self.metadata = Some(
            Metadata::builder()
                .metadata(&self.metadata.unwrap()) // Unwrap will always work because WriteOptionsBuilder starting with some metadata
                .metadata(metadata)
                .build_no_file_check() // Don't need to check consistent counts here. Builder will do it.
                .unwrap(), // Unwrap will always work nothing can go wrong
        );
        self
    }

    /// Set the path to the .fam file.
    ///
    /// If not set, the .fam file will be assumed
    /// to have the same name as the .bed file, but with the extension .fam.
    ///
    /// # Example:
    /// Write .bed, .fam, and .bim files with non-standard names.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::WriteOptions;
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.deb");
    /// let val = nd::array![[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]];

    /// WriteOptions::builder(output_file)
    ///     .fam_path(output_folder.join("small.maf"))
    ///     .bim_path(output_folder.join("small.mib"))
    ///     .write(&val)?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn fam_path(mut self, path: AnyPath) -> Self {
        self.fam_path = Some(path.to_owned());
        self
    }

    /// Set the path to the .bim file.
    ///
    /// If not set, the .bim file will be assumed
    /// to have the same name as the .bed file, but with the extension .bim.
    ///
    /// # Example:
    /// Write .bed, .fam, and .bim files with non-standard names.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.deb");
    /// let val = nd::array![[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]];

    /// WriteOptions::builder(output_file)
    ///     .fam_path(output_folder.join("small.maf"))
    ///     .bim_path(output_folder.join("small.mib"))
    ///     .write(&val)?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn bim_path(mut self, path: AnyPath) -> Self {
        self.bim_path = Some(path.to_owned());
        self
    }

    /// Value used for missing values (defaults to -127 or NaN)
    ///
    /// -127 is the default for i8 and NaN is the default for f32 and f64.
    ///
    /// # Example
    ///
    /// Extract metadata from a file.
    /// Create a random file with the same metadata.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions, sample_bed_file};
    /// use ndarray_rand::{rand::prelude::StdRng, rand::SeedableRng, rand_distr::Uniform, RandomExt};
    ///
    /// let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
    /// let metadata = bed.metadata()?;
    /// let shape = bed.dim()?;
    ///
    /// let mut rng = StdRng::seed_from_u64(0);
    /// let val = nd::Array::random_using(shape, Uniform::from(-1..3), &mut rng);
    ///
    /// let temp_out = temp_testdir::TempDir::default();
    /// let output_file = temp_out.join("random.bed");
    /// WriteOptions::builder(output_file)
    ///     .metadata(&metadata)
    ///     .missing_value(-1)
    ///     .write(&val)?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn missing_value(&mut self, missing_value: TVal) -> &mut Self {
        self.missing_value = Some(missing_value);
        self
    }

    /// Count the number allele 1 (default and PLINK standard).
    ///
    /// Also see [`is_a1_counted`](struct.WriteOptionsBuilder.html#method.is_a1_counted) and [`count_a2`](struct.WriteOptionsBuilder.html#method.count_a2).
    pub fn count_a1(&mut self) -> &mut Self {
        self.is_a1_counted = Some(true);
        self
    }

    /// Count the number allele 2.
    ///
    /// Also see [`is_a1_counted`](struct.WriteOptionsBuilder.html#method.is_a1_counted) and [`count_a1`](struct.WriteOptionsBuilder.html#method.count_a1).
    pub fn count_a2(&mut self) -> &mut Self {
        self.is_a1_counted = Some(false);
        self
    }

    /// Sets if allele 1 is counted. Default is true.
    ///
    /// Also see [`count_a1`](struct.WriteOptionsBuilder.html#method.count_a1) and [`count_a2`](struct.WriteOptionsBuilder.html#method.count_a2).    
    pub fn is_a1_counted(&mut self, is_a1_counted: bool) -> &mut Self {
        self.is_a1_counted = Some(is_a1_counted);
        self
    }

    /// Number of threads to use (defaults to all processors)
    ///
    /// Can also be set with an environment variable.
    /// See [Environment Variables](index.html#environment-variables).
    ///
    ///
    /// # Example:
    ///
    /// Write using only one thread.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::WriteOptions;
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let val = nd::array![[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]];

    /// WriteOptions::builder(output_file)
    ///     .num_threads(1)
    ///     .write(&val)?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn num_threads(&mut self, num_threads: usize) -> &mut Self {
        self.num_threads = Some(Some(num_threads));
        self
    }

    /// Skip writing .fam file.
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .i8()
    ///     .skip_fam()
    ///     .skip_bim()
    ///     .build(3, 4)?;
    /// assert!(write_options.skip_fam());
    /// assert!(write_options.skip_bim());
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn skip_fam(&mut self) -> &mut Self {
        self.skip_fam = Some(true);
        self
    }

    /// Skip writing .bim file.
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Bed, WriteOptions};
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .i8()
    ///     .skip_fam()
    ///     .skip_bim()
    ///     .build(3, 4)?;
    /// assert!(write_options.skip_fam());
    /// assert!(write_options.skip_bim());
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn skip_bim(&mut self) -> &mut Self {
        self.skip_bim = Some(true);
        self
    }

    /// Creates a new [`WriteOptions`](struct.WriteOptions.html) with the options given.
    ///
    /// > Also see [`WriteOptionsBuilder::write`](struct.WriteOptionsBuilder.html#method.write), which creates
    /// > a [`WriteOptions`](struct.WriteOptions.html) and writes to file in one step.
    ///
    /// # Example
    /// Create a new [`WriteOptions`](struct.WriteOptions.html) with some given values and some
    /// default values. Then use it to write a .bed file.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{WriteOptions, Bed};
    ///
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    /// let write_options = WriteOptions::builder(output_file)
    ///     .f64()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build(3, 4)?;
    /// println!("{0:?}", write_options.fid()); // Outputs ndarray ["0", "0", "0"]
    /// println!("{0:?}", write_options.iid()); // Outputs ndarray ["i1", "i2", "i3"]
    ///
    /// let val = nd::array![
    ///     [1.0, 0.0, f64::NAN, 0.0],
    ///     [2.0, 0.0, f64::NAN, 2.0],
    ///     [0.0, 1.0, 2.0, 0.0]
    /// ];
    /// Bed::write_with_options(&val, &write_options)?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn build(
        &self,
        iid_count: usize,
        sid_count: usize,
    ) -> Result<WriteOptions<TVal>, Box<BedErrorPlus>> {
        let path = if let Some(path) = self.path.as_ref() {
            path
        } else {
            return Err(UninitializedFieldError::new("path").into());
        };

        // unwrap always works because the metadata builder always initializes metadata
        let metadata = self.metadata.as_ref().unwrap();
        let metadata = metadata.fill(iid_count, sid_count)?;

        let write_options = WriteOptions {
            path: path.to_owned(),
            fam_path: to_metadata_path(path, &self.fam_path, "fam"),
            bim_path: to_metadata_path(path, &self.bim_path, "bim"),
            is_a1_counted: self.is_a1_counted.unwrap_or(true),
            num_threads: self.num_threads.unwrap_or(None),
            missing_value: self.missing_value.unwrap_or_else(|| TVal::missing()),
            skip_fam: self.skip_fam.unwrap_or(false),
            skip_bim: self.skip_bim.unwrap_or(false),

            metadata,
        };
        Ok(write_options)
    }

    #[anyinput]
    fn new(path: AnyPath) -> Self {
        Self {
            path: Some(path.to_owned()),
            fam_path: None,
            bim_path: None,

            metadata: Some(Metadata::new()),

            is_a1_counted: None,
            num_threads: None,
            missing_value: None,
            skip_fam: None,
            skip_bim: None,
        }
    }
}

trait FromStringArray<T> {
    fn from_string_array(
        string_array: nd::Array1<String>,
    ) -> Result<nd::Array1<Self>, Box<BedErrorPlus>>
    where
        Self: Sized;
}

impl FromStringArray<String> for String {
    fn from_string_array(
        string_array: nd::Array1<String>,
    ) -> Result<nd::Array1<String>, Box<BedErrorPlus>> {
        Ok(string_array)
    }
}

impl FromStringArray<f32> for f32 {
    fn from_string_array(
        string_array: nd::Array1<String>,
    ) -> Result<nd::Array1<f32>, Box<BedErrorPlus>> {
        let result = string_array
            .iter()
            .map(|s| s.parse::<f32>())
            .collect::<Result<nd::Array1<f32>, _>>();
        match result {
            Ok(array) => Ok(array),
            Err(e) => Err(Box::new(BedErrorPlus::ParseFloatError(e))),
        }
    }
}
impl FromStringArray<i32> for i32 {
    fn from_string_array(
        string_array: nd::Array1<String>,
    ) -> Result<nd::Array1<i32>, Box<BedErrorPlus>> {
        let result = string_array
            .iter()
            .map(|s| s.parse::<i32>())
            .collect::<Result<nd::Array1<i32>, _>>();
        match result {
            Ok(array) => Ok(array),
            Err(e) => Err(Box::new(BedErrorPlus::ParseIntError(e))),
        }
    }
}

/// Asserts two 2-D arrays are equal, treating NaNs as values.
///
/// # Example
/// ```
/// use std::f64::NAN;
/// use ndarray as nd;
/// use bed_reader::assert_eq_nan;
/// let val1 = nd::arr2(&[[1.0, 2.0], [3.0, NAN]]);
/// let val2 = nd::arr2(&[[1.0, 2.0], [3.0, NAN]]);
/// assert_eq_nan(&val1, &val2);
/// # use bed_reader::BedErrorPlus;
/// # Ok::<(), Box<BedErrorPlus>>(())
/// ```
pub fn assert_eq_nan<T: 'static + Copy + PartialEq + PartialOrd + Signed + From<i8>>(
    val: &nd::ArrayBase<nd::OwnedRepr<T>, nd::Dim<[usize; 2]>>,
    answer: &nd::ArrayBase<nd::OwnedRepr<T>, nd::Dim<[usize; 2]>>,
) {
    assert!(allclose::<T, T>(
        &val.view(),
        &answer.view(),
        0.into(),
        true
    ));
}

/// Asserts that a result is an error and that the error is of a given variant.
#[macro_export]
macro_rules! assert_error_variant {
    ($result:expr, $pattern:pat) => {
        match $result {
            Err(ref boxed_error) => match **boxed_error {
                $pattern => (),
                _ => panic!("test failure"),
            },
            _ => panic!("test failure"),
        }
    };
}

/// True if and only if two 2-D arrays are equal, within a given tolerance and possibly treating NaNs as values.
///
/// # Example
/// ```
/// use std::f64::NAN;
/// use ndarray as nd;
/// use bed_reader::allclose;
/// let val1 = nd::arr2(&[[1.0, 2.000000000001], [3.0, NAN]]);
/// let val2 = nd::arr2(&[[1.0, 2.0], [3.0, NAN]]);
/// assert!(allclose(&val1.view(), &val2.view(), 1e-08, true));
/// # use bed_reader::BedErrorPlus;
/// # Ok::<(), Box<BedErrorPlus>>(())
/// ```
pub fn allclose<
    T1: 'static + Copy + PartialEq + PartialOrd + Signed,
    T2: 'static + Copy + PartialEq + PartialOrd + Signed + Into<T1>,
>(
    val1: &nd::ArrayView2<'_, T1>,
    val2: &nd::ArrayView2<'_, T2>,
    atol: T1,
    equal_nan: bool,
) -> bool {
    assert!(val1.dim() == val2.dim());
    // Could be run in parallel

    nd::Zip::from(val1)
        .and(val2)
        .fold(true, |acc, ptr_a, ptr_b| -> bool {
            if !acc {
                return false;
            }
            // x != x is a generic nan check
            #[allow(clippy::eq_op)]
            let a_nan = *ptr_a != *ptr_a;
            #[allow(clippy::eq_op)]
            let b_nan = *ptr_b != *ptr_b;

            if a_nan || b_nan {
                if equal_nan {
                    a_nan == b_nan
                } else {
                    false
                }
            } else {
                let c: T1 = abs(*ptr_a - T2::into(*ptr_b));
                c <= atol
            }
        })
}

impl WriteOptionsBuilder<i8> {
    /// The input ndarray will be i8.
    pub fn i8(self) -> Self {
        self
    }
}

impl WriteOptionsBuilder<f32> {
    /// The input ndarray will be f32.
    pub fn f32(self) -> Self {
        self
    }
}

impl WriteOptionsBuilder<f64> {
    /// The input ndarray will be f64.
    pub fn f64(self) -> Self {
        self
    }
}

fn check_counts(
    count_vec: Vec<Option<usize>>,
    option_xid_count: &mut Option<usize>,
    prefix: &str,
) -> Result<(), Box<BedErrorPlus>> {
    for count in count_vec.into_iter().flatten() {
        if let Some(xid_count) = option_xid_count {
            if *xid_count != count {
                return Err(
                    BedError::InconsistentCount(prefix.to_string(), *xid_count, count).into(),
                );
            }
        } else {
            *option_xid_count = Some(count);
        }
    }

    Ok(())
}

// According to https://docs.rs/derive_builder/latest/derive_builder/
// "clone" is OK because "Luckily Rust is clever enough to optimize these
// clone-calls away in release builds for your every-day use cases.
// Thats quite a safe bet - we checked this for you. ;-)"
fn compute_field<T: Clone, F: Fn(usize) -> T>(
    field_name: &str,
    field: &mut Option<Rc<nd::Array1<T>>>,
    count: usize,
    lambda: F,
) -> Result<(), Box<BedErrorPlus>> {
    // let lambda = |_| "0".to_string();
    // let count = iid_count;
    // let field = &mut metadata.fid;

    if let Some(array) = field {
        if array.len() != count {
            return Err(
                BedError::InconsistentCount(field_name.to_string(), array.len(), count).into(),
            );
        }
    } else {
        let array = Rc::new((0..count).map(lambda).collect::<nd::Array1<T>>());
        *field = Some(array);
    }
    Ok(())
}

impl MetadataBuilder {
    /// Create [`Metadata`](struct.Metadata.html) from the builder.
    ///
    /// > See [`Metadata::builder()`](struct.Metadata.html#method.builder)
    pub fn build(&self) -> Result<Metadata, Box<BedErrorPlus>> {
        let metadata = self.build_no_file_check()?;

        metadata.check_counts(None, None)?;

        Ok(metadata)
    }

    /// Set the family id (fid) values.
    #[anyinput]
    pub fn fid(&mut self, fid: AnyIter<AnyString>) -> &mut Self {
        self.fid = Some(Some(Rc::new(fid.map(|s| s.as_ref().to_string()).collect())));
        self
    }

    /// Set the individual id (iid) values.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Metadata, assert_eq_nan};
    ///
    /// let metadata = Metadata::builder()
    ///    .iid(["sample1", "sample2", "sample3"])
    ///    .build()?;
    /// println!("{:?}", metadata.iid()); // Outputs ndarray Some(["sample1", "sample2", "sample3"])
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn iid(&mut self, iid: AnyIter<AnyString>) -> &mut Self {
        self.iid = Some(Some(Rc::new(iid.map(|s| s.as_ref().to_owned()).collect())));
        self
    }

    /// Set the father values.
    #[anyinput]
    pub fn father(&mut self, father: AnyIter<AnyString>) -> &mut Self {
        self.father = Some(Some(Rc::new(
            father.map(|s| s.as_ref().to_owned()).collect(),
        )));
        self
    }

    /// Override the mother values.
    #[anyinput]
    pub fn mother(&mut self, mother: AnyIter<AnyString>) -> &mut Self {
        self.mother = Some(Some(Rc::new(
            mother.map(|s| s.as_ref().to_owned()).collect(),
        )));
        self
    }

    /// Override the sex values.
    #[anyinput]
    pub fn sex(&mut self, sex: AnyIter<i32>) -> &mut Self {
        self.sex = Some(Some(Rc::new(sex.collect())));
        self
    }

    /// Override the phenotype values.
    #[anyinput]
    pub fn pheno(&mut self, pheno: AnyIter<AnyString>) -> &mut Self {
        self.pheno = Some(Some(Rc::new(
            pheno.map(|s| s.as_ref().to_owned()).collect(),
        )));
        self
    }

    /// Override the chromosome values.
    #[anyinput]
    pub fn chromosome(&mut self, chromosome: AnyIter<AnyString>) -> &mut Self {
        self.chromosome = Some(Some(Rc::new(
            chromosome.map(|s| s.as_ref().to_owned()).collect(),
        )));
        self
    }

    /// Override the SNP id (sid) values.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Metadata, assert_eq_nan};
    ///
    /// let metadata = Metadata::builder()
    ///    .sid(["SNP1", "SNP2", "SNP3", "SNP4"])
    ///    .build()?;
    /// println!("{:?}", metadata.sid()); // Outputs ndarray Some(["SNP1", "SNP2", "SNP3", "SNP4"])
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn sid(&mut self, sid: AnyIter<AnyString>) -> &mut Self {
        self.sid = Some(Some(Rc::new(
            sid.into_iter().map(|s| s.as_ref().to_owned()).collect(),
        )));
        self
    }

    /// Override the centimorgan position values.
    #[anyinput]
    pub fn cm_position(&mut self, cm_position: AnyIter<f32>) -> &mut Self {
        self.cm_position = Some(Some(Rc::new(cm_position.into_iter().collect())));
        self
    }

    /// Override the base-pair position values.
    #[anyinput]
    pub fn bp_position(&mut self, bp_position: AnyIter<i32>) -> &mut Self {
        self.bp_position = Some(Some(Rc::new(bp_position.into_iter().collect())));
        self
    }

    /// Override the allele 1 values.
    #[anyinput]
    pub fn allele_1(&mut self, allele_1: AnyIter<AnyString>) -> &mut Self {
        self.allele_1 = Some(Some(Rc::new(
            allele_1
                .into_iter()
                .map(|s| s.as_ref().to_owned())
                .collect(),
        )));
        self
    }

    /// Override the allele 2 values.
    #[anyinput]
    pub fn allele_2(&mut self, allele_2: AnyIter<AnyString>) -> &mut Self {
        self.allele_2 = Some(Some(Rc::new(
            allele_2
                .into_iter()
                .map(|s| s.as_ref().to_owned())
                .collect(),
        )));
        self
    }

    /// Merge metadata from a [`Metadata`](struct.Metadata.html).
    ///
    /// # Example
    ///
    /// In the example, we create a [`Metadata`](struct.Metadata.html) with iid
    /// and sid arrays. Next, we use another [`MetadataBuilder`](struct.MetadataBuilder.html) to set an fid array
    /// and an iid array. Then, we add the first [`Metadata`](struct.Metadata.html)
    /// to the [`MetadataBuilder`](struct.MetadataBuilder.html),
    /// overwriting iid and setting sid. Finally, we print these
    /// three arrays and chromosome. Chromosome is `None`.
    ///```
    /// use ndarray as nd;
    /// use bed_reader::Metadata;
    ///
    /// let metadata1 = Metadata::builder()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build()?;
    /// let metadata2 = Metadata::builder()
    ///     .fid(["f1", "f2", "f3"])
    ///     .iid(["x1", "x2", "x3"])
    ///     .metadata(&metadata1)
    ///     .build()?;
    ///
    /// println!("{0:?}", metadata2.fid()); // Outputs optional ndarray Some(["f1", "f2", "f3"]...)
    /// println!("{0:?}", metadata2.iid()); // Outputs optional ndarray Some(["i1", "i2", "i3"]...)
    /// println!("{0:?}", metadata2.sid()); // Outputs optional ndarray Some(["s1", "s2", "s3", "s4"]...)
    /// println!("{0:?}", metadata2.chromosome()); // Outputs None
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn metadata(&mut self, metadata: &Metadata) -> &mut Self {
        set_field(&metadata.fid, &mut self.fid);
        set_field(&metadata.iid, &mut self.iid);
        set_field(&metadata.father, &mut self.father);
        set_field(&metadata.mother, &mut self.mother);
        set_field(&metadata.sex, &mut self.sex);
        set_field(&metadata.pheno, &mut self.pheno);

        set_field(&metadata.chromosome, &mut self.chromosome);
        set_field(&metadata.sid, &mut self.sid);
        set_field(&metadata.cm_position, &mut self.cm_position);
        set_field(&metadata.bp_position, &mut self.bp_position);
        set_field(&metadata.allele_1, &mut self.allele_1);
        set_field(&metadata.allele_2, &mut self.allele_2);
        self
    }
}

impl Default for Metadata {
    fn default() -> Self {
        Self::new()
    }
}

impl Metadata {
    fn check_counts(
        &self,
        mut iid_count: Option<usize>,
        mut sid_count: Option<usize>,
    ) -> Result<(Option<usize>, Option<usize>), Box<BedErrorPlus>> {
        check_counts(
            vec![
                lazy_or_skip_count(&self.fid),
                lazy_or_skip_count(&self.iid),
                lazy_or_skip_count(&self.father),
                lazy_or_skip_count(&self.mother),
                lazy_or_skip_count(&self.sex),
                lazy_or_skip_count(&self.pheno),
            ],
            &mut iid_count,
            "iid",
        )?;
        check_counts(
            vec![
                lazy_or_skip_count(&self.chromosome),
                lazy_or_skip_count(&self.sid),
                lazy_or_skip_count(&self.cm_position),
                lazy_or_skip_count(&self.bp_position),
                lazy_or_skip_count(&self.allele_1),
                lazy_or_skip_count(&self.allele_2),
            ],
            &mut sid_count,
            "sid",
        )?;
        Ok((iid_count, sid_count))
    }

    /// Create a [`Metadata`](struct.Metadata.html) using a builder.
    ///
    /// # Example
    /// Create metadata.
    /// Create a random file with the metadata.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{Metadata, WriteOptions};
    /// use ndarray_rand::{rand::prelude::StdRng, rand::SeedableRng, rand_distr::Uniform, RandomExt};
    ///
    /// let metadata = Metadata::builder()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build()?;
    /// let mut rng = StdRng::seed_from_u64(0);
    /// let val = nd::Array::random_using((3, 4), Uniform::from(-1..3), &mut rng);

    /// let temp_out = temp_testdir::TempDir::default();
    /// let output_file = temp_out.join("random.bed");
    /// WriteOptions::builder(output_file)
    ///     .metadata(&metadata)
    ///     .missing_value(-1)
    ///     .write(&val)?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn builder() -> MetadataBuilder {
        MetadataBuilder::default()
    }

    /// Create an empty [`Metadata`](struct.Metadata.html).
    ///
    /// > See [`Metadata::builder()`](struct.Metadata.html#method.builder)
    pub fn new() -> Metadata {
        // Unwrap always works because an empty metadata builder always works.
        Metadata::builder().build().unwrap()
    }

    /// Optional family id of each of individual (sample)
    pub fn fid(&self) -> Option<&nd::Array1<String>> {
        option_rc_as_ref(&self.fid)
    }

    /// Optional individual id of each of individual (sample)
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::Metadata;
    /// let metadata = Metadata::builder().iid(["i1", "i2", "i3"]).build()?;
    /// println!("{0:?}", metadata.iid()); // Outputs optional ndarray Some(["i1", "i2", "i3"]...)
    /// println!("{0:?}", metadata.sid()); // Outputs None
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())    
    pub fn iid(&self) -> Option<&nd::Array1<String>> {
        option_rc_as_ref(&self.iid)
    }

    /// Optional father id of each of individual (sample)
    pub fn father(&self) -> Option<&nd::Array1<String>> {
        option_rc_as_ref(&self.father)
    }

    /// Optional mother id of each of individual (sample)
    pub fn mother(&self) -> Option<&nd::Array1<String>> {
        option_rc_as_ref(&self.mother)
    }

    /// Optional sex each of individual (sample)
    pub fn sex(&self) -> Option<&nd::Array1<i32>> {
        option_rc_as_ref(&self.sex)
    }

    /// Optional phenotype for each individual (seldom used)
    pub fn pheno(&self) -> Option<&nd::Array1<String>> {
        option_rc_as_ref(&self.pheno)
    }

    /// Optional chromosome of each SNP (variant)
    pub fn chromosome(&self) -> Option<&nd::Array1<String>> {
        option_rc_as_ref(&self.chromosome)
    }

    /// Optional SNP id of each SNP (variant)
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::Metadata;
    /// let metadata = Metadata::builder().iid(["i1", "i2", "i3"]).build()?;
    /// println!("{0:?}", metadata.iid()); // Outputs optional ndarray Some(["i1", "i2", "i3"]...)
    /// println!("{0:?}", metadata.sid()); // Outputs None
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())    
    pub fn sid(&self) -> Option<&nd::Array1<String>> {
        option_rc_as_ref(&self.sid)
    }

    /// Optional centimorgan position of each SNP (variant)
    pub fn cm_position(&self) -> Option<&nd::Array1<f32>> {
        option_rc_as_ref(&self.cm_position)
    }

    /// Optional base-pair position of each SNP (variant)
    pub fn bp_position(&self) -> Option<&nd::Array1<i32>> {
        option_rc_as_ref(&self.bp_position)
    }

    /// Optional first allele of each SNP (variant)
    pub fn allele_1(&self) -> Option<&nd::Array1<String>> {
        option_rc_as_ref(&self.allele_1)
    }

    /// Optional second allele of each SNP (variant)
    pub fn allele_2(&self) -> Option<&nd::Array1<String>> {
        option_rc_as_ref(&self.allele_2)
    }

    /// Create a new [`Metadata`](struct.Metadata.html) by filling in empty fields with a .fam file.
    ///
    /// # Example
    ///
    /// Read .fam and .bim information into a [`Metadata`](struct.Metadata.html).
    /// Do not skip any fields.
    /// ```
    /// use ndarray as nd;
    /// use std::collections::HashSet;
    /// use bed_reader::{Metadata, MetadataFields, sample_file};
    ///
    /// let skip_set = HashSet::<MetadataFields>::new();
    /// let metadata_empty = Metadata::new();
    /// let (metadata_fam, iid_count) =
    ///     metadata_empty.read_fam(sample_file("small.fam")?, &skip_set)?;
    /// let (metadata_bim, sid_count) =
    ///     metadata_fam.read_bim(sample_file("small.bim")?, &skip_set)?;
    /// assert_eq!(iid_count, 3);
    /// assert_eq!(sid_count, 4);
    /// println!("{0:?}", metadata_fam.iid()); // Outputs optional ndarray Some(["iid1", "iid2", "iid3"]...)
    /// println!("{0:?}", metadata_bim.sid()); // Outputs optional ndarray Some(["sid1", "sid2", "sid3", "sid4"]...)
    /// println!("{0:?}", metadata_bim.chromosome()); // Outputs optional ndarray Some(["1", "1", "5", "Y"]...)
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn read_fam(
        &self,
        path: AnyPath,
        skip_set: &HashSet<MetadataFields>,
    ) -> Result<(Metadata, usize), Box<BedErrorPlus>> {
        let mut field_vec: Vec<usize> = Vec::new();

        if self.fid.is_none() && !skip_set.contains(&MetadataFields::Fid) {
            field_vec.push(0);
        }
        if self.iid.is_none() && !skip_set.contains(&MetadataFields::Iid) {
            field_vec.push(1);
        }
        if self.father.is_none() && !skip_set.contains(&MetadataFields::Father) {
            field_vec.push(2);
        }
        if self.mother.is_none() && !skip_set.contains(&MetadataFields::Mother) {
            field_vec.push(3);
        }
        if self.sex.is_none() && !skip_set.contains(&MetadataFields::Sex) {
            field_vec.push(4);
        }
        if self.pheno.is_none() && !skip_set.contains(&MetadataFields::Pheno) {
            field_vec.push(5);
        }

        let (mut vec_of_vec, count) = self.read_fam_or_bim(&field_vec, true, path)?;

        let mut clone = self.clone();

        // unwraps are safe because we pop once for every push
        if clone.pheno.is_none() && !skip_set.contains(&MetadataFields::Pheno) {
            clone.pheno = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }
        if clone.sex.is_none() && !skip_set.contains(&MetadataFields::Sex) {
            let vec = vec_of_vec.pop().unwrap();
            let array = vec
                .iter()
                .map(|s| s.parse::<i32>())
                .collect::<Result<nd::Array1<i32>, _>>()?;
            clone.sex = Some(Rc::new(array));
        }
        if clone.mother.is_none() && !skip_set.contains(&MetadataFields::Mother) {
            clone.mother = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }
        if clone.father.is_none() && !skip_set.contains(&MetadataFields::Father) {
            clone.father = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }
        if clone.iid.is_none() && !skip_set.contains(&MetadataFields::Iid) {
            clone.iid = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }
        if clone.fid.is_none() && !skip_set.contains(&MetadataFields::Fid) {
            clone.fid = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }

        clone.check_counts(Some(count), None)?;

        Ok((clone, count))
    }

    /// cmk doc
    pub async fn read_cloud_fam<TArc>(
        &self,
        object_store: &TArc,
        path: &StorePath,
        skip_set: &HashSet<MetadataFields>,
    ) -> Result<(Metadata, usize), Box<BedErrorPlus>>
    where
        TArc: Clone + Deref + Send + Sync + 'static,
        TArc::Target: ObjectStore + Send + Sync,
    {
        let mut field_vec: Vec<usize> = Vec::new();

        if self.fid.is_none() && !skip_set.contains(&MetadataFields::Fid) {
            field_vec.push(0);
        }
        if self.iid.is_none() && !skip_set.contains(&MetadataFields::Iid) {
            field_vec.push(1);
        }
        if self.father.is_none() && !skip_set.contains(&MetadataFields::Father) {
            field_vec.push(2);
        }
        if self.mother.is_none() && !skip_set.contains(&MetadataFields::Mother) {
            field_vec.push(3);
        }
        if self.sex.is_none() && !skip_set.contains(&MetadataFields::Sex) {
            field_vec.push(4);
        }
        if self.pheno.is_none() && !skip_set.contains(&MetadataFields::Pheno) {
            field_vec.push(5);
        }

        let (mut vec_of_vec, count) = self
            .read_cloud_fam_or_bim(&field_vec, true, object_store, path)
            .await?;

        let mut clone = self.clone();

        // unwraps are safe because we pop once for every push
        if clone.pheno.is_none() && !skip_set.contains(&MetadataFields::Pheno) {
            clone.pheno = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }
        if clone.sex.is_none() && !skip_set.contains(&MetadataFields::Sex) {
            let vec = vec_of_vec.pop().unwrap();
            let array = vec
                .iter()
                .map(|s| s.parse::<i32>())
                .collect::<Result<nd::Array1<i32>, _>>()?;
            clone.sex = Some(Rc::new(array));
        }
        if clone.mother.is_none() && !skip_set.contains(&MetadataFields::Mother) {
            clone.mother = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }
        if clone.father.is_none() && !skip_set.contains(&MetadataFields::Father) {
            clone.father = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }
        if clone.iid.is_none() && !skip_set.contains(&MetadataFields::Iid) {
            clone.iid = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }
        if clone.fid.is_none() && !skip_set.contains(&MetadataFields::Fid) {
            clone.fid = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }

        clone.check_counts(Some(count), None)?;

        Ok((clone, count))
    }

    /// Create a new [`Metadata`](struct.Metadata.html) by filling in empty fields with a .bim file.
    ///
    /// # Example
    ///
    /// Read .fam and .bim information into a [`Metadata`](struct.Metadata.html).
    /// Do not skip any fields.
    /// ```
    /// use ndarray as nd;
    /// use std::collections::HashSet;
    /// use bed_reader::{Metadata, MetadataFields, sample_file};
    ///
    /// let skip_set = HashSet::<MetadataFields>::new();
    /// let metadata_empty = Metadata::new();
    /// let (metadata_fam, iid_count) =
    ///     metadata_empty.read_fam(sample_file("small.fam")?, &skip_set)?;
    /// let (metadata_bim, sid_count) =
    ///     metadata_fam.read_bim(sample_file("small.bim")?, &skip_set)?;
    /// assert_eq!(iid_count, 3);
    /// assert_eq!(sid_count, 4);
    /// println!("{0:?}", metadata_bim.iid()); // Outputs optional ndarray Some(["iid1", "iid2", "iid3"]...)
    /// println!("{0:?}", metadata_bim.sid()); // Outputs optional ndarray Some(["sid1", "sid2", "sid3", "sid4"]...)
    /// println!("{0:?}", metadata_bim.chromosome()); // Outputs optional ndarray Some(["1", "1", "5", "Y"]...)
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn read_bim(
        &self,
        path: AnyPath,
        skip_set: &HashSet<MetadataFields>,
    ) -> Result<(Metadata, usize), Box<BedErrorPlus>> {
        let mut field_vec: Vec<usize> = Vec::new();
        if self.chromosome.is_none() && !skip_set.contains(&MetadataFields::Chromosome) {
            field_vec.push(0);
        }
        if self.sid.is_none() && !skip_set.contains(&MetadataFields::Sid) {
            field_vec.push(1);
        }

        if self.cm_position.is_none() && !skip_set.contains(&MetadataFields::CmPosition) {
            field_vec.push(2);
        }
        if self.bp_position.is_none() && !skip_set.contains(&MetadataFields::BpPosition) {
            field_vec.push(3);
        }
        if self.allele_1.is_none() && !skip_set.contains(&MetadataFields::Allele1) {
            field_vec.push(4);
        }
        if self.allele_2.is_none() && !skip_set.contains(&MetadataFields::Allele2) {
            field_vec.push(5);
        }

        let mut clone = self.clone();
        let (mut vec_of_vec, count) = self.read_fam_or_bim(&field_vec, false, path)?;

        // unwraps are safe because we pop once for every push
        if clone.allele_2.is_none() && !skip_set.contains(&MetadataFields::Allele2) {
            clone.allele_2 = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }
        if clone.allele_1.is_none() && !skip_set.contains(&MetadataFields::Allele1) {
            clone.allele_1 = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }
        if clone.bp_position.is_none() && !skip_set.contains(&MetadataFields::BpPosition) {
            let vec = vec_of_vec.pop().unwrap();
            let array = vec
                .iter()
                .map(|s| s.parse::<i32>())
                .collect::<Result<nd::Array1<i32>, _>>()?;
            clone.bp_position = Some(Rc::new(array));
        }
        if clone.cm_position.is_none() && !skip_set.contains(&MetadataFields::CmPosition) {
            let vec = vec_of_vec.pop().unwrap();
            let array = vec
                .iter()
                .map(|s| s.parse::<f32>())
                .collect::<Result<nd::Array1<f32>, _>>()?;
            clone.cm_position = Some(Rc::new(array));
        }

        if clone.sid.is_none() && !skip_set.contains(&MetadataFields::Sid) {
            clone.sid = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }
        if clone.chromosome.is_none() && !skip_set.contains(&MetadataFields::Chromosome) {
            clone.chromosome = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }

        clone.check_counts(None, Some(count))?;

        Ok((clone, count))
    }

    /// cmk doc
    pub async fn read_cloud_bim<TArc>(
        &self,
        object_store: &TArc,
        path: &StorePath,
        skip_set: &HashSet<MetadataFields>,
    ) -> Result<(Metadata, usize), Box<BedErrorPlus>>
    where
        TArc: Clone + Deref + Send + Sync + 'static,
        TArc::Target: ObjectStore + Send + Sync,
    {
        let mut field_vec: Vec<usize> = Vec::new();
        if self.chromosome.is_none() && !skip_set.contains(&MetadataFields::Chromosome) {
            field_vec.push(0);
        }
        if self.sid.is_none() && !skip_set.contains(&MetadataFields::Sid) {
            field_vec.push(1);
        }

        if self.cm_position.is_none() && !skip_set.contains(&MetadataFields::CmPosition) {
            field_vec.push(2);
        }
        if self.bp_position.is_none() && !skip_set.contains(&MetadataFields::BpPosition) {
            field_vec.push(3);
        }
        if self.allele_1.is_none() && !skip_set.contains(&MetadataFields::Allele1) {
            field_vec.push(4);
        }
        if self.allele_2.is_none() && !skip_set.contains(&MetadataFields::Allele2) {
            field_vec.push(5);
        }

        let mut clone = self.clone();
        let (mut vec_of_vec, count) = self
            .read_cloud_fam_or_bim(&field_vec, false, object_store, path)
            .await?;

        // unwraps are safe because we pop once for every push
        if clone.allele_2.is_none() && !skip_set.contains(&MetadataFields::Allele2) {
            clone.allele_2 = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }
        if clone.allele_1.is_none() && !skip_set.contains(&MetadataFields::Allele1) {
            clone.allele_1 = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }
        if clone.bp_position.is_none() && !skip_set.contains(&MetadataFields::BpPosition) {
            let vec = vec_of_vec.pop().unwrap();
            let array = vec
                .iter()
                .map(|s| s.parse::<i32>())
                .collect::<Result<nd::Array1<i32>, _>>()?;
            clone.bp_position = Some(Rc::new(array));
        }
        if clone.cm_position.is_none() && !skip_set.contains(&MetadataFields::CmPosition) {
            let vec = vec_of_vec.pop().unwrap();
            let array = vec
                .iter()
                .map(|s| s.parse::<f32>())
                .collect::<Result<nd::Array1<f32>, _>>()?;
            clone.cm_position = Some(Rc::new(array));
        }

        if clone.sid.is_none() && !skip_set.contains(&MetadataFields::Sid) {
            clone.sid = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }
        if clone.chromosome.is_none() && !skip_set.contains(&MetadataFields::Chromosome) {
            clone.chromosome = Some(Rc::new(nd::Array::from_vec(vec_of_vec.pop().unwrap())));
        }

        clone.check_counts(None, Some(count))?;

        Ok((clone, count))
    }

    #[anyinput]
    fn read_fam_or_bim(
        &self,
        field_vec: &[usize],
        is_split_whitespace: bool,
        path: AnyPath,
    ) -> Result<(Vec<Vec<String>>, usize), Box<BedErrorPlus>> {
        let mut vec_of_vec = vec![vec![]; field_vec.len()];

        let file = File::open(path)?;

        let reader = BufReader::new(file);
        let mut count = 0;
        for line in reader.lines() {
            let line = line?;
            count += 1;

            let fields: Vec<&str> = if is_split_whitespace {
                line.split_whitespace().collect()
            } else {
                line.split('\t').collect()
            };

            if fields.len() != 6 {
                return Err(BedError::MetadataFieldCount(
                    6,
                    fields.len(),
                    path_ref_to_string(path),
                )
                .into());
            }

            let mut of_interest_count = 0;
            for (field_index, field) in fields.iter().enumerate() {
                if field_vec.contains(&field_index) {
                    vec_of_vec[of_interest_count].push(field.to_string());
                    of_interest_count += 1;
                }
            }
        }

        Ok((vec_of_vec, count))
    }

    async fn read_cloud_fam_or_bim<TArc>(
        &self,
        field_vec: &[usize],
        is_split_whitespace: bool,
        object_store: &TArc,
        path: &StorePath,
    ) -> Result<(Vec<Vec<String>>, usize), Box<BedErrorPlus>>
    where
        TArc: Clone + Deref + Send + Sync + 'static,
        TArc::Target: ObjectStore + Send + Sync,
    {
        let mut vec_of_vec = vec![vec![]; field_vec.len()];

        let stream = object_store
            .clone()
            .get(path)
            .await
            .map_err(BedErrorPlus::from)?
            .into_stream();

        let new_line_stream = newline_delimited_stream(stream);
        pin_mut!(new_line_stream);

        let mut count = 0;
        while let Some(chunk_result) = new_line_stream.next().await {
            let chunk = chunk_result.map_err(BedErrorPlus::from)?; // Handle the chunk result

            // Assuming chunk is a Bytes that can be converted to a string
            // Split the chunk into lines
            let lines = std::str::from_utf8(&chunk)
                .map_err(BedErrorPlus::from)?
                .split_terminator('\n');

            for line in lines {
                // let line = line?;
                count += 1;

                let fields: Vec<&str> = if is_split_whitespace {
                    line.split_whitespace().collect()
                } else {
                    line.split('\t').collect()
                };

                if fields.len() != 6 {
                    return Err(
                        BedError::MetadataFieldCount(6, fields.len(), path.to_string()).into(),
                    );
                }

                let mut of_interest_count = 0;
                for (field_index, field) in fields.iter().enumerate() {
                    if field_vec.contains(&field_index) {
                        vec_of_vec[of_interest_count].push(field.to_string());
                        of_interest_count += 1;
                    }
                }
            }
        }

        Ok((vec_of_vec, count))
    }

    fn is_some_fam(&self) -> bool {
        self.fid.is_some()
            && self.iid.is_some()
            && self.father.is_some()
            && self.mother.is_some()
            && self.sex.is_some()
            && self.pheno.is_some()
    }
    fn is_some_bim(&self) -> bool {
        self.chromosome.is_some()
            && self.sid.is_some()
            && self.cm_position.is_some()
            && self.bp_position.is_some()
            && self.allele_1.is_some()
            && self.allele_2.is_some()
    }

    /// Write the metadata related to individuals/samples to a .fam file.
    ///
    /// If any of the .fam metadata is not present, the function will return an error.
    ///
    /// # Example
    ///
    /// Create metadata with iid and sid arrays, then fill in the other
    /// fields with default arrays, finally write the .fam information
    /// to a file.
    ///```
    /// use ndarray as nd;
    /// use std::collections::HashSet;
    /// use bed_reader::Metadata;
    ///
    /// let metadata0 = Metadata::builder()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build()?;
    /// let metadata_filled = metadata0.fill(3, 4)?;

    /// let temp_out = temp_testdir::TempDir::default();
    /// let output_file = temp_out.join("no_bed.fam");
    /// metadata_filled.write_fam(output_file)?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn write_fam(&self, path: AnyPath) -> Result<(), Box<BedErrorPlus>> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        let mut result: Result<(), Box<BedErrorPlus>> = Ok(());

        if !self.is_some_fam() {
            return Err(BedError::MetadataMissingForWrite("fam".to_string()).into());
        }

        // 1st as_ref turns Option<Rc<Array>> into Option<&Rc<Array>>
        // unwrap always works because we checked that all the fields are present
        // 2nd as as_ref turns &Rc<Array> into &Array
        nd::azip!((fid in self.fid.as_ref().unwrap().as_ref(),
                   iid in self.iid.as_ref().unwrap().as_ref(),
                   father in self.father.as_ref().unwrap().as_ref(),
                   mother in self.mother.as_ref().unwrap().as_ref(),
                   sex in self.sex.as_ref().unwrap().as_ref(),
                   pheno in self.pheno.as_ref().unwrap().as_ref(),
                )
        {
            if result.is_ok() {
                if let Err(e) = writeln!(
                writer,
                "{} {} {} {} {} {}",
                *fid, *iid, *father, *mother, *sex, *pheno
            )
            {
            result = Err(Box::new(BedErrorPlus::IOError(e)));
            }
        }});
        result?;

        Ok(())
    }

    /// Write the metadata related to SNPs/variants to a .bim file.
    ///
    /// If any of the .bim metadata is not present, the function will return an error.
    ///
    /// # Example
    ///
    /// Create metadata with iid and sid arrays, then fill in the other
    /// fields with default arrays, finally write the .bim information
    /// to a file.
    ///```
    /// use ndarray as nd;
    /// use std::collections::HashSet;
    /// use bed_reader::Metadata;
    ///
    /// let metadata0 = Metadata::builder()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build()?;
    /// let metadata_filled = metadata0.fill(3, 4)?;

    /// let temp_out = temp_testdir::TempDir::default();
    /// let output_file = temp_out.join("no_bed.bim");
    /// metadata_filled.write_bim(output_file)?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn write_bim(&self, path: AnyPath) -> Result<(), Box<BedErrorPlus>> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        let mut result: Result<(), Box<BedErrorPlus>> = Ok(());

        if !self.is_some_bim() {
            return Err(BedError::MetadataMissingForWrite("bim".to_string()).into());
        }

        // 1st as_ref turns Option<Rc<Array>> into Option<&Rc<Array>>
        // unwrap always works because we checked that all the fields are present
        // 2nd as as_ref turns &Rc<Array> into &Array
        nd::azip!((
            chromosome in self.chromosome.as_ref().unwrap().as_ref(),
            sid in self.sid.as_ref().unwrap().as_ref(),
            cm_position in self.cm_position.as_ref().unwrap().as_ref(),
            bp_position in self.bp_position.as_ref().unwrap().as_ref(),
            allele_1 in self.allele_1.as_ref().unwrap().as_ref(),
            allele_2 in self.allele_2.as_ref().unwrap().as_ref(),
                )
        {
            if result.is_ok() {
                if let Err(e) = writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                *chromosome, *sid, *cm_position, *bp_position, *allele_1, *allele_2
                )
                {
                result = Err(Box::new(BedErrorPlus::IOError(e)));
                }
            }
        });
        result?;

        Ok(())
    }

    /// Create a new [`Metadata`](struct.Metadata.html) by filling in empty fields with default values.
    ///
    /// # Example
    /// ```
    /// use ndarray as nd;
    /// use std::collections::HashSet;
    /// use bed_reader::{Metadata, MetadataFields};
    ///
    /// let metadata0 = Metadata::builder()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build()?;
    /// let metadata_filled = metadata0.fill(3, 4)?;
    ///
    /// println!("{0:?}", metadata_filled.iid()); // Outputs optional ndarray Some(["i1", "i2", "i3"]...)
    /// println!("{0:?}", metadata_filled.sid()); // Outputs optional ndarray Some(["s1", "s2", "s3", "s4"]...)
    /// println!("{0:?}", metadata_filled.chromosome()); // Outputs optional ndarray Some(["0", "0", "0", "0"]...)
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    pub fn fill(&self, iid_count: usize, sid_count: usize) -> Result<Metadata, Box<BedErrorPlus>> {
        let mut metadata = self.clone();

        compute_field("fid", &mut metadata.fid, iid_count, |_| "0".to_string())?;
        compute_field("iid", &mut metadata.iid, iid_count, |i| {
            format!("iid{}", i + 1)
        })?;
        compute_field("father", &mut metadata.father, iid_count, |_| {
            "0".to_string()
        })?;
        compute_field("mother", &mut metadata.mother, iid_count, |_| {
            "0".to_string()
        })?;
        compute_field("sex", &mut metadata.sex, iid_count, |_| 0)?;
        compute_field("pheno", &mut metadata.pheno, iid_count, |_| "0".to_string())?;
        compute_field("chromosome", &mut metadata.chromosome, sid_count, |_| {
            "0".to_string()
        })?;
        compute_field("sid", &mut metadata.sid, sid_count, |i| {
            format!("sid{}", i + 1)
        })?;
        compute_field("cm_position", &mut metadata.cm_position, sid_count, |_| 0.0)?;
        compute_field("bp_position", &mut metadata.bp_position, sid_count, |_| 0)?;
        compute_field("allele_1", &mut metadata.allele_1, sid_count, |_| {
            "A1".to_string()
        })?;
        compute_field("allele_2", &mut metadata.allele_2, sid_count, |_| {
            "A2".to_string()
        })?;

        Ok(metadata)
    }

    #[anyinput]
    fn set_fid(&mut self, fid: AnyIter<AnyString>) -> &Self {
        self.fid = Some(Rc::new(
            fid.into_iter().map(|s| s.as_ref().to_owned()).collect(),
        ));
        self
    }

    #[anyinput]
    fn set_iid(&mut self, iid: AnyIter<AnyString>) -> &Self {
        self.iid = Some(Rc::new(
            iid.into_iter().map(|s| s.as_ref().to_owned()).collect(),
        ));
        self
    }

    #[anyinput]
    fn set_father(&mut self, father: AnyIter<AnyString>) -> &Self {
        self.father = Some(Rc::new(father.map(|s| s.as_ref().to_owned()).collect()));
        self
    }

    #[anyinput]
    fn set_mother(&mut self, mother: AnyIter<AnyString>) -> &Self {
        self.mother = Some(Rc::new(mother.map(|s| s.as_ref().to_owned()).collect()));
        self
    }

    #[anyinput]
    fn set_sex(&mut self, sex: AnyIter<i32>) -> &Self {
        self.sex = Some(Rc::new(sex.collect()));
        self
    }

    #[anyinput]
    fn set_pheno(&mut self, pheno: AnyIter<AnyString>) -> &Self {
        self.pheno = Some(Rc::new(pheno.map(|s| s.as_ref().to_owned()).collect()));
        self
    }

    #[anyinput]
    fn set_chromosome(&mut self, chromosome: AnyIter<AnyString>) -> &Self {
        self.chromosome = Some(Rc::new(chromosome.map(|s| s.as_ref().to_owned()).collect()));
        self
    }

    #[anyinput]
    fn set_sid(&mut self, sid: AnyIter<AnyString>) -> &Self {
        self.sid = Some(Rc::new(sid.map(|s| s.as_ref().to_owned()).collect()));
        self
    }

    #[anyinput]
    fn set_cm_position(&mut self, cm_position: AnyIter<f32>) -> &Self {
        self.cm_position = Some(Rc::new(cm_position.into_iter().collect()));
        self
    }

    #[anyinput]
    fn set_bp_position(&mut self, bp_position: AnyIter<i32>) -> &Self {
        self.bp_position = Some(Rc::new(bp_position.into_iter().collect()));
        self
    }

    #[anyinput]
    fn set_allele_1(&mut self, allele_1: AnyIter<AnyString>) -> &Self {
        self.allele_1 = Some(Rc::new(allele_1.map(|s| s.as_ref().to_owned()).collect()));
        self
    }

    #[anyinput]
    fn set_allele_2(&mut self, allele_2: AnyIter<AnyString>) -> &Self {
        self.allele_2 = Some(Rc::new(allele_2.map(|s| s.as_ref().to_owned()).collect()));
        self
    }
}

fn set_field<T>(
    field1: &Option<Rc<nd::Array1<T>>>,
    field2: &mut Option<Option<Rc<nd::Array1<T>>>>,
) {
    if let Some(array) = field1 {
        *field2 = Some(Some(array.clone()));
    }
}

fn option_rc_as_ref<T>(field: &Option<Rc<nd::Array1<T>>>) -> Option<&nd::Array1<T>> {
    match field {
        Some(array) => Some(array.as_ref()),
        None => None,
    }
}

#[allow(dead_code)]
fn matrix_subset_no_alloc<
    TIn: Copy + Default + Debug + Sync + Send + Sync + Sized,
    TOut: Copy + Default + Debug + Sync + Send + Sync + From<TIn>,
>(
    in_val: &nd::ArrayView3<'_, TIn>,
    iid_index: &[usize],
    sid_index: &[usize],
    out_val: &mut nd::ArrayViewMut3<'_, TOut>,
) -> Result<(), Box<BedErrorPlus>> {
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
        Ok(())
    } else {
        //If output is C-order, transpose input and output and recurse
        let in_val_t = in_val.view().permuted_axes([1, 0, 2]);
        let mut out_val_t = out_val.view_mut().permuted_axes([1, 0, 2]);
        matrix_subset_no_alloc(&in_val_t, sid_index, iid_index, &mut out_val_t)
    }
}

#[fetch_data::ctor]
static STATIC_FETCH_DATA: FetchData = FetchData::new(
    include_str!("../bed_reader/tests/registry.txt"),
    "https://raw.githubusercontent.com/fastlmm/bed-reader/rustybed/bed_reader/tests/data/",
    "BED_READER_DATA_DIR",
    "github.io",
    "fastlmm",
    "bed-reader",
);

/// Returns the local path to a sample .bed file. If necessary, the file will be downloaded.
///
/// The .fam and .bim files will also be downloaded, if they are not already present.
/// SHA256 hashes are used to verify that the files are correct.
/// The files will be in a directory determined by environment variable `BED_READER_DATA_DIR`.
/// If that environment variable is not set, a cache folder, appropriate to the OS, will be used.
#[anyinput]
pub fn sample_bed_file(bed_path: AnyPath) -> Result<PathBuf, Box<BedErrorPlus>> {
    let mut path_list: Vec<PathBuf> = Vec::new();
    for ext in ["bed", "bim", "fam"].iter() {
        let file_path = bed_path.with_extension(ext);
        path_list.push(file_path);
    }

    let vec = sample_files(path_list)?;
    assert!(vec.len() == 3);
    Ok(vec[0].clone())
}

/// Returns the local path to a sample file. If necessary, the file will be downloaded.
///
/// A SHA256 hash is used to verify that the file is correct.
/// The file will be in a directory determined by environment variable `BED_READER_DATA_DIR`.
/// If that environment variable is not set, a cache folder, appropriate to the OS, will be used.
#[anyinput]
pub fn sample_file(path: AnyPath) -> Result<PathBuf, Box<BedErrorPlus>> {
    match STATIC_FETCH_DATA.fetch_file(path) {
        Ok(path) => Ok(path),
        Err(e) => Err(e.into()),
    }
}

/// Returns the local paths to a list of files. If necessary, the files will be downloaded.
///
/// SHA256 hashes are used to verify that the files are correct.
/// The files will be in a directory determined by environment variable `BED_READER_DATA_DIR`.
/// If that environment variable is not set, a cache folder, appropriate to the OS, will be used.
#[anyinput]
pub fn sample_files(path_list: AnyIter<AnyPath>) -> Result<Vec<PathBuf>, Box<BedErrorPlus>>
where
{
    match STATIC_FETCH_DATA.fetch_files(path_list) {
        Ok(path) => Ok(path),
        Err(e) => Err(e.into()),
    }
}
