use anyinput::anyinput;
use bytes::Bytes;
use core::fmt;
use derive_builder::Builder;
use futures_util::{StreamExt, TryStreamExt};
use itertools::Itertools;
use nd::ShapeBuilder;
use ndarray as nd;
use object_store::delimited::newline_delimited_stream;
use object_store::local::LocalFileSystem;
use object_store::path::Path as StorePath;
use object_store::{GetOptions, GetResult};
use object_store::{ObjectMeta, ObjectStore};
use std::cmp::max;
use std::collections::HashSet;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use crate::{
    check_and_precompute_iid_index, compute_max_chunk_size, compute_max_concurrent_requests,
    set_up_two_bits_to_value, try_div_4, BedError, BedErrorPlus, BedVal, FromStringArray, Hold,
    Metadata, ReadOptions, BED_FILE_MAGIC1, BED_FILE_MAGIC2, STATIC_FETCH_DATA,
};
use crate::{MetadataFields, CB_HEADER_U64};

/// Represents a PLINK .bed file in the cloud that is open for reading genotype data and metadata.
///
/// Construct with [`BedCloud::new`](struct.BedCloud.html#method.new) or [`BedCloud::builder`](struct.BedCloud.html#method.builder).
///
/// # Example
///
/// Open a file for reading. Then, read the individual (sample) ids
/// and all the genotype data.
/// ```
/// use ndarray as nd;
/// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
/// use bed_reader::assert_eq_nan;
///
/// # Runtime::new().unwrap().block_on(async {
/// let object_path = sample_bed_object_path("small.bed")?;
/// let mut bed_cloud = BedCloud::new(object_path).await?;
/// println!("{:?}", bed_cloud.iid().await?); // Outputs ndarray ["iid1", "iid2", "iid3"]
/// let val = ReadOptions::builder().f64().read_cloud(&mut bed_cloud).await?;
///
/// assert_eq_nan(
///     &val,
///     &nd::array![
///         [1.0, 0.0, f64::NAN, 0.0],
///         [2.0, 0.0, f64::NAN, 2.0],
///         [0.0, 1.0, 2.0, 0.0]
///     ],
/// );
/// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
/// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
/// ```
#[derive(Clone, Debug, Builder)]
#[builder(build_fn(skip))]
pub struct BedCloud<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    #[builder(setter(custom))]
    object_path: ObjectPath<TObjectStore>,

    #[builder(setter(custom))]
    #[builder(default = "None")]
    fam_object_path: Option<ObjectPath<TObjectStore>>,

    #[builder(setter(custom))]
    #[builder(default = "None")]
    bim_object_path: Option<ObjectPath<TObjectStore>>,

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

// We need to define our own build_no_file_check
// because otherwise derive_builder (needlessly) requires ObjectStore: Clone
impl<TObjectStore> BedCloudBuilder<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn build_no_file_check(&self) -> Result<BedCloud<TObjectStore>, Box<BedErrorPlus>> {
        Ok(BedCloud {
            object_path: match self.object_path {
                Some(ref value) => Clone::clone(value),
                None => {
                    return Result::Err(Into::into(
                        ::derive_builder::UninitializedFieldError::from("object_path"),
                    ));
                }
            },
            fam_object_path: match self.fam_object_path {
                Some(ref value) => Clone::clone(value),
                None => None,
            },
            bim_object_path: match self.bim_object_path {
                Some(ref value) => Clone::clone(value),
                None => None,
            },
            is_checked_early: match self.is_checked_early {
                Some(ref value) => Clone::clone(value),
                None => true,
            },
            iid_count: match self.iid_count {
                Some(ref value) => Clone::clone(value),
                None => None,
            },
            sid_count: match self.sid_count {
                Some(ref value) => Clone::clone(value),
                None => None,
            },
            metadata: match self.metadata {
                Some(ref value) => Clone::clone(value),
                None => {
                    return Result::Err(Into::into(
                        ::derive_builder::UninitializedFieldError::from("metadata"),
                    ));
                }
            },
            skip_set: match self.skip_set {
                Some(ref value) => Clone::clone(value),
                None => {
                    return Result::Err(Into::into(
                        ::derive_builder::UninitializedFieldError::from("skip_set"),
                    ));
                }
            },
        })
    }
}

fn convert_negative_sid_index(
    in_sid_i_signed: isize,
    upper_sid_count: isize,
    lower_sid_count: isize,
) -> Result<u64, Box<BedErrorPlus>> {
    if (0..=upper_sid_count).contains(&in_sid_i_signed) {
        #[allow(clippy::cast_sign_loss)]
        Ok(in_sid_i_signed as u64)
    } else if (lower_sid_count..=-1).contains(&in_sid_i_signed) {
        #[allow(clippy::cast_sign_loss)]
        Ok((in_sid_i_signed - lower_sid_count) as u64) // cmk not sure about overflow
    } else {
        Err(Box::new(BedErrorPlus::BedError(BedError::SidIndexTooBig(
            in_sid_i_signed,
        ))))
    }
}

// cmk somehow we must only compile if size(usize) is 64 bits.

#[allow(clippy::too_many_arguments)]
#[allow(clippy::similar_names)]
async fn internal_read_no_alloc<TVal: BedVal, TObjectStore>(
    object_path: &ObjectPath<TObjectStore>,
    size: usize,
    in_iid_count: usize,
    in_sid_count: usize,
    is_a1_counted: bool,
    iid_index: &[isize],
    sid_index: &[isize],
    missing_value: TVal,
    max_concurrent_requests: usize,
    max_chunk_size: usize,
    out_val: &mut nd::ArrayViewMut2<'_, TVal>,
) -> Result<(), Box<BedErrorPlus>>
where
    TObjectStore: ObjectStore,
{
    // compute numbers outside of the loop
    let (in_iid_count_div4, in_iid_count_div4_u64) =
        check_file_length(in_iid_count, in_sid_count, size, object_path)?;
    let (i_div_4_array, i_mod_4_times_2_array) =
        check_and_precompute_iid_index(in_iid_count, iid_index)?;
    let chunk_size = max(1, max_chunk_size / in_iid_count_div4);
    let from_two_bits_to_value = set_up_two_bits_to_value(is_a1_counted, missing_value);
    let lower_sid_count = -(in_sid_count as isize);
    let upper_sid_count: isize = (in_sid_count as isize) - 1;

    // sid_index is a slice that tells us which columns to read from the (column-major) file.
    // out_val is a column-major array to fill the decode results.

    // For each chunk of columns to read ...

    let chunks = sid_index.iter().chunks(chunk_size);
    let iterator = chunks.into_iter().enumerate().map(|(chunk_index, chunk)| {
        let result = extract_ranges(
            chunk_size,
            chunk,
            chunk_index,
            upper_sid_count,
            lower_sid_count,
            in_iid_count_div4_u64,
        );
        async move {
            let (ranges, out_sid_i_vec) = result?;

            let vec_bytes = object_path.get_ranges(&ranges).await?;

            Result::<_, Box<BedErrorPlus>>::Ok((vec_bytes, out_sid_i_vec))
        }
    });

    let mut stream = futures_util::stream::iter(iterator).buffer_unordered(max_concurrent_requests);

    while let Some(result) = stream.next().await {
        let (vec_bytes, out_sid_i_vec) = result?;
        decode_bytes_into_columns(
            &vec_bytes,
            out_sid_i_vec,
            iid_index,
            &i_div_4_array,
            &i_mod_4_times_2_array,
            out_val,
            from_two_bits_to_value,
        );
    }

    Ok(())
}

#[inline]
#[allow(clippy::type_complexity)]
fn extract_ranges(
    chunk_size: usize,
    chunk: itertools::Chunk<'_, std::slice::Iter<'_, isize>>,
    chunk_index: usize,
    upper_sid_count: isize,
    lower_sid_count: isize,
    in_iid_count_div4_u64: u64,
) -> Result<(Vec<std::ops::Range<usize>>, Vec<usize>), Box<BedErrorPlus>> {
    let mut ranges = Vec::with_capacity(chunk_size);
    let mut out_sid_i_vec = Vec::with_capacity(chunk_size);
    for (inner_index, in_sid_i_signed) in chunk.enumerate() {
        let out_sid_i = chunk_index * chunk_size + inner_index;
        let in_sid_i =
            convert_negative_sid_index(*in_sid_i_signed, upper_sid_count, lower_sid_count)?;
        let pos: u64 = in_sid_i * in_iid_count_div4_u64 + CB_HEADER_U64; // "as" and math is safe because of early checks
        let range = pos as usize..(pos + in_iid_count_div4_u64) as usize;
        ranges.push(range);
        out_sid_i_vec.push(out_sid_i);
    }
    Ok((ranges, out_sid_i_vec))
}

#[inline]
fn decode_bytes_into_columns<TVal: BedVal>(
    bytes_slice: &[Bytes],
    out_sid_i_vec: Vec<usize>,
    iid_index: &[isize],
    i_div_4_array: &nd::prelude::ArrayBase<nd::OwnedRepr<usize>, nd::prelude::Dim<[usize; 1]>>,
    i_mod_4_times_2_array: &nd::prelude::ArrayBase<nd::OwnedRepr<u8>, nd::prelude::Dim<[usize; 1]>>,
    out_val: &mut nd::prelude::ArrayBase<nd::ViewRepr<&mut TVal>, nd::prelude::Dim<[usize; 2]>>,
    from_two_bits_to_value: [TVal; 4],
) {
    for (bytes, out_sid_i) in bytes_slice.iter().zip(out_sid_i_vec.into_iter()) {
        let mut col = out_val.column_mut(out_sid_i);
        // // cmk In parallel, decompress the iid info and put it in its column
        // // cmk .par_bridge() // This seems faster that parallel zip
        // .try_for_each(|(bytes_vector_result, mut col)| match bytes_vector_result {
        //     Err(e) => Err(e),
        //     Ok(bytes_vector) => {
        for out_iid_i in 0..iid_index.len() {
            let i_div_4 = i_div_4_array[out_iid_i];
            let i_mod_4_times_2: u8 = i_mod_4_times_2_array[out_iid_i];
            let encoded: u8 = bytes[i_div_4];
            let genotype_byte: u8 = (encoded >> i_mod_4_times_2) & 0x03;
            col[out_iid_i] = from_two_bits_to_value[genotype_byte as usize];
        }
    }
}

#[allow(clippy::similar_names)]
fn check_file_length<TObjectStore, I>(
    in_iid_count: usize,
    in_sid_count: usize,
    size: usize,
    object_path: I,
) -> Result<(usize, u64), Box<BedErrorPlus>>
where
    TObjectStore: ObjectStore,
    I: Into<ObjectPath<TObjectStore>>,
{
    let (in_iid_count_div4, in_iid_count_div4_u64) =
        try_div_4(in_iid_count, in_sid_count, CB_HEADER_U64)?;
    let file_len = size as u64;
    let file_len2 = in_iid_count_div4_u64 * (in_sid_count as u64) + CB_HEADER_U64;
    if file_len != file_len2 {
        return Err(Box::new(
            BedError::IllFormed(object_path.into().to_string()).into(),
        ));
    }
    Ok((in_iid_count_div4, in_iid_count_div4_u64))
}

#[inline]
#[allow(clippy::too_many_arguments)]
#[allow(clippy::similar_names)]
async fn read_no_alloc<TVal: BedVal, TObjectStore, I>(
    object_path: I,
    iid_count: usize,
    sid_count: usize,
    is_a1_counted: bool,
    iid_index: &[isize],
    sid_index: &[isize],
    missing_value: TVal,
    max_concurrent_requests: usize,
    max_chunk_size: usize,

    val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.
) -> Result<(), Box<BedErrorPlus>>
where
    TObjectStore: ObjectStore,
    I: Into<ObjectPath<TObjectStore>>,
{
    let object_path = object_path.into();
    let (size, bytes) = open_and_check(&object_path).await?;

    match bytes[2] {
        0 => {
            // We swap 'iid' and 'sid' and then reverse the axes.
            let mut val_t = val.view_mut().reversed_axes();

            internal_read_no_alloc(
                &object_path,
                size,
                sid_count,
                iid_count,
                is_a1_counted,
                sid_index,
                iid_index,
                missing_value,
                max_concurrent_requests,
                max_chunk_size,
                &mut val_t,
            )
            .await?;
        }
        1 => {
            internal_read_no_alloc(
                &object_path,
                size,
                iid_count,
                sid_count,
                is_a1_counted,
                iid_index,
                sid_index,
                missing_value,
                max_concurrent_requests,
                max_chunk_size,
                val,
            )
            .await?;
        }
        _ => return Err(Box::new(BedError::BadMode(object_path.to_string()).into())),
    };
    Ok(())
}

async fn open_and_check<TObjectStore, I>(
    object_path: I,
) -> Result<(usize, Bytes), Box<BedErrorPlus>>
where
    TObjectStore: ObjectStore,
    I: Into<ObjectPath<TObjectStore>>,
{
    let object_path = object_path.into();

    let object_store = object_path.object_store.clone();
    let path: &StorePath = &object_path.path;
    let object_meta = object_store
        .head(path)
        .await
        .map_err(|e| Box::new(BedErrorPlus::from(e)));
    let object_meta: ObjectMeta = object_meta?;
    let size: usize = object_meta.size;

    let get_options = GetOptions {
        range: Some(0..CB_HEADER_U64 as usize),
        ..Default::default()
    };
    let object_store = object_path.object_store.clone();
    let path: &StorePath = &object_path.path;
    let get_result = object_store
        .get_opts(path, get_options)
        .await
        .map_err(|e| Box::new(BedErrorPlus::from(e)));
    let get_result: GetResult = get_result?;
    let bytes = get_result.bytes().await.map_err(BedErrorPlus::from)?;

    if (BED_FILE_MAGIC1 != bytes[0]) || (BED_FILE_MAGIC2 != bytes[1]) {
        return Err(Box::new(
            BedError::IllFormed(object_path.to_string()).into(),
        ));
    }
    Ok((size, bytes))
}

impl<TObjectStore> BedCloudBuilder<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    // #[anyinput]
    fn new<I>(object_path: I) -> Self
    where
        I: Into<ObjectPath<TObjectStore>>,
    {
        Self {
            object_path: Some(object_path.into()),
            fam_object_path: None,
            bim_object_path: None,

            is_checked_early: None,
            iid_count: None,
            sid_count: None,

            metadata: Some(Metadata::new()),
            skip_set: Some(HashSet::new()),
        }
    }

    /// Create [`BedCloud`](struct.BedCloud.html) from the builder.
    ///
    /// > See [`BedCloud::builder`](struct.BedCloud.html#method.builder) for more details and examples.
    pub async fn build(&self) -> Result<BedCloud<TObjectStore>, Box<BedErrorPlus>> {
        let mut bed_cloud = self.build_no_file_check()?;

        // Unwrap is allowed becaue we can't construct BedCloudBuilder without object_path
        if bed_cloud.is_checked_early {
            let object_path = self.object_path.as_ref().unwrap().clone();
            open_and_check(object_path).await?;
        }

        (bed_cloud.iid_count, bed_cloud.sid_count) = bed_cloud
            .metadata
            .check_counts(bed_cloud.iid_count, bed_cloud.sid_count)?;

        Ok(bed_cloud)
    }

    /// cmk update docs
    /// Create [`BedCloud`](struct.BedCloud.html) from the builder.
    ///
    /// > See [`BedCloud::builder`](struct.BedCloud.html#method.builder) for more details and examples.
    pub fn build_no_check(&self) -> Result<BedCloud<TObjectStore>, Box<BedErrorPlus>> {
        let mut bed_cloud = self.build_no_file_check()?;

        (bed_cloud.iid_count, bed_cloud.sid_count) = bed_cloud
            .metadata
            .check_counts(bed_cloud.iid_count, bed_cloud.sid_count)?;

        Ok(bed_cloud)
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
        // Unwrap will always work because BedCloudBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_fid(fid);
        self
    }

    /// Override the individual id (iid) values found in the .fam file.
    ///
    /// By default, if iid values are needed and haven't already been found,
    /// they will be read from the .fam file.
    /// Providing them here avoids that file read and provides a way to give different values.
    /// ```
    /// # Runtime::new().unwrap().block_on(async {
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, assert_eq_nan, sample_bed_object_path};
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// use bed_reader::ReadOptions;
    ///
    /// let mut bed_cloud = BedCloud::builder(object_path)
    ///    .iid(["sample1", "sample2", "sample3"])
    ///    .build().await?;
    /// println!("{:?}", bed_cloud.iid().await?); // Outputs ndarray ["sample1", "sample2", "sample3"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    #[anyinput]
    pub fn iid(mut self, iid: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because BedCloudBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_iid(iid);
        self
    }

    /// Override the father values found in the .fam file.
    ///
    /// By default, if father values are needed and haven't already been found,
    /// they will be read from the .fam file.
    /// Providing them here avoids that file read and provides a way to gi&ve different values.
    #[anyinput]
    pub fn father(mut self, father: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because BedCloudBuilder starting with some metadata
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
        // Unwrap will always work because BedCloudBuilder starting with some metadata
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
        // Unwrap will always work because BedCloudBuilder starting with some metadata
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
        // Unwrap will always work because BedCloudBuilder starting with some metadata
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
        // Unwrap will always work because BedCloudBuilder starting with some metadata
        self.metadata.as_mut().unwrap().set_chromosome(chromosome);
        self
    }

    /// Override the SNP id (sid) values found in the .fam file.
    ///
    /// By default, if sid values are needed and haven't already been found,
    /// they will be read from the .bim file.
    /// Providing them here avoids that file read and provides a way to give different values.
    /// ```
    /// # Runtime::new().unwrap().block_on(async {
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_object_path};
    /// let object_path = sample_bed_object_path("small.bed")?;
    ///
    /// let mut bed_cloud = BedCloud::builder(object_path)
    ///    .sid(["SNP1", "SNP2", "SNP3", "SNP4"])
    ///    .build().await?;
    /// println!("{:?}", bed_cloud.sid().await?); // Outputs ndarray ["SNP1", "SNP2", "SNP3", "SNP4"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
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
        // Unwrap will always work because BedCloudBuilder starting with some metadata
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
        // Unwrap will always work because BedCloudBuilder starting with some metadata
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
        // Unwrap will always work because BedCloudBuilder starting with some metadata
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
        // Unwrap will always work because BedCloudBuilder starting with some metadata
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
    /// By default, when a [`BedCloud`](struct.BedCloud.html) struct is created, the .bed
    /// file header is checked. This stops that early check.
    /// ```
    /// # Runtime::new().unwrap().block_on(async {
    /// # use ndarray as nd;
    /// # use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_object_path};
    /// # let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::builder(&object_path).skip_early_check().build().await?;
    /// let val = bed_cloud.read::<f64>().await?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1.0, 0.0, f64::NAN, 0.0],
    ///         [2.0, 0.0, f64::NAN, 2.0],
    ///         [0.0, 1.0, 2.0, 0.0]
    ///     ],
    /// );
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    pub fn skip_early_check(mut self) -> Self {
        self.is_checked_early = Some(false);
        self
    }

    /// Set the cloud location of the .fam file.
    ///
    /// If not set, the .fam file will be assumed
    /// to have the same location as the .bed file, but with the extension .fam.
    ///
    /// # Example:
    /// Read .bed, .fam, and .bim files with non-standard names.
    /// ```
    /// use bed_reader::{BedCloud, ReadOptions, sample_object_paths};
    /// # Runtime::new().unwrap().block_on(async {
    /// let deb_maf_mib = sample_object_paths(["small.deb", "small.maf", "small.mib"])?;
    /// let mut bed_cloud = BedCloud::builder(&deb_maf_mib[0])
    ///    .fam_object_path(&deb_maf_mib[1])
    ///    .bim_object_path(&deb_maf_mib[2])
    ///    .build().await?;
    /// println!("{:?}", bed_cloud.iid().await?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// println!("{:?}", bed_cloud.sid().await?); // Outputs ndarray ["sid1", "sid2", "sid3", "sid4"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    pub fn fam_object_path<I>(mut self, object_path: I) -> Self
    where
        I: Into<ObjectPath<TObjectStore>>,
    {
        self.fam_object_path = Some(Some(object_path.into()));
        self
    }

    /// Set the cloud location of the .bim file.
    ///
    /// If not set, the .bim file will be assumed
    /// to have the same location as the .bed file, but with the extension .bim.
    ///
    /// # Example:
    /// Read .bed, .fam, and .bim files with non-standard names.
    /// ```
    /// # Runtime::new().unwrap().block_on(async {
    /// use bed_reader::{BedCloud, ReadOptions, sample_object_paths};
    /// let deb_maf_mib = sample_object_paths(["small.deb", "small.maf", "small.mib"])?;
    /// let mut bed_cloud = BedCloud::builder(&deb_maf_mib[0])
    ///    .fam_object_path(&deb_maf_mib[1])
    ///    .bim_object_path(&deb_maf_mib[2])
    ///    .build().await?;
    /// println!("{:?}", bed_cloud.iid().await?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// println!("{:?}", bed_cloud.sid().await?); // Outputs ndarray ["sid1", "sid2", "sid3", "sid4"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    // #[anyinput]
    pub fn bim_object_path<I>(mut self, object_path: I) -> Self
    where
        I: Into<ObjectPath<TObjectStore>>,
    {
        let object_path = object_path.into();
        self.bim_object_path = Some(Some(object_path));
        self
    }

    /// Don't read the fid information from the .fam file.
    ///
    /// By default, when the .fam is read, the fid (the family id) is recorded.
    /// This stops that recording. This is useful if the fid is not needed.
    /// Asking for the fid after skipping it results in an error.    
    pub fn skip_fid(mut self) -> Self {
        // Unwrap will always work because BedCloudBuilder starting with some skip_set
        self.skip_set.as_mut().unwrap().insert(MetadataFields::Fid);
        self
    }

    /// Don't read the iid information from the .fam file.
    ///
    /// By default, when the .fam is read, the iid (the individual id) is recorded.
    /// This stops that recording. This is useful if the iid is not needed.
    /// Asking for the iid after skipping it results in an error.
    pub fn skip_iid(mut self) -> Self {
        // Unwrap will always work because BedCloudBuilder starting with some skip_set
        self.skip_set.as_mut().unwrap().insert(MetadataFields::Iid);
        self
    }

    /// Don't read the father information from the .fam file.
    ///
    /// By default, when the .fam is read, the father id is recorded.
    /// This stops that recording. This is useful if the father id is not needed.
    /// Asking for the father id after skipping it results in an error.    
    pub fn skip_father(mut self) -> Self {
        // Unwrap will always work because BedCloudBuilder starting with some skip_set
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
        // Unwrap will always work because BedCloudBuilder starting with some skip_set
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
        // Unwrap will always work because BedCloudBuilder starting with some skip_set
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
        // Unwrap will always work because BedCloudBuilder starting with some skip_set
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
        // Unwrap will always work because BedCloudBuilder starting with some skip_set
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
        // Unwrap will always work because BedCloudBuilder starting with some skip_set
        self.skip_set.as_mut().unwrap().insert(MetadataFields::Sid);
        self
    }

    /// Don't read the centimorgan position information from the .bim file.
    ///
    /// By default, when the .bim is read, the cm position is recorded.
    /// This stops that recording. This is useful if the cm position is not needed.
    /// Asking for the cm position after skipping it results in an error.    
    pub fn skip_cm_position(mut self) -> Self {
        // Unwrap will always work because BedCloudBuilder starting with some skip_set
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
        // Unwrap will always work because BedCloudBuilder starting with some skip_set
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
        // Unwrap will always work because BedCloudBuilder starting with some skip_set
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
        // Unwrap will always work because BedCloudBuilder starting with some skip_set
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
    /// and sid arrays. Next, we use [`BedCloudBuilder`](struct.BedCloudBuilder.html) to override the fid array
    /// and an iid array. Then, we add the metadata to the [`BedCloudBuilder`](struct.BedCloudBuilder.html),
    /// overwriting iid (again) and overriding sid. Finally, we print these
    /// three arrays and chromosome. Chromosome was never overridden so
    /// it is read from the *.bim file.
    ///```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, Metadata, sample_bed_object_path};
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let metadata = Metadata::builder()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build()?;
    /// let mut bed_cloud = BedCloud::builder(object_path)
    ///     .fid(["f1", "f2", "f3"])
    ///     .iid(["x1", "x2", "x3"])
    ///     .metadata(&metadata)
    ///     .build().await?;
    /// println!("{0:?}", bed_cloud.fid().await?);  // Outputs ndarray ["f1", "f2", "f3"]
    /// println!("{0:?}", bed_cloud.iid().await?);  // Outputs ndarray ["i1", "i2", "i3"]
    /// println!("{0:?}", bed_cloud.sid().await?);  // Outputs ndarray ["s1", "s2", "s3", "s4"]
    /// println!("{0:?}", bed_cloud.chromosome().await?);  // Outputs ndarray ["1", "1", "5", "Y"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    pub fn metadata(mut self, metadata: &Metadata) -> Self {
        self.metadata = Some(
            Metadata::builder()
                .metadata(&self.metadata.unwrap()) // unwrap is ok because we know we have metadata
                .metadata(metadata) // consistent counts will be check later by the BedCloudBuilder
                .build_no_file_check()
                .unwrap(), // unwrap is ok because nothing can go wrong
        );

        self
    }
}

impl<TObjectStore> BedCloud<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    /// Attempts to open a PLINK .bed file in the cloud for reading. Supports options.
    ///
    /// > Also see [`BedCloud::new`](struct.BedCloud.html#method.new), which does not support options.
    ///
    /// The options, [listed here](struct.BedCloudBuilder.html#implementations), can:
    ///  * set the cloud location of the .fam and/or .bim file
    ///  * override some metadata, for example, replace the individual ids.
    ///  * set the number of individuals (samples) or SNPs (variants)
    ///  * control checking the validity of the .bed file's header
    ///  * skip reading selected metadata
    ///
    /// # Errors
    /// By default, this method will return an error if the file is missing or its header
    /// is ill-formed. It will also return an error if the options contradict each other.
    /// See [`BedError`](enum.BedError.html) and [`BedErrorPlus`](enum.BedErrorPlus.html)
    /// for all possible errors.
    ///
    /// # Examples
    /// List individual (sample) [`iid`](struct.BedCloud.html#method.iid) and
    /// SNP (variant) [`sid`](struct.BedCloud.html#method.sid),
    /// then [`read`](struct.BedCloud.html#method.read) the whole file.
    ///
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, assert_eq_nan, sample_bed_object_path};
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::builder(object_path).build().await?;
    /// println!("{:?}", bed_cloud.iid().await?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// println!("{:?}", bed_cloud.sid().await?); // Outputs ndarray ["snp1", "snp2", "snp3", "snp4"]
    /// let val = bed_cloud.read::<f64>().await?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1.0, 0.0, f64::NAN, 0.0],
    ///         [2.0, 0.0, f64::NAN, 2.0],
    ///         [0.0, 1.0, 2.0, 0.0]
    ///     ],
    /// );
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    ///
    /// Replace [`iid`](struct.BedCloud.html#method.iid).
    /// ```
    /// # Runtime::new().unwrap().block_on(async {
    /// # use ndarray as nd;
    /// # use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_object_path};
    /// # let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::builder(object_path)
    ///    .iid(["sample1", "sample2", "sample3"])
    ///    .build().await?;
    /// println!("{:?}", bed_cloud.iid().await?); // Outputs ndarray ["sample1", "sample2", "sample3"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    /// Give the number of individuals (samples) and SNPs (variants) so that the .fam and
    /// .bim files need never be opened. Use `.skip_early_check()` to avoid opening the
    /// .bed before the first read.
    /// ```
    /// # Runtime::new().unwrap().block_on(async {
    /// # use ndarray as nd;
    /// # use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_object_path};
    /// # let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::builder(&object_path)
    ///     .iid_count(3)
    ///     .sid_count(4)
    ///     .skip_early_check()
    ///     .build()
    ///     .await?;
    /// let val = bed_cloud.read::<f64>().await?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1.0, 0.0, f64::NAN, 0.0],
    ///         [2.0, 0.0, f64::NAN, 2.0],
    ///         [0.0, 1.0, 2.0, 0.0]
    ///     ],
    /// );
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    /// Mark some properties as "don’t read or offer".
    /// ```
    /// # Runtime::new().unwrap().block_on(async {
    /// # use ndarray as nd;
    /// # use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_object_path};
    /// # let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::builder(object_path)
    ///     .skip_father()
    ///     .skip_mother()
    ///     .skip_sex()
    ///     .skip_pheno()
    ///     .skip_allele_1()
    ///     .skip_allele_2()
    ///     .build().await?;
    /// println!("{:?}", bed_cloud.iid().await?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// bed_cloud.allele_2().await.expect_err("Can't be read");
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    ///
    pub fn builder<I>(object_path: I) -> BedCloudBuilder<TObjectStore>
    where
        I: Into<ObjectPath<TObjectStore>>,
    {
        let object_path = object_path.into();
        BedCloudBuilder::new(object_path)
    }

    /// Attempts to open a PLINK .bed file in the cloud for reading. Does not support options.
    ///
    /// > Also see [`BedCloud::builder`](struct.BedCloud.html#method.builder), which does support options,
    /// > including [`skip_early_check`](struct.BedCloudBuilder.html#method.skip_early_check).
    ///
    /// # Errors
    /// By default, this method will return an error if the file is missing or its header
    /// is ill-formed. See [`BedError`](enum.BedError.html) and [`BedErrorPlus`](enum.BedErrorPlus.html)
    /// for all possible errors.
    ///
    /// # Examples
    /// List individual (sample) [`iid`](struct.BedCloud.html#method.iid) and
    /// SNP (variant) [`sid`](struct.BedCloud.html#method.sid),
    /// then [`read`](struct.BedCloud.html#method.read) the whole file.
    ///
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, assert_eq_nan, sample_bed_object_path};
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// println!("{:?}", bed_cloud.iid().await?); // Outputs ndarray: ["iid1", "iid2", "iid3"]
    /// println!("{:?}", bed_cloud.sid().await?); // Outputs ndarray: ["sid1", "sid2", "sid3", "sid4"]
    /// let val = bed_cloud.read::<f64>().await?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1.0, 0.0, f64::NAN, 0.0],
    ///         [2.0, 0.0, f64::NAN, 2.0],
    ///         [0.0, 1.0, 2.0, 0.0]
    ///     ],
    /// );
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    ///
    /// Open the file and read data for one SNP (variant)
    /// at index position 2.
    /// ```
    /// # use ndarray as nd;
    /// # use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_object_path};
    /// # Runtime::new().unwrap().block_on(async {
    /// # let object_path = sample_bed_object_path("small.bed")?;
    ///
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let val = ReadOptions::builder().sid_index(2).f64().read_cloud(&mut bed_cloud).await?;
    ///
    /// assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    pub async fn new<I>(object_path: I) -> Result<Self, Box<BedErrorPlus>>
    where
        I: Into<ObjectPath<TObjectStore>>,
    {
        let object_path = object_path.into();
        BedCloud::builder(object_path).build().await
    }

    /// Number of individuals (samples)
    ///
    /// If this number is needed, it will be found
    /// by opening the .fam file and quickly counting the number
    /// of lines. Once found, the number will be remembered.
    /// The file read can be avoided by setting the
    /// number with [`BedCloudBuilder::iid_count`](struct.BedCloudBuilder.html#method.iid_count)
    /// or, for example, [`BedCloudBuilder::iid`](struct.BedCloudBuilder.html#method.iid).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_object_path};
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let iid_count = bed_cloud.iid_count().await?;
    ///
    /// assert!(iid_count == 3);
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    pub async fn iid_count(&mut self) -> Result<usize, Box<BedErrorPlus>> {
        if let Some(iid_count) = self.iid_count {
            Ok(iid_count)
        } else {
            let fam_object_path = self.fam_object_path()?;
            let iid_count = count_lines(&fam_object_path).await?;
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
    /// number with [`BedCloudBuilder::sid_count`](struct.BedCloudBuilder.html#method.sid_count)
    /// or, for example, [`BedCloudBuilder::sid`](struct.BedCloudBuilder.html#method.sid).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_object_path};
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let sid_count = bed_cloud.sid_count().await?;
    ///
    /// assert!(sid_count == 4);
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    pub async fn sid_count(&mut self) -> Result<usize, Box<BedErrorPlus>> {
        if let Some(sid_count) = self.sid_count {
            Ok(sid_count)
        } else {
            let bim_object_path = self.bim_object_path()?;
            let sid_count = count_lines(&bim_object_path).await?;
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
    /// number with [`BedCloudBuilder::iid_count`](struct.BedCloudBuilder.html#method.iid_count)
    /// and [`BedCloudBuilder::sid_count`](struct.BedCloudBuilder.html#method.sid_count).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let dim = bed_cloud.dim().await?;
    ///
    /// assert!(dim == (3,4));
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    // cmk call these at the same time?
    pub async fn dim(&mut self) -> Result<(usize, usize), Box<BedErrorPlus>> {
        Ok((self.iid_count().await?, self.sid_count().await?))
    }

    /// Family id of each of individual (sample)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .fam file. Once found, this ndarray
    /// and other information in the .fam file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedCloudBuilder::fid`](struct.BedCloudBuilder.html#method.fid).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let fid = bed_cloud.fid().await?;
    /// println!("{fid:?}"); // Outputs ndarray ["fid1", "fid1", "fid2"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    pub async fn fid(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_fam::<String>(self.metadata.fid.is_none(), MetadataFields::Fid, "fid")
            .await?;
        Ok(self.metadata.fid.as_ref().unwrap()) //unwrap always works because of lazy_fam
    }

    /// Individual id of each of individual (sample)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .fam file. Once found, this ndarray
    /// and other information in the .fam file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedCloudBuilder::iid`](struct.BedCloudBuilder.html#method.iid).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let iid = bed_cloud.iid().await?;    ///
    /// println!("{iid:?}"); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    pub async fn iid(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_fam::<String>(self.metadata.iid.is_none(), MetadataFields::Iid, "iid")
            .await?;
        Ok(self.metadata.iid.as_ref().unwrap()) //unwrap always works because of lazy_fam
    }

    /// Father id of each of individual (sample)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .fam file. Once found, this ndarray
    /// and other information in the .fam file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedCloudBuilder::father`](struct.BedCloudBuilder.html#method.father).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let father = bed_cloud.father().await?;
    /// println!("{father:?}"); // Outputs ndarray ["iid23", "iid23", "iid22"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    pub async fn father(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_fam::<String>(
            self.metadata.father.is_none(),
            MetadataFields::Father,
            "father",
        )
        .await?;
        Ok(self.metadata.father.as_ref().unwrap()) //unwrap always works because of lazy_fam
    }

    /// Mother id of each of individual (sample)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .fam file. Once found, this ndarray
    /// and other information in the .fam file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedCloudBuilder::mother`](struct.BedCloudBuilder.html#method.mother).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let mother = bed_cloud.mother().await?;
    /// println!("{mother:?}"); // Outputs ndarray ["iid34", "iid34", "iid33"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    pub async fn mother(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_fam::<String>(
            self.metadata.mother.is_none(),
            MetadataFields::Mother,
            "mother",
        )
        .await?;
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
    /// array with [`BedCloudBuilder::sex`](struct.BedCloudBuilder.html#method.sex).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let sex = bed_cloud.sex().await?;
    /// println!("{sex:?}"); // Outputs ndarray [1, 2, 0]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    pub async fn sex(&mut self) -> Result<&nd::Array1<i32>, Box<BedErrorPlus>> {
        self.unlazy_fam::<String>(self.metadata.sex.is_none(), MetadataFields::Sex, "sex")
            .await?;
        Ok(self.metadata.sex.as_ref().unwrap()) //unwrap always works because of lazy_fam
    }

    /// A phenotype for each individual (seldom used)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .fam file. Once found, this ndarray
    /// and other information in the .fam file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedCloudBuilder::pheno`](struct.BedCloudBuilder.html#method.pheno).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let pheno = bed_cloud.pheno().await?;
    /// println!("{pheno:?}"); // Outputs ndarray ["red", "red", "blue"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    pub async fn pheno(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_fam::<String>(
            self.metadata.pheno.is_none(),
            MetadataFields::Pheno,
            "pheno",
        )
        .await?;
        Ok(self.metadata.pheno.as_ref().unwrap()) //unwrap always works because of lazy_fam
    }

    /// Chromosome of each SNP (variant)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .bim file. Once found, this ndarray
    /// and other information in the .bim file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedCloudBuilder::chromosome`](struct.BedCloudBuilder.html#method.chromosome).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let chromosome = bed_cloud.chromosome().await?;
    /// println!("{chromosome:?}"); // Outputs ndarray ["1", "1", "5", "Y"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    pub async fn chromosome(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_bim::<String>(
            self.metadata.chromosome.is_none(),
            MetadataFields::Chromosome,
            "chromosome",
        )
        .await?;
        Ok(self.metadata.chromosome.as_ref().unwrap()) //unwrap always works because of lazy_bim
    }

    /// SNP id of each SNP (variant)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .bim file. Once found, this ndarray
    /// and other information in the .bim file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedCloudBuilder::sid`](struct.BedCloudBuilder.html#method.sid).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let sid = bed_cloud.sid().await?;
    /// println!("{sid:?}"); // Outputs ndarray "sid1", "sid2", "sid3", "sid4"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    pub async fn sid(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_bim::<String>(self.metadata.sid.is_none(), MetadataFields::Sid, "sid")
            .await?;
        Ok(self.metadata.sid.as_ref().unwrap()) //unwrap always works because of lazy_bim
    }

    /// Centimorgan position of each SNP (variant)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .bim file. Once found, this ndarray
    /// and other information in the .bim file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedCloudBuilder::cm_position`](struct.BedCloudBuilder.html#method.cm_position).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let cm_position = bed_cloud.cm_position().await?;
    /// println!("{cm_position:?}"); // Outputs ndarray [100.4, 2000.5, 4000.7, 7000.9]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    pub async fn cm_position(&mut self) -> Result<&nd::Array1<f32>, Box<BedErrorPlus>> {
        self.unlazy_bim::<String>(
            self.metadata.cm_position.is_none(),
            MetadataFields::CmPosition,
            "cm_position",
        )
        .await?;
        Ok(self.metadata.cm_position.as_ref().unwrap()) //unwrap always works because of lazy_bim
    }

    /// Base-pair position of each SNP (variant)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .bim file. Once found, this ndarray
    /// and other information in the .bim file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedCloudBuilder::bp_position`](struct.BedCloudBuilder.html#method.bp_position).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let bp_position = bed_cloud.bp_position().await?;
    /// println!("{bp_position:?}"); // Outputs ndarray [1, 100, 1000, 1004]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    pub async fn bp_position(&mut self) -> Result<&nd::Array1<i32>, Box<BedErrorPlus>> {
        self.unlazy_bim::<String>(
            self.metadata.bp_position.is_none(),
            MetadataFields::BpPosition,
            "bp_position",
        )
        .await?;
        Ok(self.metadata.bp_position.as_ref().unwrap()) //unwrap always works because of lazy_bim
    }

    /// First allele of each SNP (variant)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .bim file. Once found, this ndarray
    /// and other information in the .bim file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedCloudBuilder::allele_1`](struct.BedCloudBuilder.html#method.allele_1).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    /// # Runtime::new().unwrap().block_on(async {
    ///
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let allele_1 = bed_cloud.allele_1().await?;
    /// println!("{allele_1:?}"); // Outputs ndarray ["A", "T", "A", "T"]
    /// # let object_path = sample_bed_object_path("small.bed")?;
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    pub async fn allele_1(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_bim::<String>(
            self.metadata.allele_1.is_none(),
            MetadataFields::Allele1,
            "allele_1",
        )
        .await?;
        Ok(self.metadata.allele_1.as_ref().unwrap()) //unwrap always works because of lazy_bim
    }

    /// Second allele of each SNP (variant)
    ///
    /// If this ndarray is needed, it will be found
    /// by reading the .bim file. Once found, this ndarray
    /// and other information in the .bim file will be remembered.
    /// The file read can be avoided by setting the
    /// array with [`BedCloudBuilder::allele_2`](struct.BedCloudBuilder.html#method.allele_2).
    ///
    /// # Example:
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    /// # Runtime::new().unwrap().block_on(async {
    ///
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let allele_2 = bed_cloud.allele_2().await?;
    /// println!("{allele_2:?}"); // Outputs ndarray ["A", "C", "C", "G"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```        
    pub async fn allele_2(&mut self) -> Result<&nd::Array1<String>, Box<BedErrorPlus>> {
        self.unlazy_bim::<String>(
            self.metadata.allele_2.is_none(),
            MetadataFields::Allele2,
            "allele_2",
        )
        .await?;
        Ok(self.metadata.allele_2.as_ref().unwrap()) //unwrap always works because of lazy_bim
    }

    /// [`Metadata`](struct.Metadata.html) for this dataset, for example, the individual (sample) Ids.
    ///
    /// This returns a struct with 12 fields. Each field is a ndarray.
    /// The struct will always be new, but the 12 ndarrays will be
    /// shared with this [`BedCloud`](struct.BedCloud.html).
    ///
    /// If the needed, the metadata will be read from the .fam and/or .bim files.
    /// ```
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, sample_bed_object_path};
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let metadata = bed_cloud.metadata().await?;
    /// println!("{0:?}", metadata.iid()); // Outputs Some(["iid1", "iid2", "iid3"] ...)
    /// println!("{0:?}", metadata.sid()); // Outputs Some(["sid1", "sid2", "sid3", "sid4"] ...)
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    pub async fn metadata(&mut self) -> Result<Metadata, Box<BedErrorPlus>> {
        self.fam().await?;
        self.bim().await?;
        Ok(self.metadata.clone())
    }

    /// Return the `ObjectPath` of the .bed file.
    #[must_use]
    pub fn object_path(&self) -> &ObjectPath<TObjectStore> {
        &self.object_path
    }

    /// Return the cloud location of the .fam file.
    pub fn fam_object_path(&mut self) -> Result<ObjectPath<TObjectStore>, Box<BedErrorPlus>> {
        // We need to clone the object_path because self might mutate later
        if let Some(fam_object_path) = &self.fam_object_path {
            Ok(fam_object_path.clone())
        } else {
            let fam_object_path =
                to_metadata_path(&self.object_path, &self.fam_object_path, "fam")?;
            self.fam_object_path = Some(fam_object_path.clone());
            Ok(fam_object_path)
        }
    }

    /// Return the cloud location of the .bim file.
    pub fn bim_object_path(&mut self) -> Result<ObjectPath<TObjectStore>, Box<BedErrorPlus>> {
        // We need to clone the object_path because self might mutate later
        if let Some(bim_object_path) = &self.bim_object_path {
            Ok(bim_object_path.clone())
        } else {
            let bim_object_path =
                to_metadata_path(&self.object_path, &self.bim_object_path, "bim")?;
            self.bim_object_path = Some(bim_object_path.clone());
            Ok(bim_object_path)
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let val = bed_cloud.read::<f64>().await?;
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
    /// let val = bed_cloud.read::<i8>().await?;
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1, 0, -127, 0],
    ///         [2, 0, -127, 2],
    ///         [0, 1, 2, 0]
    ///     ],
    /// );
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```    
    pub async fn read<TVal: BedVal>(&mut self) -> Result<nd::Array2<TVal>, Box<BedErrorPlus>> {
        let read_options = ReadOptions::<TVal>::builder().build()?;
        self.read_with_options(&read_options).await
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// // Read the SNPs indexed by 2.
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let read_options = ReadOptions::builder().sid_index(2).build()?;
    /// let mut val = nd::Array2::<f64>::default((3, 1));
    /// bed_cloud.read_and_fill_with_options(&mut val.view_mut(), &read_options).await?;
    ///
    /// assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```  
    #[allow(clippy::similar_names)]
    pub async fn read_and_fill_with_options<TVal: BedVal>(
        &mut self,
        val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.,
        read_options: &ReadOptions<TVal>,
    ) -> Result<(), Box<BedErrorPlus>> {
        // must do these one-at-a-time because they mutate self to cache the results
        let iid_count = self.iid_count().await?;
        let sid_count = self.sid_count().await?;

        let max_concurrent_requests =
            compute_max_concurrent_requests(read_options.max_concurrent_requests)?;

        let max_chunk_size = compute_max_chunk_size(read_options.max_chunk_size)?;

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
            &self.object_path,
            iid_count,
            sid_count,
            read_options.is_a1_counted,
            iid_index,
            sid_index,
            read_options.missing_value,
            max_concurrent_requests,
            max_chunk_size,
            &mut val.view_mut(),
        )
        .await
    }

    /// cmk doc
    // have read_and_fill_with_options call this
    pub async fn read_and_fill_with_options_no_mut<TVal: BedVal>(
        &self,
        iid_count: usize,
        sid_count: usize,
        val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.,
        read_options: &ReadOptions<TVal>,
    ) -> Result<(), Box<BedErrorPlus>> {
        // // must do these one-at-a-time because they mutate self to cache the results
        // let iid_count = self.iid_count().await?;
        // let sid_count = self.sid_count().await?;

        let max_concurrent_requests =
            compute_max_concurrent_requests(read_options.max_concurrent_requests)?;

        let max_chunk_size = compute_max_chunk_size(read_options.max_chunk_size)?;

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
            &self.object_path,
            iid_count,
            sid_count,
            read_options.is_a1_counted,
            iid_index,
            sid_index,
            read_options.missing_value,
            max_concurrent_requests,
            max_chunk_size,
            &mut val.view_mut(),
        )
        .await
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let mut val = nd::Array2::<i8>::default(bed_cloud.dim().await?);
    /// bed_cloud.read_and_fill(&mut val.view_mut()).await?;
    ///
    /// assert_eq_nan(
    ///     &val,
    ///     &nd::array![
    ///         [1, 0, -127, 0],
    ///         [2, 0, -127, 2],
    ///         [0, 1, 2, 0]
    ///     ],
    /// );
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    pub async fn read_and_fill<TVal: BedVal>(
        &mut self,
        val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.,
    ) -> Result<(), Box<BedErrorPlus>> {
        let read_options = ReadOptions::<TVal>::builder().build()?;
        self.read_and_fill_with_options(val, &read_options).await
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_object_path};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// // Read the SNPs indexed by 2.
    /// let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(object_path).await?;
    /// let read_options = ReadOptions::builder().sid_index(2).f64().build()?;
    /// let val = bed_cloud.read_with_options(&read_options).await?;
    ///
    /// assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```  
    pub async fn read_with_options<TVal: BedVal>(
        &mut self,
        read_options: &ReadOptions<TVal>,
    ) -> Result<nd::Array2<TVal>, Box<BedErrorPlus>> {
        let iid_count_in = self.iid_count().await?;
        let sid_count_in = self.sid_count().await?;
        let iid_count_out = read_options.iid_index.len(iid_count_in)?;
        let sid_count_out = read_options.sid_index.len(sid_count_in)?;
        let shape = ShapeBuilder::set_f((iid_count_out, sid_count_out), read_options.is_f);
        let mut val = nd::Array2::<TVal>::default(shape);

        self.read_and_fill_with_options(&mut val.view_mut(), read_options)
            .await?;

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
    /// ```ignore // cmk
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, WriteOptions};
    ///
    /// let output_folder = temp_testdir::TempDir::default();
    /// let output_file = output_folder.join("small.bed");
    ///
    /// let val = nd::array![[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]];
    /// BedCloud::write(&val, &output_file)?;
    ///
    /// // If we then read the new file and list the chromosome property,
    /// // it is an array of zeros, the default chromosome value.
    /// let mut bed_cloud2 = BedCloud::new(&output_file)?;
    /// println!("{:?}", bed_cloud2.chromosome().await?); // Outputs ndarray ["0", "0", "0", "0"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    // cmk need to do 'write'
    // pub fn write<S: nd::Data<Elem = TVal>, TVal: BedVal>(
    //     val: &nd::ArrayBase<S, nd::Ix2>,
    //     path: &Path,
    // ) -> Result<(), Box<BedErrorPlus>> {
    //     WriteOptions::builder(path).write(val)
    // }

    /// Given an 2D array of genotype data and a [`WriteOptions`](struct.WriteOptionsBuilder.html), write to a .bed file.
    ///
    /// > Also see [`WriteOptionsBuilder::write`](struct.WriteOptionsBuilder.html#method.write), which creates
    /// > a [`WriteOptions`](struct.WriteOptionsBuilder.html) and writes to file in one step.
    ///
    /// # Example
    /// ```ignore // cmk
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, WriteOptions};
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
    /// BedCloud::write_with_options(&val, &write_options)?;
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    // cmk need to do 'write_with_options'
    // pub fn write_with_options<S, TVal>(
    //     val: &nd::ArrayBase<S, nd::Ix2>,
    //     write_options: &WriteOptions<TVal>,
    // ) -> Result<(), Box<BedErrorPlus>>
    // where
    //     S: nd::Data<Elem = TVal>,
    //     TVal: BedVal,
    // {
    //     let (iid_count, sid_count) = val.dim();
    //     if iid_count != write_options.iid_count() {
    //         return Err(BedError::InconsistentCount(
    //             "iid".to_string(),
    //             write_options.iid_count(),
    //             iid_count,
    //         )
    //         .into());
    //     }
    //     if sid_count != write_options.sid_count() {
    //         return Err(BedError::InconsistentCount(
    //             "sid".to_string(),
    //             write_options.sid_count(),
    //             sid_count,
    //         )
    //         .into());
    //     }

    //     let num_threads = compute_num_threads(write_options.num_threads)?;
    //     write_val(
    //         &write_options.path,
    //         val,
    //         write_options.is_a1_counted,
    //         write_options.missing_value,
    //         num_threads,
    //     )?;

    //     if !write_options.skip_fam() {
    //         if let Err(e) = write_options.metadata.write_fam(write_options.fam_object_path()) {
    //             // Clean up the file
    //             let _ = fs::remove_file(&write_options.fam_object_path);
    //             return Err(e);
    //         }
    //     }

    //     if !write_options.skip_bim() {
    //         if let Err(e) = write_options.metadata.write_bim(write_options.bim_object_path()) {
    //             // Clean up the file
    //             let _ = fs::remove_file(&write_options.bim_object_path);
    //             return Err(e);
    //         }
    //     }

    //     Ok(())
    // }

    async fn unlazy_fam<T: FromStringArray<T>>(
        &mut self,
        is_none: bool,
        field_index: MetadataFields,
        name: &str,
    ) -> Result<(), Box<BedErrorPlus>> {
        if self.skip_set.contains(&field_index) {
            return Err(BedError::CannotUseSkippedMetadata(name.to_string()).into());
        }
        if is_none {
            self.fam().await?;
        }
        Ok(())
    }

    async fn unlazy_bim<T: FromStringArray<T>>(
        &mut self,
        is_none: bool,
        field_index: MetadataFields,
        name: &str,
    ) -> Result<(), Box<BedErrorPlus>> {
        if self.skip_set.contains(&field_index) {
            return Err(BedError::CannotUseSkippedMetadata(name.to_string()).into());
        }
        if is_none {
            self.bim().await?;
        }
        Ok(())
    }

    async fn fam(&mut self) -> Result<(), Box<BedErrorPlus>> {
        let fam_object_path = self.fam_object_path()?;

        let (metadata, count) = self
            .metadata
            .read_fam_cloud(&fam_object_path, &self.skip_set)
            .await?;
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

    async fn bim(&mut self) -> Result<(), Box<BedErrorPlus>> {
        let bim_object_path = self.bim_object_path()?;

        let (metadata, count) = self
            .metadata
            .read_bim_cloud(&bim_object_path, &self.skip_set)
            .await?;
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

/// Returns the cloud location of a sample .bed file.
///
/// Behind the scenes, the "cloud location" will actually be local.
/// If necessary, the file will be downloaded.
/// The .fam and .bim files will also be downloaded, if they are not already present.
/// SHA256 hashes are used to verify that the files are correct.
/// The files will be in a directory determined by environment variable `BED_READER_DATA_DIR`.
/// If that environment variable is not set, a cache folder, appropriate to the OS, will be used.
#[anyinput]
pub fn sample_bed_object_path(
    bed_path: AnyPath,
) -> Result<ObjectPath<LocalFileSystem>, Box<BedErrorPlus>> {
    let mut path_list: Vec<PathBuf> = Vec::new();
    for ext in &["bed", "bim", "fam"] {
        let file_path = bed_path.with_extension(ext);
        path_list.push(file_path);
    }

    let mut vec = sample_object_paths(path_list)?;
    debug_assert!(vec.len() == 3);
    Ok(vec.swap_remove(0))
}

/// Returns the cloud location of a sample file.
///
/// Behind the scenes, the "cloud location" will actually be local.
/// If necessary, the file will be downloaded.
/// A SHA256 hash is used to verify that the file is correct.
/// The file will be in a directory determined by environment variable `BED_READER_DATA_DIR`.
/// If that environment variable is not set, a cache folder, appropriate to the OS, will be used.
#[anyinput]
pub fn sample_object_path(path: AnyPath) -> Result<ObjectPath<LocalFileSystem>, Box<BedErrorPlus>> {
    let object_store = Arc::new(LocalFileSystem::new());

    let file_path = STATIC_FETCH_DATA
        .fetch_file(path)
        .map_err(BedErrorPlus::from)?;
    let store_path = StorePath::from_filesystem_path(file_path).map_err(BedErrorPlus::from)?;
    let object_path = ObjectPath::new(object_store, store_path);
    Ok(object_path)
}

/// Returns the cloud locations of a list of files.
///
/// Behind the scenes, the "cloud location" will actually be local.
/// If necessary, the file will be downloaded.
/// SHA256 hashes are used to verify that the files are correct.
/// The files will be in a directory determined by environment variable `BED_READER_DATA_DIR`.
/// If that environment variable is not set, a cache folder, appropriate to the OS, will be used.
#[anyinput]
pub fn sample_object_paths(
    path_list: AnyIter<AnyPath>,
) -> Result<Vec<ObjectPath<LocalFileSystem>>, Box<BedErrorPlus>> {
    let object_store = Arc::new(LocalFileSystem::new());

    let file_paths = STATIC_FETCH_DATA
        .fetch_files(path_list)
        .map_err(BedErrorPlus::from)?;
    file_paths
        .iter()
        .map(|file_path| {
            let path = StorePath::from_filesystem_path(file_path).map_err(BedErrorPlus::from)?;
            Ok(ObjectPath::new(object_store.clone(), path))
        })
        .collect()
}

#[derive(Debug)]
/// The location of a file in the cloud.
///
/// The location is made up of of two parts, an `Arc`-wrapped [`ObjectStore`] and a [`StorePath`].
/// The [`ObjectStore`] is a file server, for example, AWS S3, Azure, the local file system, etc.
/// The [`StorePath`] is the path to the file on the file server.
///
///
/// [`ObjectStore`]: object_store::ObjectStore
/// [`StorePath`]: object_store::path::Path
///
/// # Examples
///
/// You can create an `ObjectPath` from a tuple -- with or without any references.
/// ```
/// use std::sync::Arc;
/// use object_store::{local::LocalFileSystem, path::Path as StorePath};
/// use bed_reader::{ObjectPath, BedErrorPlus, sample_bed_file};
///
/// # Runtime::new().unwrap().block_on(async {
/// let object_store = Arc::new(LocalFileSystem::new()); // Arc-wrapped ObjectStore
/// let file_path = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?; // regular Rust PathBuf
/// let store_path = StorePath::from_filesystem_path(&file_path).map_err(BedErrorPlus::from)?; // StorePath
///
/// let object_path0  = ObjectPath::<_>::from(&(&object_store, &store_path)); // ObjectPath from references
/// let object_path1: ObjectPath<_> = (object_store, store_path).into(); // ObjectPath from owned values
///
/// assert_eq!(object_path0.size().await?, object_path1.size().await?);
/// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
/// # use {tokio::runtime::Runtime};
/// ```
///
/// You can skip `Arc`-wrapping, but then the [`ObjectStore`] must be owned.
/// ```
/// # use std::sync::Arc;
/// # use object_store::{local::LocalFileSystem, path::Path as StorePath};
/// # use bed_reader::{ObjectPath, BedErrorPlus, sample_bed_file};
/// # Runtime::new().unwrap().block_on(async {
/// let object_store = LocalFileSystem::new(); // ObjectStore
/// let file_path = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?; // regular Rust PathBuf
/// let store_path = StorePath::from_filesystem_path(&file_path).map_err(BedErrorPlus::from)?; // StorePath
///
/// let object_path: ObjectPath<_> = (object_store, &store_path).into(); // ObjectPath from owned object_store
/// assert_eq!(object_path.size().await?, 303);
/// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
/// # use {tokio::runtime::Runtime};
/// ```
///
/// Alternatively, you can use [`ObjectPath::new`](struct.ObjectPath.html#method.new), but both parts must be owned
/// and the [`ObjectStore`] must be `Arc`-wrapped.
/// ```
/// # use std::sync::Arc;
/// # use object_store::{local::LocalFileSystem, path::Path as StorePath};
/// # use bed_reader::{ObjectPath, BedErrorPlus, sample_bed_file};
/// # Runtime::new().unwrap().block_on(async {
/// let object_store = Arc::new(LocalFileSystem::new()); // Arc-wrapped ObjectStore
/// let file_path = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?; // regular Rust PathBuf
/// let store_path = StorePath::from_filesystem_path(&file_path).map_err(BedErrorPlus::from)?; // StorePath
///
/// let object_path: ObjectPath<_> = ObjectPath::new(object_store, store_path); // ObjectPath from owned values
/// assert_eq!(object_path.size().await?, 303);
/// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
/// # use {tokio::runtime::Runtime};
/// ```
pub struct ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    /// cmk doc
    pub object_store: Arc<TObjectStore>,
    /// cmk doc
    pub path: StorePath,
}

impl<TObjectStore> Clone for ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn clone(&self) -> Self {
        ObjectPath {
            object_store: self.object_store.clone(),
            path: self.path.clone(),
        }
    }
}

impl<TObjectStore> ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    /// Create a new [`ObjectPath`] from an `Arc`-wrapped [`ObjectStore`] and a [`StorePath`].
    ///
    /// Both parts must be owned, but see [`ObjectPath`] for examples of creating from a tuple with references.
    ///
    /// # Example
    /// ```
    /// use std::sync::Arc;
    /// use object_store::{local::LocalFileSystem, path::Path as StorePath};
    /// use bed_reader::{ObjectPath, BedErrorPlus, sample_bed_file};
    /// # Runtime::new().unwrap().block_on(async {
    /// let object_store = Arc::new(LocalFileSystem::new()); // Arc-wrapped ObjectStore
    /// let file_path = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?; // regular Rust PathBuf
    /// let store_path = StorePath::from_filesystem_path(&file_path).map_err(BedErrorPlus::from)?; // StorePath
    ///
    /// let object_path: ObjectPath<_> = ObjectPath::new(object_store, store_path); // ObjectPath from owned values
    /// assert_eq!(object_path.size().await?, 303);
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime};
    /// ```

    pub fn new(object_store: Arc<TObjectStore>, path: StorePath) -> Self {
        ObjectPath { object_store, path }
    }

    /// Return the size of a file stored in the cloud.
    ///
    /// # Example
    /// ```
    /// use bed_reader::{sample_bed_object_path};
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let mut object_path = sample_bed_object_path("plink_sim_10s_100v_10pmiss.bed")?;
    /// assert_eq!(object_path.size().await?, 303);
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    pub async fn size(&self) -> Result<usize, Box<BedErrorPlus>> {
        let get_result = self.get().await?;
        let object_meta = &get_result.meta; // cmk good idea?
        Ok(object_meta.size)
    }

    /// Return the bytes that are stored at the specified location in the given byte ranges
    pub async fn get_ranges(
        &self,
        ranges: &[core::ops::Range<usize>],
    ) -> Result<Vec<Bytes>, Box<BedErrorPlus>> {
        self.object_store
            .get_ranges(&self.path, ranges)
            .await
            .map_err(|e| Box::new(BedErrorPlus::from(e)))
    }

    /// Perform a get request with options
    pub async fn get_opts(&self, get_options: GetOptions) -> Result<GetResult, Box<BedErrorPlus>> {
        self.object_store
            .get_opts(&self.path, get_options)
            .await
            .map_err(|e| Box::new(BedErrorPlus::from(e)))
    }

    /// Return the bytes that are stored at the specified location.
    pub async fn get(&self) -> Result<GetResult, Box<BedErrorPlus>> {
        self.object_store
            .get(&self.path)
            .await
            .map_err(|e| Box::new(BedErrorPlus::from(e)))
    }

    /// Updates the [`ObjectPath`] to have the given extension.
    ///
    /// It removes the current extension, if any.
    /// It appends the given extension, if any.
    ///
    /// # Example
    /// ```
    /// use bed_reader::{sample_bed_object_path};
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let mut object_path = sample_bed_object_path("plink_sim_10s_100v_10pmiss.bed")?;
    /// assert_eq!(object_path.size().await?, 303);
    /// object_path.set_extension("fam")?;
    /// assert_eq!(object_path.size().await?, 130);
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    pub fn set_extension(&mut self, extension: &str) -> Result<(), Box<BedErrorPlus>> {
        let mut path_str = self.path.to_string();

        // Find the last dot in the object path
        if let Some(dot_index) = path_str.rfind('.') {
            // Remove the current extension
            path_str.truncate(dot_index);
        }

        if !extension.is_empty() {
            // Append the new extension
            path_str.push('.');
            path_str.push_str(extension);
        }

        // Parse the string back to StorePath
        self.path = StorePath::parse(&path_str).map_err(BedErrorPlus::from)?;
        Ok(())
    }
}

// Implementing From trait for ObjectPath to allow tuple conversions.
impl<TObjectStore> From<(Arc<TObjectStore>, StorePath)> for ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn from(tuple: (Arc<TObjectStore>, StorePath)) -> Self {
        ObjectPath {
            object_store: tuple.0,
            path: tuple.1,
        }
    }
}
impl<TObjectStore> From<(&Arc<TObjectStore>, &StorePath)> for ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn from(tuple: (&Arc<TObjectStore>, &StorePath)) -> Self {
        ObjectPath {
            object_store: tuple.0.clone(),
            path: tuple.1.clone(),
        }
    }
}
impl<TObjectStore> From<(Arc<TObjectStore>, &StorePath)> for ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn from(tuple: (Arc<TObjectStore>, &StorePath)) -> Self {
        ObjectPath {
            object_store: tuple.0,
            path: tuple.1.clone(),
        }
    }
}
impl<TObjectStore> From<(&Arc<TObjectStore>, StorePath)> for ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn from(tuple: (&Arc<TObjectStore>, StorePath)) -> Self {
        ObjectPath {
            object_store: tuple.0.clone(),
            path: tuple.1,
        }
    }
}

impl<TObjectStore> From<&(Arc<TObjectStore>, StorePath)> for ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn from(tuple: &(Arc<TObjectStore>, StorePath)) -> Self {
        ObjectPath {
            object_store: tuple.0.clone(),
            path: tuple.1.clone(),
        }
    }
}
impl<TObjectStore> From<&(&Arc<TObjectStore>, &StorePath)> for ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn from(tuple: &(&Arc<TObjectStore>, &StorePath)) -> Self {
        ObjectPath {
            object_store: tuple.0.clone(),
            path: tuple.1.clone(),
        }
    }
}
impl<TObjectStore> From<&(Arc<TObjectStore>, &StorePath)> for ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn from(tuple: &(Arc<TObjectStore>, &StorePath)) -> Self {
        ObjectPath {
            object_store: tuple.0.clone(),
            path: tuple.1.clone(),
        }
    }
}
impl<TObjectStore> From<&(&Arc<TObjectStore>, StorePath)> for ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn from(tuple: &(&Arc<TObjectStore>, StorePath)) -> Self {
        ObjectPath {
            object_store: tuple.0.clone(),
            path: tuple.1.clone(),
        }
    }
}

// Implementing From trait for ObjectPath to allow tuple conversions.
impl<TObjectStore> From<(TObjectStore, StorePath)> for ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn from(tuple: (TObjectStore, StorePath)) -> Self {
        ObjectPath {
            object_store: Arc::new(tuple.0),
            path: tuple.1,
        }
    }
}
impl<TObjectStore> From<(TObjectStore, &StorePath)> for ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn from(tuple: (TObjectStore, &StorePath)) -> Self {
        ObjectPath {
            object_store: Arc::new(tuple.0),
            path: tuple.1.clone(),
        }
    }
}
impl<TObjectStore> From<&ObjectPath<TObjectStore>> for ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn from(ref_thing: &ObjectPath<TObjectStore>) -> Self {
        ref_thing.clone()
    }
}

impl<TObjectStore> fmt::Display for ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ObjectPath: {:?}", self.path)
    }
}

fn to_metadata_path<TObjectStore>(
    bed_object_path: &ObjectPath<TObjectStore>,
    metadata_object_path: &Option<ObjectPath<TObjectStore>>,
    extension: &str,
) -> Result<ObjectPath<TObjectStore>, Box<BedErrorPlus>>
where
    TObjectStore: ObjectStore,
{
    if let Some(metadata_object_path) = metadata_object_path {
        Ok(metadata_object_path.clone())
    } else {
        let mut meta_object_path = bed_object_path.clone();
        meta_object_path.set_extension(extension)?;
        Ok(meta_object_path)
    }
}

async fn count_lines<TObjectStore, I>(object_path: I) -> Result<usize, Box<BedErrorPlus>>
where
    TObjectStore: ObjectStore,
    I: Into<ObjectPath<TObjectStore>>,
{
    let stream = object_path.into().get().await?.into_stream();

    let new_line_stream = newline_delimited_stream(stream);

    let newline_count = AtomicUsize::new(0);
    new_line_stream
        .try_for_each(|bytes| {
            let count = bytecount::count(&bytes, b'\n');
            newline_count.fetch_add(count, Ordering::SeqCst);
            async { Ok(()) } // Return Ok(()) for each successful iteration
        })
        .await
        .map_err(BedErrorPlus::from)?; // Convert the error and propagate it if present

    Ok(newline_count.load(Ordering::SeqCst))
}
