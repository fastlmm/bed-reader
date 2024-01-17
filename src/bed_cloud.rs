#[cfg(not(target_pointer_width = "64"))]
compile_error!(
    "This code requires a 64-bit target architecture because the cloud-library assumes it."
);

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
use object_store::{GetOptions, GetRange, GetResult};
use object_store::{ObjectMeta, ObjectStore};
use std::cmp::max;
use std::collections::HashSet;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use url::Url;

use crate::{
    check_and_precompute_iid_index, compute_max_chunk_size, compute_max_concurrent_requests,
    set_up_two_bits_to_value, try_div_4, BedError, BedErrorPlus, BedVal, FromStringArray, Hold,
    Metadata, ReadOptions, BED_FILE_MAGIC1, BED_FILE_MAGIC2, STATIC_FETCH_DATA,
};
use crate::{MetadataFields, CB_HEADER_U64};

/// Represents a PLINK .bed file in the cloud that is open for reading genotype data and metadata.
///
/// Construct with [`BedCloud::new`](struct.BedCloud.html#method.new), [`BedCloud::builder`](struct.BedCloud.html#method.builder),
/// [`BedCloud::from_object_path`](struct.BedCloud.html#method.from_object_path), or
/// [`BedCloud::builder_from_object_path`](struct.BedCloud.html#method.builder_from_object_path).
///
/// > For reading local files, see [`Bed`](struct.Bed.html).
///
/// # Example
///
/// Open a file for reading. Then, read the individual (sample) ids
/// and all the genotype data.
/// ```
/// use ndarray as nd;
/// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
/// use bed_reader::assert_eq_nan;
///
/// # Runtime::new().unwrap().block_on(async {
/// let url = sample_bed_url("small.bed")?;
/// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
                None => Err(BedError::UninitializedField("object_path"))?,
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
                None => Err(BedError::UninitializedField("metadata"))?,
            },
            skip_set: match self.skip_set {
                Some(ref value) => Clone::clone(value),
                None => Err(BedError::UninitializedField("skip_set"))?,
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
        Ok((in_sid_i_signed - lower_sid_count) as u64)
    } else {
        Err(Box::new(BedErrorPlus::BedError(BedError::SidIndexTooBig(
            in_sid_i_signed,
        ))))
    }
}

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
        // LATER: Consider doing this in parallel as in the non-cloud version.
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
fn check_file_length<TObjectStore>(
    in_iid_count: usize,
    in_sid_count: usize,
    size: usize,
    object_path: &ObjectPath<TObjectStore>,
) -> Result<(usize, u64), Box<BedErrorPlus>>
where
    TObjectStore: ObjectStore,
{
    let (in_iid_count_div4, in_iid_count_div4_u64) =
        try_div_4(in_iid_count, in_sid_count, CB_HEADER_U64)?;
    let file_len = size as u64;
    let file_len2 = in_iid_count_div4_u64 * (in_sid_count as u64) + CB_HEADER_U64;
    if file_len != file_len2 {
        Err(BedError::IllFormed(object_path.to_string()))?;
    }
    Ok((in_iid_count_div4, in_iid_count_div4_u64))
}

#[inline]
#[allow(clippy::too_many_arguments)]
#[allow(clippy::similar_names)]
async fn read_no_alloc<TVal: BedVal, TObjectStore>(
    object_path: &ObjectPath<TObjectStore>,
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
{
    let (size, bytes) = open_and_check(object_path).await?;

    match bytes[2] {
        0 => {
            // We swap 'iid' and 'sid' and then reverse the axes.
            let mut val_t = val.view_mut().reversed_axes();

            internal_read_no_alloc(
                object_path,
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
                object_path,
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
        _ => Err(BedError::BadMode(object_path.to_string()))?,
    };
    Ok(())
}

async fn open_and_check<TObjectStore>(
    object_path: &ObjectPath<TObjectStore>,
) -> Result<(usize, Bytes), Box<BedErrorPlus>>
where
    TObjectStore: ObjectStore,
{
    let object_store = object_path.object_store.clone();
    let path: &StorePath = &object_path.path;
    let object_meta: ObjectMeta = object_store.head(path).await?;
    let size: usize = object_meta.size;

    let get_options = GetOptions {
        range: Some(GetRange::Bounded(0..CB_HEADER_U64 as usize)),
        ..Default::default()
    };
    let object_store = object_path.object_store.clone();
    let path: &StorePath = &object_path.path;
    let get_result = object_store.get_opts(path, get_options).await?;
    let bytes = get_result.bytes().await?;

    if (BED_FILE_MAGIC1 != bytes[0]) || (BED_FILE_MAGIC2 != bytes[1]) {
        Err(BedError::IllFormed(object_path.to_string()))?;
    }
    Ok((size, bytes))
}

impl BedCloudBuilder<Box<dyn ObjectStore>> {
    fn new<S, I, K, V>(url: S, options: I) -> Result<Self, Box<BedErrorPlus>>
    where
        S: AsRef<str>,
        I: IntoIterator<Item = (K, V)>,
        K: AsRef<str>,
        V: Into<String>,
    {
        let object_path = ObjectPath::from_url(url, options)?;
        Ok(BedCloudBuilder::from(object_path))
    }

    /// Set the cloud location of the .fam file. Specify the file with a URL string.
    ///
    /// If not set, the .fam file will be assumed
    /// to have the same location as the .bed file, but with the extension .fam.
    ///
    /// > See [`BedCloudBuider::fam_object_path`](struct.BedCloudBuilder.html#method.fam_object_path) to specify the file with an [`ObjectPath`](struct.ObjectPath.html)
    /// > instead of a URL string.
    ///
    /// # Example:
    /// Read .bed, .fam, and .bim files with non-standard names.
    /// ```
    /// use bed_reader::{BedCloud, ReadOptions, sample_urls, EMPTY_OPTIONS};
    /// # Runtime::new().unwrap().block_on(async {
    /// let deb_maf_mib = sample_urls(["small.deb", "small.maf", "small.mib"])?;
    /// let mut bed_cloud = BedCloud::builder(&deb_maf_mib[0], EMPTY_OPTIONS)?
    ///    .fam_from_url(&deb_maf_mib[1], EMPTY_OPTIONS)?
    ///    .bim_from_url(&deb_maf_mib[2], EMPTY_OPTIONS)?
    ///    .build().await?;
    /// println!("{:?}", bed_cloud.iid().await?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// println!("{:?}", bed_cloud.sid().await?); // Outputs ndarray ["sid1", "sid2", "sid3", "sid4"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    pub fn fam_from_url<S, I, K, V>(mut self, url: S, options: I) -> Result<Self, Box<BedErrorPlus>>
    where
        S: AsRef<str>,
        I: IntoIterator<Item = (K, V)>,
        K: AsRef<str>,
        V: Into<String>,
    {
        let object_path = ObjectPath::from_url(url, options)?;
        self.fam_object_path = Some(Some(object_path));
        Ok(self)
    }

    /// Set the cloud location of the .bim file. Specify the file with a URL string.
    ///
    /// If not set, the .bim file will be assumed
    /// to have the same location as the .bed file, but with the extension .bim.
    ///
    /// > See [`BedCloudBuider::fam_object_path`](struct.BedCloudBuilder.html#method.bim_object_path) to specify the file with an [`ObjectPath`](struct.ObjectPath.html)
    /// > instead of a URL string.
    ///
    /// # Example:
    /// Read .bed, .fam, and .bim files with non-standard names.
    /// ```
    /// use bed_reader::{BedCloud, ReadOptions, sample_urls, EMPTY_OPTIONS};
    /// # Runtime::new().unwrap().block_on(async {
    /// let deb_maf_mib = sample_urls(["small.deb", "small.maf", "small.mib"])?;
    /// let mut bed_cloud = BedCloud::builder(&deb_maf_mib[0], EMPTY_OPTIONS)?
    ///    .fam_from_url(&deb_maf_mib[1], EMPTY_OPTIONS)?
    ///    .bim_from_url(&deb_maf_mib[2], EMPTY_OPTIONS)?
    ///    .build().await?;
    /// println!("{:?}", bed_cloud.iid().await?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// println!("{:?}", bed_cloud.sid().await?); // Outputs ndarray ["sid1", "sid2", "sid3", "sid4"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    pub fn bim_from_url<S, I, K, V>(mut self, url: S, options: I) -> Result<Self, Box<BedErrorPlus>>
    where
        S: AsRef<str>,
        I: IntoIterator<Item = (K, V)>,
        K: AsRef<str>,
        V: Into<String>,
    {
        let object_path = ObjectPath::from_url(url, options)?;
        self.bim_object_path = Some(Some(object_path));
        Ok(self)
    }
}

impl<TObjectStore> From<&ObjectPath<TObjectStore>> for BedCloudBuilder<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn from(object_path: &ObjectPath<TObjectStore>) -> Self {
        Self {
            object_path: Some(object_path.clone()), // Cloned here.
            fam_object_path: None,
            bim_object_path: None,

            is_checked_early: None,
            iid_count: None,
            sid_count: None,

            metadata: Some(Metadata::new()),
            skip_set: Some(HashSet::new()),
        }
    }
}

impl<TObjectStore> From<ObjectPath<TObjectStore>> for BedCloudBuilder<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    fn from(object_path: ObjectPath<TObjectStore>) -> Self {
        Self {
            object_path: Some(object_path), // Cloned here.
            fam_object_path: None,
            bim_object_path: None,

            is_checked_early: None,
            iid_count: None,
            sid_count: None,

            metadata: Some(Metadata::new()),
            skip_set: Some(HashSet::new()),
        }
    }
}

impl<TObjectStore> BedCloudBuilder<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    /// Create a [`BedCloud`](struct.BedCloud.html) from the builder.
    ///
    /// > See [`BedCloud::builder`](struct.BedCloud.html#method.builder) for more details and examples.
    pub async fn build(&self) -> Result<BedCloud<TObjectStore>, Box<BedErrorPlus>> {
        let mut bed_cloud = self.build_no_file_check()?;

        // Unwrap is allowed becaue we can't construct BedCloudBuilder without object_path
        if bed_cloud.is_checked_early {
            let object_path = self.object_path.as_ref().unwrap().clone();
            open_and_check(&object_path).await?;
        }

        (bed_cloud.iid_count, bed_cloud.sid_count) = bed_cloud
            .metadata
            .check_counts(bed_cloud.iid_count, bed_cloud.sid_count)?;

        Ok(bed_cloud)
    }

    /// Override the family id (fid) values found in the .fam file.
    ///
    /// By default, if fid values are needed and haven't already been found,
    /// they will be read from the .fam file.
    /// Providing them here avoids that file read and provides a way to give different values.
    #[anyinput]
    #[must_use]
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
    /// use bed_reader::{BedCloud, assert_eq_nan, sample_bed_url, EMPTY_OPTIONS};
    /// let url = sample_bed_url("small.bed")?;
    /// use bed_reader::ReadOptions;
    ///
    /// let mut bed_cloud = BedCloud::builder(url, EMPTY_OPTIONS)?
    ///    .iid(["sample1", "sample2", "sample3"])
    ///    .build().await?;
    /// println!("{:?}", bed_cloud.iid().await?); // Outputs ndarray ["sample1", "sample2", "sample3"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    #[anyinput]
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    /// use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_url, EMPTY_OPTIONS};
    /// let url = sample_bed_url("small.bed")?;
    ///
    /// let mut bed_cloud = BedCloud::builder(url, EMPTY_OPTIONS)?
    ///    .sid(["SNP1", "SNP2", "SNP3", "SNP4"])
    ///    .build().await?;
    /// println!("{:?}", bed_cloud.sid().await?); // Outputs ndarray ["SNP1", "SNP2", "SNP3", "SNP4"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    #[anyinput]
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    /// # use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_url, EMPTY_OPTIONS};
    /// # let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::builder(&url, EMPTY_OPTIONS)?.skip_early_check().build().await?;
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
    #[must_use]
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_urls, EMPTY_OPTIONS};
    /// # Runtime::new().unwrap().block_on(async {
    /// let deb_maf_mib = sample_urls(["small.deb", "small.maf", "small.mib"])?;
    /// let mut bed_cloud = BedCloud::builder(&deb_maf_mib[0], EMPTY_OPTIONS)?
    ///    .fam_from_url(&deb_maf_mib[1], EMPTY_OPTIONS)?
    ///    .bim_from_url(&deb_maf_mib[2], EMPTY_OPTIONS)?
    ///    .build().await?;
    /// println!("{:?}", bed_cloud.iid().await?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// println!("{:?}", bed_cloud.sid().await?); // Outputs ndarray ["sid1", "sid2", "sid3", "sid4"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    #[must_use]
    pub fn fam_object_path(mut self, object_path: &ObjectPath<TObjectStore>) -> Self {
        self.fam_object_path = Some(Some(object_path.clone()));
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
    /// let mut bed_cloud = BedCloud::builder_from_object_path(&deb_maf_mib[0])
    ///    .fam_object_path(&deb_maf_mib[1])
    ///    .bim_object_path(&deb_maf_mib[2])
    ///    .build().await?;
    /// println!("{:?}", bed_cloud.iid().await?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// println!("{:?}", bed_cloud.sid().await?); // Outputs ndarray ["sid1", "sid2", "sid3", "sid4"]
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    #[must_use]
    pub fn bim_object_path(mut self, object_path: &ObjectPath<TObjectStore>) -> Self {
        let object_path = object_path.clone();
        self.bim_object_path = Some(Some(object_path));
        self
    }

    /// Don't read the fid information from the .fam file.
    ///
    /// By default, when the .fam is read, the fid (the family id) is recorded.
    /// This stops that recording. This is useful if the fid is not needed.
    /// Asking for the fid after skipping it results in an error.    
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
    /// use bed_reader::{BedCloud, Metadata, sample_bed_url, EMPTY_OPTIONS};
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let metadata = Metadata::builder()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build()?;
    /// let mut bed_cloud = BedCloud::builder(url, EMPTY_OPTIONS)?
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
    #[must_use]
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

impl BedCloud<Box<dyn ObjectStore>> {
    #[allow(clippy::doc_link_with_quotes)]
    /// Attempts to open a PLINK .bed file in the cloud for reading. The file is specified with a URL string.
    ///
    /// See ["Cloud URLs and `ObjectPath` Examples"](supplemental_document_cloud_urls/index.html) for details specifying a file.
    ///
    /// You may give [cloud options](supplemental_document_options/index.html#cloud-options) but not
    /// [`BedCloud` options](supplemental_document_options/index.html#bedbedcloud-options) or
    /// [`ReadOptions`](supplemental_document_options/index.html#readoptions).
    /// See ["Options, Options, Options"](supplemental_document_options/index.html) for details.
    ///
    /// > Also see [`BedCloud::builder`](struct.BedCloud.html#method.builder), which does support
    /// > `BedCloud` options.
    /// > Alternatively, you can use [`BedCloud::builder_from_object_path`](struct.BedCloud.html#method.builder_from_object_path)
    /// > to specify the cloud file via an [`ObjectPath`](struct.ObjectPath.html). For reading local files,
    /// > see [`Bed`](struct.Bed.html).
    ///
    /// # Errors
    /// URL parsing may return an error.
    /// Also, by default, this method will return an error if the file is missing or its header
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
    /// use bed_reader::{BedCloud, assert_eq_nan, sample_bed_url, EMPTY_OPTIONS};
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// println!("{url:?}"); // For example, "file:///C:/Users/carlk/AppData/Local/bed_reader/bed_reader/Cache/small.bed"
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// # use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_url, EMPTY_OPTIONS};
    /// # Runtime::new().unwrap().block_on(async {
    /// # let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
    /// let val = ReadOptions::builder().sid_index(2).f64().read_cloud(&mut bed_cloud).await?;
    ///
    /// assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    pub async fn new<S, I, K, V>(url: S, cloud_options: I) -> Result<Self, Box<BedErrorPlus>>
    where
        S: AsRef<str>,
        I: IntoIterator<Item = (K, V)>,
        K: AsRef<str>,
        V: Into<String>,
    {
        let object_path = ObjectPath::from_url(url, cloud_options)?;
        let bed_cloud = BedCloud::from_object_path(&object_path).await?;
        Ok(bed_cloud)
    }

    #[allow(clippy::doc_link_with_quotes)]
    /// Attempts to open a PLINK .bed file in the cloud for reading. The file is specified with a URL string.
    /// Supports both [cloud options](supplemental_document_options/index.html#cloud-options) and
    /// [`BedCloud` options](supplemental_document_options/index.html#bedbedcloud-options).
    ///
    /// See ["Cloud URLs and `ObjectPath` Examples"](supplemental_document_cloud_urls/index.html) for details of specifying a file.
    /// See ["Options, Options, Options"](supplemental_document_options/index.html) for an overview of options types.
    ///
    /// > Also see [`BedCloud::new`](struct.BedCloud.html#method.url),
    /// > which does not support `BedCloud` options.
    /// > Alternatively, you can use [`BedCloud::builder_from_object_path`](struct.BedCloud.html#method.builder_from_object_path)
    /// > to specify the cloud file via an [`ObjectPath`](struct.ObjectPath.html). For reading local files,
    /// > see [`Bed`](struct.Bed.html).
    ///
    /// The `BedCloud` options, [listed here](struct.BedCloudBuilder.html#implementations), can:
    ///  * set the cloud location of the .fam and/or .bim file
    ///  * override some metadata, for example, replace the individual ids.
    ///  * set the number of individuals (samples) or SNPs (variants)
    ///  * control checking the validity of the .bed file's header
    ///  * skip reading selected metadata
    ///
    /// # Errors
    /// URL parsing may return an error.
    /// Also, by default, this method will return an error if the file is missing or its header
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
    /// use bed_reader::{BedCloud, assert_eq_nan, sample_bed_url, EMPTY_OPTIONS};
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// println!("{url:?}"); // For example, "file:///C:/Users/carlk/AppData/Local/bed_reader/bed_reader/Cache/small.bed"
    /// let mut bed_cloud = BedCloud::builder(url, EMPTY_OPTIONS)?.build().await?;
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
    /// # use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_url, EMPTY_OPTIONS};
    /// # let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::builder(url, EMPTY_OPTIONS)?
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
    /// # use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_url, EMPTY_OPTIONS};
    /// # let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::builder(&url, EMPTY_OPTIONS)?
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
    /// # use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_url, EMPTY_OPTIONS};
    /// # let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::builder(url, EMPTY_OPTIONS)?
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
    pub fn builder<S, I, K, V>(
        url: S,
        options: I,
    ) -> Result<BedCloudBuilder<Box<dyn ObjectStore>>, Box<BedErrorPlus>>
    where
        S: AsRef<str>,
        I: IntoIterator<Item = (K, V)>,
        K: AsRef<str>,
        V: Into<String>,
    {
        BedCloudBuilder::new(url, options)
    }
}

impl<TObjectStore> BedCloud<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    /// Attempts to open a PLINK .bed file in the cloud for reading. Specify the file with an [`ObjectPath`].
    /// Supports `BedCloud` options.
    ///
    /// See ["Cloud URLs and `ObjectPath` Examples"](supplemental_document_cloud_urls/index.html) for details of specifying a file.
    ///
    /// > Also see [`BedCloud::from_object_path`](struct.BedCloud.html#method.from_object_path)
    /// > which does not support `BedCloud` options.
    /// > Alternatively, you can use [`BedCloud::builder`](struct.BedCloud.html#method.builder)
    /// > to specify the cloud file via a URL string. For reading local files,
    /// > see [`Bed`](struct.Bed.html).
    ///
    /// The `BedCloud` options, [listed here](struct.BedCloudBuilder.html#implementations), can:
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
    /// let mut bed_cloud = BedCloud::builder_from_object_path(&object_path).build().await?;
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
    /// # use bed_reader::{BedCloud, assert_eq_nan, sample_bed_object_path};
    /// # let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
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
    /// # use bed_reader::{BedCloud, assert_eq_nan, sample_bed_object_path};
    /// # let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
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
    /// # use bed_reader::{BedCloud, assert_eq_nan, sample_bed_object_path};
    /// # let object_path = sample_bed_object_path("small.bed")?;
    /// let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
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
    #[must_use]
    pub fn builder_from_object_path(
        object_path: &ObjectPath<TObjectStore>,
    ) -> BedCloudBuilder<TObjectStore> {
        BedCloudBuilder::from(object_path)
    }

    /// Attempts to open a PLINK .bed file in the cloud for reading. Specify the file with an [`ObjectPath`].
    ///
    /// See ["Cloud URLs and `ObjectPath` Examples"](supplemental_document_cloud_urls/index.html) for details specifying a file.
    ///
    /// You may not give
    /// [`BedCloud` options](supplemental_document_options/index.html#bedbedcloud-options).
    /// See [`BedCloud::builder_from_object_path`](struct.BedCloud.html#method.builder_from_object_path), which does support
    /// `BedCloud` options.
    ///
    /// > Also see, [`BedCloud::builder`](struct.BedCloud.html#method.builder)
    /// > to specify the cloud file via a URL string. For reading local files,
    /// > see [`Bed`](struct.Bed.html).
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
    /// let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
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
    /// let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    /// let val = ReadOptions::builder().sid_index(2).f64().read_cloud(&mut bed_cloud).await?;
    ///
    /// assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    /// ```
    pub async fn from_object_path(
        object_path: &ObjectPath<TObjectStore>,
    ) -> Result<Self, Box<BedErrorPlus>> {
        BedCloudBuilder::from(object_path).build().await
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
    /// use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_url, EMPTY_OPTIONS};
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_url, EMPTY_OPTIONS};
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
    /// let dim = bed_cloud.dim().await?;
    ///
    /// assert!(dim == (3,4));
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
    // LATER: Could these be called at the same time, async?
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    /// # Runtime::new().unwrap().block_on(async {
    ///
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
    /// let allele_1 = bed_cloud.allele_1().await?;
    /// println!("{allele_1:?}"); // Outputs ndarray ["A", "T", "A", "T"]
    /// # let url = sample_bed_url("small.bed")?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    /// # Runtime::new().unwrap().block_on(async {
    ///
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, sample_bed_url, EMPTY_OPTIONS};
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    pub fn object_path(&self) -> ObjectPath<TObjectStore> {
        self.object_path.clone()
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// // Read the SNPs indexed by 2.
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
            Err(BedError::InvalidShape(
                iid_index.len(),
                sid_index.len(),
                dim.0,
                dim.1,
            ))?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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
    /// use bed_reader::{BedCloud, ReadOptions, sample_bed_url, EMPTY_OPTIONS};
    /// use bed_reader::assert_eq_nan;
    ///
    /// # Runtime::new().unwrap().block_on(async {
    /// // Read the SNPs indexed by 2.
    /// let url = sample_bed_url("small.bed")?;
    /// let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
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

    // LATER: Support writing to a BedCloud

    async fn unlazy_fam<T: FromStringArray<T>>(
        &mut self,
        is_none: bool,
        field_index: MetadataFields,
        name: &str,
    ) -> Result<(), Box<BedErrorPlus>> {
        if self.skip_set.contains(&field_index) {
            Err(BedError::CannotUseSkippedMetadata(name.into()))?;
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
            Err(BedError::CannotUseSkippedMetadata(name.into()))?;
        }
        if is_none {
            self.bim().await?;
        }
        Ok(())
    }

    async fn fam(&mut self) -> Result<(), Box<BedErrorPlus>> {
        let fam_object_path = self.fam_object_path()?.clone();

        let (metadata, count) = self
            .metadata
            .read_fam_cloud(&fam_object_path, &self.skip_set)
            .await?;
        self.metadata = metadata;

        match self.iid_count {
            Some(iid_count) => {
                if iid_count != count {
                    Err(BedError::InconsistentCount("iid".into(), iid_count, count))?;
                }
            }
            None => {
                self.iid_count = Some(count);
            }
        }
        Ok(())
    }

    async fn bim(&mut self) -> Result<(), Box<BedErrorPlus>> {
        let bim_object_path = self.bim_object_path()?.clone();

        let (metadata, count) = self
            .metadata
            .read_bim_cloud(&bim_object_path, &self.skip_set)
            .await?;
        self.metadata = metadata;

        match self.sid_count {
            Some(sid_count) => {
                if sid_count != count {
                    Err(BedError::InconsistentCount("sid".into(), sid_count, count))?;
                }
            }
            None => {
                self.sid_count = Some(count);
            }
        }
        Ok(())
    }
}

/// Returns the cloud locations of a .bed file as an [`ObjectPath`](struct.ObjectPath.html).
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

/// Returns the cloud locations of a file as an [`ObjectPath`](struct.ObjectPath.html).
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
        .map_err(|e| BedError::SampleFetch(e.to_string()))?;
    let store_path = StorePath::from_filesystem_path(file_path)?;
    let object_path = ObjectPath::new(object_store, store_path);
    Ok(object_path)
}

/// Returns the cloud locations of a list of files as [`ObjectPath`](struct.ObjectPath.html)s.
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
    let arc_object_store = Arc::new(LocalFileSystem::new());

    let file_paths = STATIC_FETCH_DATA
        .fetch_files(path_list)
        .map_err(|e| BedError::SampleFetch(e.to_string()))?;
    file_paths
        .iter()
        .map(|file_path| {
            let path = StorePath::from_filesystem_path(file_path)?;
            Ok(ObjectPath::new(arc_object_store.clone(), path))
        })
        .collect()
}

/// Returns the cloud location of a sample .bed file as a URL string.
///
/// Behind the scenes, the "cloud location" will actually be local.
/// If necessary, the file will be downloaded.
/// The .fam and .bim files will also be downloaded, if they are not already present.
/// SHA256 hashes are used to verify that the files are correct.
/// The files will be in a directory determined by environment variable `BED_READER_DATA_DIR`.
/// If that environment variable is not set, a cache folder, appropriate to the OS, will be used.
#[anyinput]
pub fn sample_bed_url(bed_path: AnyPath) -> Result<String, Box<BedErrorPlus>> {
    let mut path_list: Vec<PathBuf> = Vec::new();
    for ext in &["bed", "bim", "fam"] {
        let file_path = bed_path.with_extension(ext);
        path_list.push(file_path);
    }

    let mut vec = sample_urls(path_list)?;
    Ok(vec.swap_remove(0))
}

/// Returns a local file's path as a URL string, taking care of any needed encoding.
///
/// ```
/// use bed_reader::{path_to_url_string, sample_bed_file};
/// let path = sample_bed_file("small.bed")?;
/// let url: String = path_to_url_string(path)?;
/// assert!(url.starts_with("file:///") && url.ends_with("/small.bed"));
/// # use bed_reader::BedErrorPlus;
/// # Ok::<(), Box<BedErrorPlus>>(())
/// ```
#[anyinput]
pub fn path_to_url_string(path: AnyPath) -> Result<String, Box<BedErrorPlus>> {
    let url = Url::from_file_path(path)
        .map_err(|_e| BedError::CannotCreateUrlFromFilePath(path.to_string_lossy().to_string()))?
        .to_string();
    Ok(url)
}

/// Returns the cloud location of a sample file as a URL string.
///
/// Behind the scenes, the "cloud location" will actually be local.
/// If necessary, the file will be downloaded.
/// A SHA256 hash is used to verify that the file is correct.
/// The file will be in a directory determined by environment variable `BED_READER_DATA_DIR`.
/// If that environment variable is not set, a cache folder, appropriate to the OS, will be used.
#[anyinput]
pub fn sample_url(path: AnyPath) -> Result<String, Box<BedErrorPlus>> {
    let file_path = STATIC_FETCH_DATA
        .fetch_file(path)
        .map_err(|e| BedError::SampleFetch(e.to_string()))?;
    let url = path_to_url_string(file_path)?;
    Ok(url)
}

/// Returns the cloud locations of a list of files as URL strings.
///
/// Behind the scenes, the "cloud location" will actually be local.
/// If necessary, the file will be downloaded.
/// SHA256 hashes are used to verify that the files are correct.
/// The files will be in a directory determined by environment variable `BED_READER_DATA_DIR`.
/// If that environment variable is not set, a cache folder, appropriate to the OS, will be used.
#[anyinput]
pub fn sample_urls(path_list: AnyIter<AnyPath>) -> Result<Vec<String>, Box<BedErrorPlus>> {
    let file_paths = STATIC_FETCH_DATA
        .fetch_files(path_list)
        .map_err(|e| BedError::SampleFetch(e.to_string()))?;
    file_paths
        .iter()
        .map(|file_path| {
            let url = path_to_url_string(file_path)?;
            Ok(url)
        })
        .collect()
}

#[derive(Debug)]
/// The location of a file in the cloud.
///
/// The location is made up of of two parts, an `Arc`-wrapped [`ObjectStore`](https://docs.rs/object_store/latest/object_store/trait.ObjectStore.html) and an [`object_store::path::Path as StorePath`](https://docs.rs/object_store/latest/object_store/path/struct.Path.html).
/// The [`ObjectStore`](https://docs.rs/object_store/latest/object_store/trait.ObjectStore.html) is a cloud service, for example, AWS S3, Azure, the local file system, etc.
/// The `StorePath` is the path to the file on the cloud service.
///
/// See ["Cloud URLs and `ObjectPath` Examples"](supplemental_document_cloud_urls/index.html) for details specifying a file.
///
/// An `ObjectPath` can be efficiently cloned because the `ObjectStore` is `Arc`-wrapped.
///
/// # Examples
///
/// ```
/// use std::sync::Arc;
/// use object_store::{local::LocalFileSystem, path::Path as StorePath};
/// use bed_reader::{ObjectPath, BedErrorPlus, sample_bed_file};
///
/// # Runtime::new().unwrap().block_on(async {
/// let arc_object_store = Arc::new(LocalFileSystem::new()); // Arc-wrapped ObjectStore
/// let file_path = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?; // regular Rust PathBuf
/// let store_path = StorePath::from_filesystem_path(&file_path)?; // StorePath
/// let object_path = ObjectPath::new(arc_object_store, store_path); // ObjectPath
/// assert_eq!(object_path.size().await?, 303);
/// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
/// # use {tokio::runtime::Runtime};
/// ```
pub struct ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    /// An `Arc`-wrapped [`ObjectStore`](https://docs.rs/object_store/latest/object_store/trait.ObjectStore.html) cloud service, for example, AWS S3, Azure, the local file system, etc.
    pub object_store: Arc<TObjectStore>,
    /// A [`object_store::path::Path as StorePath`](https://docs.rs/object_store/latest/object_store/path/struct.Path.html) that points to a file on the [`ObjectStore`](https://docs.rs/object_store/latest/object_store/trait.ObjectStore.html)
    /// that gives the path to the file on the cloud service.
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

/// An empty set of [cloud options](supplemental_document_options/index.html#cloud-options)
///
/// See ["Cloud URLs and `ObjectPath` Examples"](supplemental_document_cloud_urls/index.html) for examples.
pub const EMPTY_OPTIONS: [(&str, String); 0] = [];

impl ObjectPath<Box<dyn ObjectStore>> {
    /// Create a new [`ObjectPath`] from a URL string and [cloud options](supplemental_document_options/index.html#cloud-options).
    ///
    /// See ["Cloud URLs and `ObjectPath` Examples"](supplemental_document_cloud_urls/index.html) for details specifying a file.
    ///
    /// # Example
    /// ```
    /// use std::sync::Arc;
    /// use object_store::{local::LocalFileSystem, path::Path as StorePath};
    /// use bed_reader::{ObjectPath, BedErrorPlus, sample_bed_url, EMPTY_OPTIONS};
    /// # Runtime::new().unwrap().block_on(async {
    /// let url: String = sample_bed_url("plink_sim_10s_100v_10pmiss.bed")?;
    /// let object_path: ObjectPath<_> = ObjectPath::from_url(url, EMPTY_OPTIONS)?;
    /// assert_eq!(object_path.size().await?, 303);
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime};
    /// ```
    pub fn from_url<I, K, V, S>(
        location: S,
        options: I,
    ) -> Result<ObjectPath<Box<dyn ObjectStore>>, Box<BedErrorPlus>>
    where
        I: IntoIterator<Item = (K, V)>,
        K: AsRef<str>,
        V: Into<String>,
        S: AsRef<str>,
    {
        let location = location.as_ref();
        let url = Url::parse(location)
            .map_err(|e| BedError::CannotParseUrl(location.to_string(), e.to_string()))?;

        let (object_store, store_path): (Box<dyn ObjectStore>, StorePath) =
            object_store::parse_url_opts(&url, options)?;
        let object_path = ObjectPath::new(Arc::new(object_store), store_path);
        Ok(object_path)
    }
}

impl<TObjectStore> ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    /// Create a new [`ObjectPath`] from an `Arc`-wrapped [`ObjectStore`](https://docs.rs/object_store/latest/object_store/trait.ObjectStore.html) and an [`object_store::path::Path as StorePath`](https://docs.rs/object_store/latest/object_store/path/struct.Path.html).
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
    /// let store_path = StorePath::from_filesystem_path(&file_path)?; // StorePath
    ///
    /// let object_path: ObjectPath<_> = ObjectPath::new(object_store, store_path); // ObjectPath from owned values
    /// assert_eq!(object_path.size().await?, 303);
    /// # Ok::<(), Box<BedErrorPlus>>(())}).unwrap();
    /// # use {tokio::runtime::Runtime};
    /// ```
    pub fn new(arc_object_store: Arc<TObjectStore>, path: StorePath) -> Self {
        ObjectPath {
            object_store: arc_object_store,
            path,
        }
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
        // LATER: See if https://github.com/apache/arrow-rs/issues/5272 if fixed in
        // a way so that only one read is needed.
        let object_meta = &get_result.meta;
        Ok(object_meta.size)
    }

    /// Return the bytes that are stored at the specified location in the given byte ranges
    pub async fn get_ranges(
        &self,
        ranges: &[core::ops::Range<usize>],
    ) -> Result<Vec<Bytes>, Box<BedErrorPlus>> {
        Ok(self.object_store.get_ranges(&self.path, ranges).await?)
    }

    /// Perform a get request with options
    pub async fn get_opts(&self, get_options: GetOptions) -> Result<GetResult, Box<BedErrorPlus>> {
        Ok(self.object_store.get_opts(&self.path, get_options).await?)
    }

    /// Return the bytes that are stored at the specified location.
    pub async fn get(&self) -> Result<GetResult, Box<BedErrorPlus>> {
        Ok(self.object_store.get(&self.path).await?)
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
        self.path = StorePath::parse(&path_str)?;
        Ok(())
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

async fn count_lines<TObjectStore>(
    object_path: &ObjectPath<TObjectStore>,
) -> Result<usize, Box<BedErrorPlus>>
where
    TObjectStore: ObjectStore,
{
    let stream = object_path.get().await?.into_stream();

    let new_line_stream = newline_delimited_stream(stream);

    let newline_count = AtomicUsize::new(0);
    new_line_stream
        .try_for_each(|bytes| {
            let count = bytecount::count(&bytes, b'\n');
            newline_count.fetch_add(count, Ordering::SeqCst);
            async { Ok(()) } // Return Ok(()) for each successful iteration
        })
        .await?; // Convert the error and propagate it if present

    Ok(newline_count.load(Ordering::SeqCst))
}
