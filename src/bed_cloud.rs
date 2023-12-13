use anyinput::anyinput;
use derive_builder::Builder;
use nd::ShapeBuilder;
use ndarray as nd;
use std::collections::HashSet;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use tokio::sync::Semaphore;
use tokio::task::{self, JoinHandle};

// cmk really need futures_util and bytes...???
use bytes::Bytes;
use futures_util::TryStreamExt;
use itertools::Itertools;
use object_store::delimited::newline_delimited_stream;
use object_store::path::Path as StorePath;
use object_store::ObjectStore;
use object_store::{GetOptions, ObjectMeta};

use crate::{
    check_and_precompute_iid_index, compute_max_chunk_size, compute_max_concurrent_requests,
    set_up_two_bits_to_value, try_div_4, BedError, BedErrorPlus, BedVal, Hold, Metadata,
    ReadOptions, BED_FILE_MAGIC1, BED_FILE_MAGIC2,
};
use crate::{MetadataFields, CB_HEADER_U64};

/// cmk doc
#[derive(Clone, Debug, Builder)]
#[builder(build_fn(private, name = "build_no_file_check", error = "BedErrorPlus"))]
pub struct BedCloud<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    #[builder(setter(custom))]
    object_store: Arc<TObjectStore>,

    #[builder(setter(custom))]
    path: StorePath,

    // cmk do we want to cache the ObjectMeta?
    #[builder(setter(custom))]
    #[builder(default = "None")]
    fam_path: Option<StorePath>,

    // // cmk needed? combine with path?
    // #[builder(setter(custom))]
    // #[builder(default = "None")]
    // fam_object_meta: Option<ObjectMeta>,
    #[builder(setter(custom))]
    #[builder(default = "None")]
    bim_path: Option<StorePath>,

    // #[builder(setter(custom))]
    // #[builder(default = "None")]
    // bim_object_meta: Option<ObjectMeta>,
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

impl<TObjectStore> BedCloud<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    /// cmk doc
    #[anyinput]
    pub fn builder(
        object_store: &Arc<TObjectStore>,
        path: &StorePath,
    ) -> BedCloudBuilder<TObjectStore> {
        BedCloudBuilder::new(object_store, path)
    }

    /// cmk doc
    // cmk #[anyinput]
    pub async fn new(
        object_store: &Arc<TObjectStore>,
        path: &StorePath,
    ) -> Result<Self, Box<BedErrorPlus>> {
        // cmk do this?? let path = path.into();
        BedCloud::builder(object_store, path).build().await
    }

    // #[anyinput]
    async fn count_lines(&self, path: &StorePath) -> Result<usize, Box<BedErrorPlus>> {
        let stream = self
            .object_store
            .clone()
            .get(path)
            .await
            .map_err(BedErrorPlus::from)?
            .into_stream();

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

    // cmk #[anyinput]
    fn to_metadata_path(
        bed_path: &StorePath,
        metadata_path: &Option<StorePath>,
        extension: &str,
    ) -> Result<StorePath, Box<BedErrorPlus>> {
        if let Some(metadata_path) = metadata_path {
            Ok(metadata_path.to_owned())
        } else {
            change_extension(bed_path, extension)
        }
    }

    /// Return the path of the .fam file.
    pub fn fam_path(&mut self) -> Result<StorePath, Box<BedErrorPlus>> {
        // We need to clone the path because self might mutate later
        if let Some(path) = &self.fam_path {
            Ok(path.clone())
        } else {
            let path =
                BedCloud::<TObjectStore>::to_metadata_path(&self.path, &self.fam_path, "fam")?;
            self.fam_path = Some(path.clone());
            Ok(path)
        }
    }

    /// Return the path of the .bim file.
    pub fn bim_path(&mut self) -> Result<StorePath, Box<BedErrorPlus>> {
        // We need to clone the path because self might mutate later
        if let Some(path) = &self.bim_path {
            Ok(path.clone())
        } else {
            let path =
                BedCloud::<TObjectStore>::to_metadata_path(&self.path, &self.bim_path, "bim")?;
            self.bim_path = Some(path.clone());
            Ok(path)
        }
    }

    /// cmk doc
    pub async fn iid_count(&mut self) -> Result<usize, Box<BedErrorPlus>> {
        if let Some(iid_count) = self.iid_count {
            Ok(iid_count)
        } else {
            let fam_path = self.fam_path()?;
            let iid_count = self.count_lines(&fam_path).await?;
            self.iid_count = Some(iid_count);
            Ok(iid_count)
        }
    }

    /// cmk doc
    pub async fn sid_count(&mut self) -> Result<usize, Box<BedErrorPlus>> {
        if let Some(sid_count) = self.sid_count {
            Ok(sid_count)
        } else {
            let bim_path = self.bim_path()?;
            let sid_count = self.count_lines(&bim_path).await?;
            self.sid_count = Some(sid_count);
            Ok(sid_count)
        }
    }

    /// cmk doc
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
            &self.object_store, // cmk rename "object_store" ??
            &self.path,
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

    /// cmk doc
    pub async fn read<TVal: BedVal>(&mut self) -> Result<nd::Array2<TVal>, Box<BedErrorPlus>> {
        let read_options = ReadOptions::<TVal>::builder().build()?;
        self.read_with_options(&read_options).await
    }
}

fn change_extension(
    store_path: &StorePath,
    new_extension: &str,
) -> Result<StorePath, Box<BedErrorPlus>> {
    let mut path_str = store_path.to_string();

    // Find the last dot in the path
    if let Some(dot_index) = path_str.rfind('.') {
        // Remove the current extension
        path_str.truncate(dot_index);
    }

    // Append the new extension
    path_str.push('.');
    path_str.push_str(new_extension);

    // Parse the string back to StorePath
    StorePath::parse(&path_str).map_err(|e| BedErrorPlus::from(e).into())
}

// cmk #[anyinput]
// cmk needed?
fn store_path_to_string(store_path: StorePath) -> String {
    StorePath::to_string(&store_path)
}

#[allow(clippy::too_many_arguments)]
// cmk #[anyinput]
async fn internal_read_no_alloc<TVal: BedVal, TStore: ObjectStore>(
    object_store: &Arc<TStore>,
    path: &StorePath,
    object_meta: &ObjectMeta,
    in_iid_count: usize,
    in_sid_count: usize,
    is_a1_counted: bool,
    iid_index: &[isize],
    sid_index: &[isize],
    missing_value: TVal,
    max_concurrent_requests: usize,
    max_chunk_size: usize,
    out_val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.
) -> Result<(), Box<BedErrorPlus>> {
    // Check the file length

    let (in_iid_count_div4, in_iid_count_div4_u64) =
        try_div_4(in_iid_count, in_sid_count, CB_HEADER_U64)?;
    // "as" and math is safe because of early checks
    let file_len = object_meta.size as u64;
    let file_len2 = in_iid_count_div4_u64 * (in_sid_count as u64) + CB_HEADER_U64;
    if file_len != file_len2 {
        return Err(Box::new(
            BedError::IllFormed(store_path_to_string(path.clone())).into(),
        ));
    }

    // Check and precompute for each iid_index

    let (i_div_4_array, i_mod_4_times_2_array) =
        check_and_precompute_iid_index(in_iid_count, iid_index)?;

    let chunk_size = std::cmp::max(1, max_chunk_size / in_iid_count_div4);

    // Check and compute work for each sid_index

    let from_two_bits_to_value = set_up_two_bits_to_value(is_a1_counted, missing_value);
    let lower_sid_count = -(in_sid_count as isize);
    let upper_sid_count: isize = (in_sid_count as isize) - 1;

    let semaphore = Arc::new(Semaphore::new(max_concurrent_requests));

    for chunk in &sid_index.iter().chunks(chunk_size) {
        // ======== prepare ranges and cols ================ cmk refactor
        let mut ranges = Vec::with_capacity(chunk_size);
        // let mut cols = Vec::with_capacity(chunk_size);
        let mut in_sid_i_vec = Vec::with_capacity(chunk_size);
        for in_sid_i_signed in chunk {
            // cmk similar code elsewhere
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

            let pos: u64 = in_sid_i * in_iid_count_div4_u64 + CB_HEADER_U64; // "as" and math is safe because of early checks
            let range = pos as usize..(pos + in_iid_count_div4_u64) as usize;
            ranges.push(range);
            // cols.push(col);
            in_sid_i_vec.push(in_sid_i);
        }

        let permit = semaphore.clone().acquire_owned().await.unwrap(); // cmk unwrap
        let path_clone = path.clone(); // cmk fast enough?
        let object_store_clone = object_store.clone(); // This is Arc fast
        let handle: JoinHandle<Result<_, BedErrorPlus>> = task::spawn(async move {
            // cmk somehow we must only compile is size(usize) is 64 bits.
            // cmk we should turn sid_index into a slice of ranges.
            let vec_bytes = object_store_clone
                .get_ranges(&path_clone, &ranges)
                .await
                .map_err(BedErrorPlus::from)?; // cmk unwrap
            Ok((vec_bytes, in_sid_i_vec))
        });

        match handle.await {
            Ok(Ok((vec_bytes, in_sid_i_vec))) => {
                drop(permit);
                for (bytes, in_sid_i) in vec_bytes.into_iter().zip(in_sid_i_vec.into_iter()) {
                    let mut col = out_val.column_mut(in_sid_i as usize);
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
            Ok(Err(e)) => return Err(e.into()),
            Err(e) => return Err(BedErrorPlus::from(e).into()),
        }
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
// cmk #[anyinput]
async fn read_no_alloc<TVal: BedVal, TStore: ObjectStore>(
    object_store: &Arc<TStore>,
    path: &StorePath,
    iid_count: usize,
    sid_count: usize,
    is_a1_counted: bool,
    iid_index: &[isize],
    sid_index: &[isize],
    missing_value: TVal,
    max_concurrent_requests: usize,
    max_chunk_size: usize,

    val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.
) -> Result<(), Box<BedErrorPlus>> {
    let (object_meta, bytes) = open_and_check(object_store, path).await?;

    match bytes[2] {
        0 => {
            // We swap 'iid' and 'sid' and then reverse the axes.
            let mut val_t = val.view_mut().reversed_axes();

            internal_read_no_alloc(
                object_store,
                path,
                &object_meta,
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
                object_store,
                path,
                &object_meta,
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
            .await?
        }
        _ => {
            return Err(Box::new(
                BedError::BadMode(store_path_to_string(path.clone())).into(),
            ))
        }
    };
    Ok(())
}

// cmk #[anyinput]
async fn open_and_check<TStore>(
    object_store: &Arc<TStore>,
    path: &StorePath,
) -> Result<(ObjectMeta, Bytes), Box<BedErrorPlus>>
where
    TStore: ObjectStore,
{
    let get_options = GetOptions {
        range: Some(0..CB_HEADER_U64 as usize),
        ..Default::default()
    };
    let get_result = object_store
        .get_opts(path, get_options)
        .await
        .map_err(BedErrorPlus::from)?;

    let object_meta = get_result.meta.clone(); // cmk good idea?
    let bytes = get_result.bytes().await.map_err(BedErrorPlus::from)?;

    if (BED_FILE_MAGIC1 != bytes[0]) || (BED_FILE_MAGIC2 != bytes[1]) {
        return Err(Box::new(
            BedError::IllFormed(store_path_to_string(path.clone())).into(),
        ));
    }
    Ok((object_meta, bytes))
}

impl<TObjectStore> BedCloudBuilder<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    // #[anyinput]
    fn new(object_store: &Arc<TObjectStore>, path: &StorePath) -> Self {
        Self {
            object_store: Some(object_store.clone()),
            path: Some(path.clone()),
            fam_path: None,
            bim_path: None,

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

        // cmk is this unwrap OK?
        if bed_cloud.is_checked_early {
            open_and_check(&self.object_store.unwrap(), &bed_cloud.path).await?;
        }

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
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, assert_eq_nan, sample_bed_file};
    /// let file_name = sample_bed_file("small.bed_cloud")?;
    /// use bed_reader::ReadOptions;
    ///
    /// let mut bed_cloud = BedCloud::builder(file_name)
    ///    .iid(["sample1", "sample2", "sample3"])
    ///    .build()?;
    /// println!("{:?}", bed_cloud.iid()?); // Outputs ndarray ["sample1", "sample2", "sample3"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    #[anyinput]
    pub fn iid(mut self, iid: AnyIter<AnyString>) -> Self {
        // Unwrap will always work because BedCloudBuilder starting with some metadata
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
    /// use ndarray as nd;
    /// use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_file};
    /// let file_name = sample_bed_file("small.bed_cloud")?;
    ///
    /// let mut bed_cloud = BedCloud::builder(file_name)
    ///    .sid(["SNP1", "SNP2", "SNP3", "SNP4"])
    ///    .build()?;
    /// println!("{:?}", bed_cloud.sid()?); // Outputs ndarray ["SNP1", "SNP2", "SNP3", "SNP4"]
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

    /// Don't check the header of the .bed_cloud file until and unless the file is actually read.
    ///
    /// By default, when a [`BedCloud`](struct.BedCloud.html) struct is created, the .bed_cloud
    /// file header is checked. This stops that early check.
    pub fn skip_early_check(mut self) -> Self {
        self.is_checked_early = Some(false);
        self
    }

    /// Set the path to the .fam file.
    ///
    /// If not set, the .fam file will be assumed
    /// to have the same name as the .bed_cloud file, but with the extension .fam.
    ///
    /// # Example:
    /// Read .bed_cloud, .fam, and .bim files with non-standard names.
    /// ```
    /// use bed_reader::{BedCloud, ReadOptions, sample_files};
    /// let deb_maf_mib = sample_files(["small.deb", "small.maf", "small.mib"])?;
    /// let mut bed_cloud = BedCloud::builder(&deb_maf_mib[0])
    ///    .fam_path(&deb_maf_mib[1])
    ///    .bim_path(&deb_maf_mib[2])
    ///    .build()?;
    /// println!("{:?}", bed_cloud.iid()?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// println!("{:?}", bed_cloud.sid()?); // Outputs ndarray ["sid1", "sid2", "sid3", "sid4"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    // cmk #[anyinput]
    pub fn fam_path(mut self, path: StorePath) -> Self {
        self.fam_path = Some(Some(path.to_owned()));
        self
    }

    /// Set the path to the .bim file.
    ///
    /// If not set, the .bim file will be assumed
    /// to have the same name as the .bed_cloud file, but with the extension .bim.
    ///
    /// # Example:
    /// Read .bed_cloud, .fam, and .bim files with non-standard names.
    /// ```
    /// use bed_reader::{BedCloud, ReadOptions, sample_files};
    /// let deb_maf_mib = sample_files(["small.deb", "small.maf", "small.mib"])?;
    /// let mut bed_cloud = BedCloud::builder(&deb_maf_mib[0])
    ///    .fam_path(&deb_maf_mib[1])
    ///    .bim_path(&deb_maf_mib[2])
    ///    .build()?;
    /// println!("{:?}", bed_cloud.iid()?); // Outputs ndarray ["iid1", "iid2", "iid3"]
    /// println!("{:?}", bed_cloud.sid()?); // Outputs ndarray ["sid1", "sid2", "sid3", "sid4"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
    /// ```
    // #[anyinput]
    pub fn bim_path(mut self, path: StorePath) -> Self {
        self.bim_path = Some(Some(path.to_owned()));
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
    /// use bed_reader::{BedCloud, Metadata, sample_bed_file};
    ///
    /// let file_name = sample_bed_file("small.bed_cloud")?;
    /// let metadata = Metadata::builder()
    ///     .iid(["i1", "i2", "i3"])
    ///     .sid(["s1", "s2", "s3", "s4"])
    ///     .build()?;
    /// let mut bed_cloud = BedCloud::builder(file_name)
    ///     .fid(["f1", "f2", "f3"])
    ///     .iid(["x1", "x2", "x3"])
    ///     .metadata(&metadata)
    ///     .build()?;
    /// println!("{0:?}", bed_cloud.fid()?);  // Outputs ndarray ["f1", "f2", "f3"]
    /// println!("{0:?}", bed_cloud.iid()?);  // Outputs ndarray ["i1", "i2", "i3"]
    /// println!("{0:?}", bed_cloud.sid()?);  // Outputs ndarray ["s1", "s2", "s3", "s4"]
    /// println!("{0:?}", bed_cloud.chromosome()?);  // Outputs ndarray ["1", "1", "5", "Y"]
    /// # use bed_reader::BedErrorPlus;
    /// # Ok::<(), Box<BedErrorPlus>>(())
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
