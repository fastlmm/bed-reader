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
use object_store::ObjectMeta;
use object_store::ObjectStore;

use crate::{
    check_and_precompute_iid_index, compute_max_chunk_size, compute_max_concurrent_requests,
    set_up_two_bits_to_value, try_div_4, BedError, BedErrorPlus, BedVal, Hold, Metadata,
    ReadOptions, BED_FILE_MAGIC1, BED_FILE_MAGIC2,
};
use crate::{MetadataFields, CB_HEADER_U64};

/// cmk doc
pub struct BedCloud<T>
where
    T: ObjectStore,
{
    pub(crate) store: Arc<T>,

    pub(crate) path: StorePath,

    pub(crate) object_meta: ObjectMeta,

    pub(crate) fam_path: Option<StorePath>,

    // cmk needed? combine with path?
    pub(crate) fam_object_meta: Option<ObjectMeta>,

    pub(crate) bim_path: Option<StorePath>,

    pub(crate) bim_object_meta: Option<ObjectMeta>,

    pub(crate) is_checked_early: bool,

    pub(crate) iid_count: Option<usize>,

    pub(crate) sid_count: Option<usize>,

    pub(crate) metadata: Metadata,

    pub(crate) skip_set: HashSet<MetadataFields>,
}

impl<T> BedCloud<T>
where
    T: ObjectStore,
{
    // #[anyinput]
    pub(crate) async fn count_lines(&self, path: &StorePath) -> Result<usize, Box<BedErrorPlus>> {
        let stream = self
            .store
            .clone()
            // cmk !!!! does this get get the whole file or just the first chunk?
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
            let path = BedCloud::<T>::to_metadata_path(&self.path, &self.fam_path, "fam")?;
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
            let path = BedCloud::<T>::to_metadata_path(&self.path, &self.bim_path, "bim")?;
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
            self.store.clone(), // cmk rename "object_store" ??
            self.path.clone(),
            &self.object_meta,
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
        .await?;

        Ok(())
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
pub(crate) async fn internal_read_no_alloc<TVal: BedVal, TStore: ObjectStore>(
    object_store: Arc<TStore>,
    path: StorePath,
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
            BedError::IllFormed(store_path_to_string(path)).into(),
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
pub(crate) async fn read_no_alloc<TVal: BedVal, TStore: ObjectStore>(
    object_store: Arc<TStore>,
    path: StorePath,
    object_meta: &ObjectMeta,
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
    let bytes = open_and_check(object_store.clone(), &path).await?;

    match bytes[2] {
        0 => {
            // We swap 'iid' and 'sid' and then reverse the axes.
            let mut val_t = val.view_mut().reversed_axes();

            internal_read_no_alloc(
                object_store,
                path,
                object_meta,
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
                object_meta,
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
                BedError::BadMode(store_path_to_string(path)).into(),
            ))
        }
    };
    Ok(())
}

// cmk #[anyinput]
async fn open_and_check<TStore>(
    object_store: Arc<TStore>,
    path: &StorePath,
) -> Result<Bytes, Box<BedErrorPlus>>
where
    TStore: ObjectStore,
{
    let bytes = object_store
        .clone()
        .get_range(path, 0..CB_HEADER_U64 as usize)
        .await
        .map_err(BedErrorPlus::from)?;

    if (BED_FILE_MAGIC1 != bytes[0]) || (BED_FILE_MAGIC2 != bytes[1]) {
        return Err(Box::new(
            BedError::IllFormed(store_path_to_string(path.clone())).into(),
        ));
    }
    Ok(bytes)
}
