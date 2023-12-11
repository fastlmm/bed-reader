use ndarray as nd;
use std::collections::HashSet;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

// cmk really need futures_util and bytes...???
use futures_util::TryStreamExt;
use object_store::delimited::newline_delimited_stream;
use object_store::path::Path as StorePath;
use object_store::ObjectMeta;
use object_store::ObjectStore;

use crate::{
    check_and_precompute_iid_index, set_up_two_bits_to_value, try_div_4, BedError, BedErrorPlus,
    BedVal, Metadata,
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
            change_extension(&bed_path, extension)
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
    // is_a1_counted: bool,
    // iid_index: &[isize],
    // sid_index: &[isize],
    // missing_value: TVal,
    // out_val: &mut nd::ArrayViewMut2<'_, TVal>, //mutable slices additionally allow to modify elements. But slices cannot grow - they are just a view into some vector.
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

    // // Check and precompute for each iid_index

    // let (i_div_4_array, i_mod_4_times_2_array) =
    //     check_and_precompute_iid_index(in_iid_count, iid_index)?;

    // // Check and compute work for each sid_index

    // let from_two_bits_to_value = set_up_two_bits_to_value(is_a1_counted, missing_value);
    // let lower_sid_count = -(in_sid_count as isize);
    // let upper_sid_count: isize = (in_sid_count as isize) - 1;
    // // See https://morestina.net/blog/1432/parallel-stream-processing-with-rayon
    // // Possible optimization: We could try to read only the iid info needed
    // // Possible optimization: We could read snp in their input order instead of their output order
    // sid_index
    //     .iter()
    //     .map(|in_sid_i_signed| {
    //         // Turn signed sid_index into unsigned sid_index (or error)
    //         let in_sid_i = if (0..=upper_sid_count).contains(in_sid_i_signed) {
    //             *in_sid_i_signed as u64
    //         } else if (lower_sid_count..=-1).contains(in_sid_i_signed) {
    //             (in_sid_count - ((-in_sid_i_signed) as usize)) as u64
    //         } else {
    //             return Err(Box::new(BedErrorPlus::BedError(BedError::SidIndexTooBig(
    //                 *in_sid_i_signed,
    //             ))));
    //         };

    //         // Read the iid info for one snp from the disk
    //         let mut bytes_vector: Vec<u8> = vec![0; in_iid_count_div4];
    //         let pos: u64 = in_sid_i * in_iid_count_div4_u64 + CB_HEADER_U64; // "as" and math is safe because of early checks
    //         buf_reader.seek(SeekFrom::Start(pos))?;
    //         buf_reader.read_exact(&mut bytes_vector)?;
    //         Ok(bytes_vector)
    //     })
    //     // Zip in the column of the output array
    //     .zip(out_val.axis_iter_mut(nd::Axis(1)))
    //     // In parallel, decompress the iid info and put it in its column
    //     .par_bridge() // This seems faster that parallel zip
    //     .try_for_each(|(bytes_vector_result, mut col)| match bytes_vector_result {
    //         Err(e) => Err(e),
    //         Ok(bytes_vector) => {
    //             for out_iid_i in 0..iid_index.len() {
    //                 let i_div_4 = i_div_4_array[out_iid_i];
    //                 let i_mod_4_times_2 = i_mod_4_times_2_array[out_iid_i];
    //                 let genotype_byte: u8 = (bytes_vector[i_div_4] >> i_mod_4_times_2) & 0x03;
    //                 col[out_iid_i] = from_two_bits_to_value[genotype_byte as usize];
    //             }
    //             Ok(())
    //         }
    //     })?;

    Ok(())
}
