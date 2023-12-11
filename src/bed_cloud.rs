use std::collections::HashSet;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use futures_util::TryStreamExt;
use object_store::delimited::newline_delimited_stream;
use object_store::path::Path as StorePath;
use object_store::ObjectMeta;
use object_store::ObjectStore;

use crate::MetadataFields;
use crate::{BedErrorPlus, Metadata};

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
        bed_path: StorePath,
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
    pub fn fam_path(&mut self) -> StorePath {
        // We need to clone the path because self might mutate later
        if let Some(path) = &self.fam_path {
            path.clone()
        } else {
            let path = &self.to_metadata_path(&self.fam_path, "fam");
            self.fam_path = Some(path.clone());
            path
        }
    }

    pub async fn iid_count(&mut self) -> Result<usize, Box<BedErrorPlus>> {
        if let Some(iid_count) = self.iid_count {
            Ok(iid_count)
        } else {
            let fam_path = self.fam_path();
            let iid_count = self.count_lines(&fam_path).await?;
            self.iid_count = Some(iid_count);
            Ok(iid_count)
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
    path_str.push_str(".");
    path_str.push_str(new_extension);

    // Parse the string back to StorePath
    StorePath::parse(&path_str).map_err(Into::into)
}
