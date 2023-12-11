use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::{collections::HashSet, path::PathBuf};

use futures_util::TryStreamExt;
use object_store::delimited::newline_delimited_stream;
use object_store::ObjectMeta;
use object_store::{path::Path, ObjectStore};

use crate::MetadataFields;
use crate::{BedErrorPlus, Metadata};

/// cmk doc
pub struct BedCloud<T>
where
    T: ObjectStore,
{
    pub(crate) store: Arc<T>,

    pub(crate) path: Path,

    pub(crate) object_meta: ObjectMeta,

    pub(crate) fam_path: Option<PathBuf>,

    // cmk needed? combine with path?
    pub(crate) fam_object_meta: Option<ObjectMeta>,

    pub(crate) bim_path: Option<PathBuf>,

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
    pub(crate) async fn count_lines(&self, path: &Path) -> Result<usize, Box<BedErrorPlus>> {
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
}
