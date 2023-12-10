use std::sync::Arc;
use std::{collections::HashSet, path::PathBuf};

use object_store::buffered::BufReader;
use object_store::{path::Path, ObjectStore};

use crate::Metadata;
use crate::MetadataFields;

/// cmk doc
pub struct BedCloud<T>
where
    T: ObjectStore,
{
    pub(crate) store: Arc<T>,

    pub(crate) path: Path,

    pub(crate) fam_path: Option<PathBuf>,

    pub(crate) bim_path: Option<PathBuf>,

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
    async fn count_lines(&self, path: Path) -> Result<usize, anyhow::Error> {
        let data = self.store.get(&path).await?;

        // Create an async reader from the bytes stream
        let reader = BufReader::new(data);

        // Count lines
        let mut line_count = 0;
        let mut line = String::new();
        while reader.read_line(&mut line).await? != 0 {
            line_count += 1;
            line.clear();
        }

        Ok(line_count)
    }
}
