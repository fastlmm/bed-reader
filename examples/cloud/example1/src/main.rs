use anyhow::anyhow;
use futures::pin_mut;
use futures::StreamExt;
use object_store::http::HttpBuilder;
use object_store::{delimited::newline_delimited_stream, path::Path as StorePath};
use object_store::{GetResult, ObjectStore};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::sync::Arc;
use url::Url;

async fn count_lines<TObjectStore>(
    object_path: &ObjectPath<TObjectStore>,
) -> Result<usize, anyhow::Error>
where
    TObjectStore: ObjectStore,
{
    let mut stream = object_path.get().await?.into_stream();
    let mut newline_count: usize = 0;
    while let Some(bytes) = stream.next().await {
        let bytes = bytes?;
        let count = bytecount::count(&bytes, b'\n');
        newline_count += count;
    }

    Ok(newline_count)
}

async fn read_random_line(
    object_path: &ObjectPath<&dyn ObjectStore>,
    seed: u64,
) -> Result<String, anyhow::Error> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut selected_line = None;

    let stream = object_path.get().await?.into_stream();
    let line_chunk_stream = newline_delimited_stream(stream);
    pin_mut!(line_chunk_stream);

    while let Some(line_chunk) = line_chunk_stream.next().await {
        let line_chunk = line_chunk?;
        let lines = std::str::from_utf8(&line_chunk)?.split_terminator('\n');

        for (index, line) in lines.enumerate() {
            // Reservoir sampling: replace the selected line with probability 1/(index+1)
            if rng.gen_range(0..=index) == 0 {
                selected_line = Some(line.to_string());
            }
        }
    }

    selected_line.ok_or_else(|| anyhow!("No lines found in the file"))
}

pub struct ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    /// An `Arc`-wrapped [`ObjectStore`](https://docs.rs/object_store/latest/object_store/trait.ObjectStore.html) cloud service, for example, Http, AWS S3,
    /// Azure, the local file system, etc.
    pub object_store: Arc<TObjectStore>,
    /// A [`object_store::path::Path as StorePath`](https://docs.rs/object_store/latest/object_store/path/struct.Path.html) that points to a file on
    /// the [`ObjectStore`](https://docs.rs/object_store/latest/object_store/trait.ObjectStore.html)
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

impl<TObjectStore> ObjectPath<TObjectStore>
where
    TObjectStore: ObjectStore,
{
    pub fn new(arc_object_store: Arc<TObjectStore>, path: StorePath) -> Self {
        ObjectPath {
            object_store: arc_object_store,
            path,
        }
    }

    /// Return the bytes that are stored at the specified location.
    pub async fn get(&self) -> Result<GetResult, anyhow::Error> {
        Ok(self.object_store.get(&self.path).await?)
    }
}

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
        // cmk should we call this 'new'?
        location: S,
        options: I,
    ) -> Result<ObjectPath<Box<dyn ObjectStore>>, anyhow::Error>
    where
        I: IntoIterator<Item = (K, V)>,
        K: AsRef<str>,
        V: Into<String>,
        S: AsRef<str>,
    {
        let location = location.as_ref();
        let url =
            Url::parse(location).map_err(|e| anyhow!("Cannot parse url: {} {}", location, e))?;

        let (object_store, store_path): (Box<dyn ObjectStore>, StorePath) =
            parse_url_opts_work_around(&url, options)?;
        let object_path = ObjectPath::new(Arc::new(object_store), store_path);
        Ok(object_path)
    }
}

#[allow(clippy::match_bool)]
fn parse_work_around(url: &Url) -> Result<(bool, StorePath), object_store::Error> {
    let strip_bucket = || Some(url.path().strip_prefix('/')?.split_once('/')?.1);

    let (scheme, path) = match (url.scheme(), url.host_str()) {
        ("http", Some(_)) => (true, url.path()),
        ("https", Some(host)) => {
            if host.ends_with("dfs.core.windows.net")
                || host.ends_with("blob.core.windows.net")
                || host.ends_with("dfs.fabric.microsoft.com")
                || host.ends_with("blob.fabric.microsoft.com")
            {
                (false, url.path())
            } else if host.ends_with("amazonaws.com") {
                match host.starts_with("s3") {
                    true => (false, strip_bucket().unwrap_or_default()),
                    false => (false, url.path()),
                }
            } else if host.ends_with("r2.cloudflarestorage.com") {
                (false, strip_bucket().unwrap_or_default())
            } else {
                (true, url.path())
            }
        }
        _ => (false, url.path()),
    };

    Ok((scheme, StorePath::from_url_path(path)?))
}

// LATER when https://github.com/apache/arrow-rs/issues/5310 gets fixed, can remove work around
pub fn parse_url_opts_work_around<I, K, V>(
    url: &Url,
    options: I,
) -> Result<(Box<dyn ObjectStore>, StorePath), object_store::Error>
where
    I: IntoIterator<Item = (K, V)>,
    K: AsRef<str>,
    V: Into<String>,
{
    let (is_http, path) = parse_work_around(url)?;
    if is_http {
        let url = &url[..url::Position::BeforePath];
        let path = StorePath::parse(path)?;
        let builder = options.into_iter().fold(
            <HttpBuilder>::new().with_url(url),
            |builder, (key, value)| match key.as_ref().parse() {
                Ok(k) => builder.with_config(k, value),
                Err(_) => builder,
            },
        );
        let store = Box::new(builder.build()?) as _;
        Ok((store, path))
    } else {
        object_store::parse_url_opts(url, options)
    }
}
pub const EMPTY_OPTIONS: [(&str, String); 0] = [];

#[tokio::main]
async fn main() -> Result<(), anyhow::Error> {
    let object_path = ObjectPath::from_url(
        "https://raw.githubusercontent.com/fastlmm/bed-sample-files/main/small.fam",
        EMPTY_OPTIONS,
    )?;
    let line_count = count_lines(&object_path).await?;
    println!("line_count: {}", line_count);
    Ok(())
}
