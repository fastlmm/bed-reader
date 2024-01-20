# Cloud URLs and ObjectPath Examples

cmk update with this about how splitting url in two:
    let option_store = HttpBuilder::new()
        .with_url("https://www.ebi.ac.uk/")
        .build()?;
    let store_path =
        StorePath::parse("biostudies/files/S-BSST936/example/synthetic_small_v1_chr-10.bed")?;
    let object_path = ObjectPath::new(Arc::new(option_store), store_path);

To specify a file in the cloud, you must specify either

* a URL string plus options, or
* an [`ObjectPath`](../struct.ObjectPath.html)

<!-- cmk: and http -->
Let's look at how to do this [for a local file](#local-file) and [for AWS S3](#aws-s3).

## Local File

We can specify a local file as if it is in the cloud. This is a great way to test our cloud functions.

### Local File URL

The URL for a local file takes the form `file:///{encoded_file_name}`. No cloud options are needed, so we use `EMPTY_OPTIONS`.

*Example:*

```rust
use ndarray as nd;
use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_file, EMPTY_OPTIONS, path_to_url_string};
# use {bed_reader::BedErrorPlus, tokio::runtime::Runtime}; // '#' needed for doctest
# Runtime::new().unwrap().block_on(async {

let file_name = sample_bed_file("small.bed")?.to_string_lossy().to_string();
println!("{file_name:?}"); // For example, "C:\\Users\\carlk\\AppData\\Local\\fastlmm\\bed-reader\\cache\\small.bed"
let url: String = path_to_url_string(file_name)?;
println!("{url:?}"); // For example, "file:///C:/Users/carlk/AppData/Local/bed_reader/bed_reader/Cache/small.bed"

let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
let val = ReadOptions::builder().sid_index(2).f64().read_cloud(&mut bed_cloud).await?;
assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
# Ok::<(), Box<dyn std::error::Error>>(())
# }).unwrap();
```

### Local File ObjectPath

If we want to work with structs instead of `String`, we can specify a file via an [`ObjectPath`](../struct.ObjectPath.html).

We first create an [`ObjectStore`](https://docs.rs/object_store/latest/object_store/trait.ObjectStore.html) that tells which cloud service we are using.
Here we use [`LocalFileSystem`](https://docs.rs/object_store/latest/object_store/local/struct.LocalFileSystem.html).
We, next, wrap the result in an `Arc` to faciliate efficient cloning of something that would otherwise be un-clonable.
We then create an [`object_store::path::Path as StorePath`](https://docs.rs/object_store/latest/object_store/path/struct.Path.html)
using [`StorePath::from_filesystem_path`](https://docs.rs/object_store/latest/object_store/path/struct.Path.html#method.from_filesystem_path).
It tells where on the cloud service your file is located.
Finally, we put the `ObjectStore` and `StorePath` together, creating an [`ObjectPath`](../struct.ObjectPath.html).

[`BedCloud`](../struct.BedCloud.html) can work via an `ObjectPath` slightly more efficiently than a string URL.
However, I generally use string URLs because I find them easier.

This example puts all the steps together:

```rust
use std::sync::Arc;
use ndarray as nd;
use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_file, ObjectPath};
use object_store::{local::LocalFileSystem, path::Path as StorePath};
# use {bed_reader::BedErrorPlus, tokio::runtime::Runtime}; // '#' needed for doctest
# Runtime::new().unwrap().block_on(async {
let file_name = sample_bed_file("small.bed")?.to_string_lossy().to_string();

let arc_object_store = Arc::new(LocalFileSystem::new());
let path = StorePath::from_filesystem_path(&file_name)?;
let object_path = ObjectPath::new(arc_object_store, path);

let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
let val = ReadOptions::builder().sid_index(2).f64().read_cloud(&mut bed_cloud).await?;
assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
# Ok::<(), Box<dyn std::error::Error>>(())
# }).unwrap();
```

## AWS S3

Let's look next at a more useful scenerio: Reading a file (or part of a file) from AWS S3.

### AWS S3 URL

The URL for an AWS S3 file takes the form `s3://{bucket_name}/{s3_path}`.

AWS forbids putting some needed information in the URL. Instead, that information must
go into a string-to-string map of options. Specifically, we'll put `"aws_region"`, `"aws_access_key_id"`, and `"aws_secret_access_key"` in the options.
For security, we pull the last two option values from a file rather than hard-coding them into the program.

*Example:*

> **Note:** I can run this, but others can't because of the authentication checks.

```rust
use bed_reader::{BedCloud,BedErrorPlus};
use tokio::runtime::Runtime;
use rusoto_credential::{CredentialsError, ProfileProvider, ProvideAwsCredentials};

Runtime::new().unwrap().block_on(async {
    // Read my AWS credentials from file ~/.aws/credentials
    let credentials = if let Ok(provider) = ProfileProvider::new() {
        provider.credentials().await
    } else {
        Err(CredentialsError::new("No credentials found"))
    };

    let Ok(credentials) = credentials else {
        eprintln!("Skipping test because no AWS credentials found");
        return Ok(());
    };

    let url = "s3://bedreader/v1/toydata.5chrom.bed";
    let options = [
        ("aws_region", "us-west-2"),
        ("aws_access_key_id", credentials.aws_access_key_id()),
        ("aws_secret_access_key", credentials.aws_secret_access_key()),
    ];

    let mut bed_cloud = BedCloud::new(url, options).await?;
    let val = bed_cloud.read::<i8>().await?;
    assert_eq!(val.shape(), &[500, 10_000]);
    Ok::<(), Box<BedErrorPlus>>(())
});
Ok::<(), Box<BedErrorPlus>>(())
```

### AWS S3 ObjectPath

Again, suppose we want to work with structs instead of `String`. Again, we can specify a file via an [`ObjectPath`](../struct.ObjectPath.html).

We first create an [`ObjectStore`](https://docs.rs/object_store/latest/object_store/trait.ObjectStore.html) that tells which cloud service we are using.
We use [`AmazonS3Builder`](https://docs.rs/object_store/latest/object_store/aws/struct.AmazonS3Builder.html).
We again wrap the result in an `Arc` to faciliate efficient cloning of something that would otherwise be un-clonable.
We then create an [`object_store::path::Path as StorePath`](https://docs.rs/object_store/latest/object_store/path/struct.Path.html)
using [`StorePath::parse`](https://docs.rs/object_store/latest/object_store/path/struct.Path.html#method.parse).
Finally, we put the `ObjectStore` and `StorePath` together, creating an [`ObjectPath`](../struct.ObjectPath.html).

As before, [`BedCloud`](../struct.BedCloud.html) can work via an `ObjectPath` slightly more efficiently than a string URL.
However, I generally use string URLs because I find them easier.

Here is an AWS Object Path example:

> **Note:** I can run this, but others can't because of the authentication checks.

```rust
use std::sync::Arc;
use ndarray as nd;
use bed_reader::{BedCloud, ObjectPath};
use object_store::{aws::AmazonS3Builder, path::Path as StorePath};
# use {bed_reader::BedErrorPlus, tokio::runtime::Runtime}; // '#' needed for doctest
use rusoto_credential::{CredentialsError, ProfileProvider, ProvideAwsCredentials};

# Runtime::new().unwrap().block_on(async {
// Read my AWS credentials from file ~/.aws/credentials
let credentials = if let Ok(provider) = ProfileProvider::new() {
    provider.credentials().await
} else {
    Err(CredentialsError::new("No credentials found"))
};

let Ok(credentials) = credentials else {
    eprintln!("Skipping test because no AWS credentials found");
    return Ok(());
};

let arc_s3 = Arc::new(
    AmazonS3Builder::new()
        .with_region("us-west-2")
        .with_bucket_name("bedreader")
        .with_access_key_id(credentials.aws_access_key_id())
        .with_secret_access_key(credentials.aws_secret_access_key())
        .build()?,
);
let store_path = StorePath::parse("/v1/toydata.5chrom.bed")?;
let object_path = ObjectPath::new(arc_s3, store_path);

let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
let val = bed_cloud.read::<i8>().await?;
assert_eq!(val.shape(), &[500, 10_000]);
# Ok::<(), Box<dyn std::error::Error>>(())
# }).unwrap();
```
