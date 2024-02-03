# Cloud URL Examples

> *Table of Contents*:
>
> * [Http](#http)
> * [local file](#local-file)
> * [AWS S3](#aws-s3)
>
> *The `bed-reader` crate also supports Azure and GCP, but we don't have examples.*

To specify a file in the cloud, you must specify a URL string plus options.

The exact details depend on the cloud service. We'll look at [http](#http), at [local files](#local-file), and at [AWS S3](#aws-s3).

## Http

You can read \*.bed files from web sites directly. For small files, access will be fast.
For medium-sized files, you may need to extend the default `timeout`.

Reading from large files can also be practical and even fast under these conditions:

* You need only some of the information
* (Optional, but helpful) You can provide some metadata about individuals (samples) and SNPs (variants) locally.

Let's first look at reading a small or medium-sized dataset using a URL string.

*Example:*

Read an entire file and find the fraction of missing values.

```rust
use ndarray as nd;
use bed_reader::{BedCloud};

# #[cfg(feature = "tokio")] Runtime::new().unwrap().block_on(async { // '#' needed for doctest
let url = "https://raw.githubusercontent.com/fastlmm/bed-sample-files/main/small.bed";
let mut bed_cloud = BedCloud::new(url).await?;
let val: nd::Array2<f32> = bed_cloud.read().await?;
let missing_count = val.iter().filter(|x| x.is_nan()).count();
let missing_fraction = missing_count as f32 / val.len() as f32;
println!("{missing_fraction:.2}"); // Outputs 0.17
assert_eq!(missing_count, 2);
# Ok::<(), Box<dyn std::error::Error>>(()) }).unwrap()};  // '#' needed for doctest
# #[cfg(feature = "tokio")] use {tokio::runtime::Runtime, bed_reader::BedErrorPlus}; // '#' needed for doctest
```

When reading a medium-sized file, you may need to set a `timeout` in your options.
With a `timeout`, you can give your code more than the default 30 seconds to read
metadata from the \*.fam and \*.bim files (or genomic data from \*.bed).

See [`ClientConfigKey`](https://docs.rs/object_store/latest/object_store/enum.ClientConfigKey.html) for a list of cloud options, such as `timeout`,
that you can use with `Http`.

You may also wish to use [`.skip_early_check()`](../struct.BedCloudBuilder.html#method.skip_early_check)
to avoid a fast, early check of the \*.bed file's header.

Here we print the first five iids (individual or sample ids) and first five sids (SNP or variant ids).
We then, print all unique chromosome values. Finally, we read all data from chromosome 5 and print its dimensions.

```rust
use ndarray::{self as nd, s};
use bed_reader::{BedCloud, ReadOptions};
use std::collections::BTreeSet;

# use {bed_reader::BedErrorPlus, tokio::runtime::Runtime}; // '#' needed for doctest
# Runtime::new().unwrap().block_on(async {
let mut bed_cloud = BedCloud::builder(
    "https://raw.githubusercontent.com/fastlmm/bed-sample-files/main/toydata.5chrom.bed",
    [("timeout", "100s")],
)?.skip_early_check().build().await?;
println!("{:?}", bed_cloud.iid().await?.slice(s![..5])); // Outputs ndarray: ["per0", "per1", "per2", "per3", "per4"]
println!("{:?}", bed_cloud.sid().await?.slice(s![..5])); // Outputs ndarray: ["null_0", "null_1", "null_2", "null_3", "null_4"]
println!("{:?}", bed_cloud.chromosome().await?
            .iter().collect::<BTreeSet<_>>()); // Outputs: {"1", "2", "3", "4", "5"}
let val = ReadOptions::builder()
    .sid_index(bed_cloud.chromosome().await?.map(|elem| elem == "5"))
    .f32()
    .read_cloud(&mut bed_cloud)
    .await?;
assert_eq!(val.dim(), (500, 440));
# Ok::<(), Box<dyn std::error::Error>>(())
# }).unwrap();
```

Now, let's read from a large file containing data from over 1 million individuals (samples) and over 300,000
SNPs (variants). The file size is 91 GB. In this example, we read data for just one SNP (variant). If
we know the number of individuals (samples) and SNPs (variants) exactly, we can read this SNP quickly and with
just one file access.

What is the mean value of the SNP (variant) at index position 100,000?

```rust
use ndarray as nd;
use bed_reader::{BedCloud, ReadOptions};

# use {bed_reader::BedErrorPlus, tokio::runtime::Runtime}; // '#' needed for doctest
# Runtime::new().unwrap().block_on(async {
let mut bed_cloud = BedCloud::builder(
        "https://www.ebi.ac.uk/biostudies/files/S-BSST936/genotypes/synthetic_v1_chr-10.bed",
        EMPTY_OPTIONS)?
    .skip_early_check()
    .iid_count(1_008_000)
    .sid_count(361_561)
    .build()
    .await?;
let val = ReadOptions::builder()
    .sid_index(100_000)
    .f32()
    .read_cloud(&mut bed_cloud)
    .await?;
assert_eq!(val.mean(), Some(0.03391369));
# Ok::<(), Box<dyn std::error::Error>>(())
# }).unwrap();
```

## Local File

We can specify a local file as if it is in the cloud. This is a great way to test cloud functions. For real work and better efficiency,
however, use [`Bed`](../struct.Bed.html) instead of [`BedCloud`](../struct.BedCloud.html).

### Local File URL

The URL for a local file takes the form `file:///{encoded_file_name}`. We can use the [`abs_path_to_url_string`](../fn.abs_path_to_url_string.html)
function to do this encoding. When it was the url to `BedCloud`,  no cloud options are needed, so we use `EMPTY_OPTIONS`.

*Example:*

```rust
use ndarray as nd;
use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_bed_file};
use cloud_file::abs_path_to_url_string;

# #[cfg(feature = "tokio")] Runtime::new().unwrap().block_on(async { // '#' needed for doc test
let file_name = sample_bed_file("small.bed")?;
println!("{file_name:?}"); // For example, "C:\\Users\\carlk\\AppData\\Local\\fastlmm\\bed-reader\\cache\\small.bed"
let url: String = abs_path_to_url_string(file_name)?;
println!("{url:?}"); // For example, "file:///C:/Users/carlk/AppData/Local/bed_reader/bed_reader/Cache/small.bed"

let mut bed_cloud = BedCloud::new(url).await?;
let val = ReadOptions::builder().sid_index(2).f64().read_cloud(&mut bed_cloud).await?;
assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
# Ok::<(), Box<dyn std::error::Error>>(()) }).unwrap()};  // '#' needed for doctest
# #[cfg(feature = "tokio")] use {tokio::runtime::Runtime, bed_reader::BedErrorPlus}; // '#' needed for doctest
```

## AWS S3

Let's look next at reading a file (or part of a file) from AWS S3.

The URL for an AWS S3 file takes the form `s3://{bucket_name}/{s3_path}`.

AWS forbids putting some needed information in the URL. Instead, that information must
go into a string-to-string map of options. Specifically, we'll put `"aws_region"`, `"aws_access_key_id"`, and `"aws_secret_access_key"` in the options.
For security, we pull the last two option values from a file rather than hard-coding them into the program.

See [`ClientConfigKey`](https://docs.rs/object_store/latest/object_store/enum.ClientConfigKey.html) for a list of cloud options, such as `timeout`,
that you can always use.
See [`AmazonS3ConfigKey`](https://docs.rs/object_store/latest/object_store/aws/enum.AmazonS3ConfigKey.html) for a list of AWS-specific options.
See [`AzureConfigKey`](https://docs.rs/object_store/latest/object_store/azure/enum.AzureConfigKey.html) for a list of Azure-specific options.
See [`GoogleConfigKey`](https://docs.rs/object_store/latest/object_store/gcp/enum.GoogleConfigKey.html) for a list of Google-specific options.

*Example:*

> **Note:** I can run this, but others can't because of the authentication checks.

```rust
use bed_reader::{BedCloud,BedErrorPlus};
use rusoto_credential::{CredentialsError, ProfileProvider, ProvideAwsCredentials};

# #[cfg(feature = "tokio")] Runtime::new().unwrap().block_on(async { // '#' needed for doc test
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

let mut bed_cloud = BedCloud::new_with_options(url, options).await?;
let val = bed_cloud.read::<i8>().await?;
assert_eq!(val.shape(), &[500, 10_000]);
# Ok::<(), Box<dyn std::error::Error>>(()) }).unwrap()};
# #[cfg(feature = "tokio")] use {tokio::runtime::Runtime, bed_reader::BedErrorPlus};
```
