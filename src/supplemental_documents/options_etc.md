# Options, Options, Options

Within this crate, the term "options" can refer to three levels of options: [Cloud](#cloud-options), [Bed/BedCloud](#bedbedcloud-options), and [ReadOptions](#readoptions).

## Cloud options

When specifying a file in the cloud via a URL, we use methods [`BedCloud::new(url)`](../struct.BedCloud.html#method.new),
[`BedCloud::new_with_options(url, options)`](../struct.BedCloud.html#method.new_with_options), and
[`BedCloud::builder(url, options)`](../struct.BedCloud.html#method.builder).

The cloud providers forbid putting some needed information in the URL. Instead, that information must
go into `options`. For example, AWS S3 requires that information
about `"aws_region"`, `"aws_access_key_id"`, and `"aws_secret_access_key"` be placed in the options.

See [`ClientConfigKey`](https://docs.rs/object_store/latest/object_store/enum.ClientConfigKey.html) for a list of cloud options, such as `timeout`, that you can always use. See [`AmazonS3ConfigKey`](https://docs.rs/object_store/latest/object_store/aws/enum.AmazonS3ConfigKey.html) for a list of AWS-specific options.
See [`AzureConfigKey`](https://docs.rs/object_store/latest/object_store/azure/enum.AzureConfigKey.html) for a list of Azure-specific options.
See [`GoogleConfigKey`](https://docs.rs/object_store/latest/object_store/gcp/enum.GoogleConfigKey.html) for a list of Google-specific options.

Here is an AWS example:

> **Note:** I can run this, but others can't because of the authentication checks.

```rust
use bed_reader::{BedCloud,BedErrorPlus};
use rusoto_credential::{CredentialsError, ProfileProvider, ProvideAwsCredentials};

# use tokio::runtime::Runtime;
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

let url = "s3://bedreader/v1/toydata.5chrom.bed";
let options = [
    ("aws_region", "us-west-2"),
    ("aws_access_key_id", credentials.aws_access_key_id()),
    ("aws_secret_access_key", credentials.aws_secret_access_key()),
];

let mut bed_cloud = BedCloud::new_with_options(url, options).await?;
let val = bed_cloud.read::<i8>().await?;
assert_eq!(val.shape(), &[500, 10_000]);
# Ok::<(), Box<dyn std::error::Error>>(())
# }).unwrap();
```

We can also read from local files as though they are in the cloud. In that case, no cloud options are needed, so we use `EMPTY_OPTIONS`:

```rust
use ndarray as nd;
use bed_reader::{BedCloud, ReadOptions, assert_eq_nan, sample_url};
# use {bed_reader::BedErrorPlus, tokio::runtime::Runtime}; // '#' needed for doctest
# Runtime::new().unwrap().block_on(async {
let url = sample_url("small.bed")?;
println!("{url:?}"); // For example, "file:///C:/Users/carlk/AppData/Local/bed_reader/bed_reader/Cache/small.bed"
let options = EMPTY_OPTIONS; // map of authentication keys, etc., if needed.
let mut bed_cloud = BedCloud::new_with_options(url, options).await?;
let val = ReadOptions::builder().sid_index(2).f64().read_cloud(&mut bed_cloud).await?;
assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
# Ok::<(), Box<dyn std::error::Error>>(())
# }).unwrap();
```

Other cloud services, for example, Azure and Google Cloud, also need cloud options. Their options are similar to, but not identical to, the options for AWS S3. You will need to research the details to use any cloud service.

## Bed/BedCloud-Options

When you open a local file for reading, you can set options via [`Bed::builder`](../struct.Bed.html#method.builder). When you open a cloud file for reading, you can set options via [`BedCloud::builder`](../struct.BedCloud.html#method.builder).

The options, [listed here](../struct.BedBuilder.html#implementations), can:

* set the path of the .fam and/or .bim file
* override some metadata, for example, replace the individual ids.
* set the number of individuals (samples) or SNPs (variants)
* control checking the validity of the .bed file’s header
* skip reading selected metadata

For example, here we replace the `iid` in the file with our own list:

```rust
use bed_reader::{Bed, sample_bed_file};
let file_name = sample_bed_file("small.bed")?;
let mut bed = Bed::builder(file_name)
   .iid(["sample1", "sample2", "sample3"])
   .build()?;
println!("{:?}", bed.iid()?); // Outputs ndarray ["sample1", "sample2", "sample3"]
# use bed_reader::BedErrorPlus;
# Ok::<(), Box<BedErrorPlus>>(())
```

## ReadOptions

When reading read genotype data, use [`ReadOptions::builder`](../struct.ReadOptions.html#method.builder) to specify:

* a desired numeric type,
* which individuals (samples) to read,
* which SNPs (variants) to read,
* etc.

See this crate's introductory material for [a complete table](../index.html#readoptions).
