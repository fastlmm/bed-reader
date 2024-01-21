#![cfg(feature = "cloud")]

use bed_reader::allclose;
use bed_reader::assert_eq_nan;
use bed_reader::assert_error_variant;
use bed_reader::path_to_url_string;
use bed_reader::sample_bed_file;
use bed_reader::sample_bed_object_path;
use bed_reader::sample_bed_url;
use bed_reader::sample_file;
use bed_reader::sample_object_path;
use bed_reader::sample_object_paths;
use bed_reader::BedCloud;
use bed_reader::BedError;
use bed_reader::BedErrorPlus;
use bed_reader::Metadata;
use bed_reader::MetadataFields;
use bed_reader::ObjectPath;
use bed_reader::ReadOptions;
use bed_reader::SliceInfo1;
use bed_reader::EMPTY_OPTIONS;
use bed_reader::{sample_url, sample_urls};
use ndarray as nd;
use ndarray::s;
use object_store::aws::AmazonS3Builder;
use object_store::local::LocalFileSystem;
use object_store::path::Path as StorePath;
use object_store::ObjectStore;
use std::collections::HashSet;
use std::panic::catch_unwind;
use std::sync::Arc;
use thousands::Separable;
use tokio::runtime;
use url::Url;

#[tokio::test]
async fn rusty_cloud_bed0() -> Result<(), Box<BedErrorPlus>> {
    let file_path = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?;
    use bed_reader::Bed;

    let mut bed = Bed::new(&file_path)?;
    let val0 = bed.read::<i8>()?;
    println!("{val0:?}");

    let url = sample_bed_url("plink_sim_10s_100v_10pmiss.bed")?;

    let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
    let val = bed_cloud.read::<i8>().await?;
    println!("{val:?}");

    assert!(allclose(&val0.view(), &val.view(), 0, true));
    Ok(())
}

#[tokio::test]
async fn rusty_cloud_bed0_url() -> Result<(), Box<BedErrorPlus>> {
    let file_path = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?;
    use bed_reader::Bed;

    let mut bed = Bed::new(&file_path)?;
    let val0 = bed.read::<i8>()?;
    println!("{val0:?}");

    let url = sample_url("plink_sim_10s_100v_10pmiss.bed")?;
    let object_path = ObjectPath::from_url(&url, EMPTY_OPTIONS)?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    let val = bed_cloud.read::<i8>().await?;
    println!("{val:?}");
    assert!(allclose(&val0.view(), &val.view(), 0, true));

    let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
    let val = bed_cloud.read::<i8>().await?;
    println!("{val:?}");
    assert!(allclose(&val0.view(), &val.view(), 0, true));
    Ok(())
}

#[tokio::test]
async fn rusty_cloud_bed1() -> Result<(), Box<BedErrorPlus>> {
    let file_path = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?;
    use bed_reader::Bed;

    let mut bed = Bed::new(&file_path)?;
    let val = bed.read::<i8>()?;
    println!("{val:?}");
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    assert!(mean == -13.142); // really shouldn't do mean on data where -127 represents missing

    let mut bed = Bed::new(&file_path)?;
    let val = ReadOptions::builder().count_a2().i8().read(&mut bed)?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    assert!(mean == -13.274); // really shouldn't do mean on data where -127 represents missing

    let url = path_to_url_string(file_path)?;

    let mut bed_cloud = BedCloud::new(&url, EMPTY_OPTIONS).await?;
    let val = bed_cloud.read::<i8>().await?;
    println!("{val:?}");
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    assert!(mean == -13.142); // really shouldn't do mean on data where -127 represents missing

    let mut bed_cloud = BedCloud::new(&url, EMPTY_OPTIONS).await?;
    let val = ReadOptions::builder()
        .count_a2()
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    assert!(mean == -13.274); // really shouldn't do mean on data where -127 represents missing
    Ok(())
}

#[tokio::test]
async fn rusty_cloud_bed1_url() -> Result<(), Box<BedErrorPlus>> {
    let url = sample_bed_url("plink_sim_10s_100v_10pmiss.bed")?;

    let mut bed_cloud = BedCloud::new(&url, EMPTY_OPTIONS).await?;
    let val = bed_cloud.read::<i8>().await?;
    println!("{val:?}");
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    assert!(mean == -13.142); // really shouldn't do mean on data where -127 represents missing

    let mut bed_cloud = BedCloud::new(&url, EMPTY_OPTIONS).await?;
    let val = ReadOptions::builder()
        .count_a2()
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    assert!(mean == -13.274); // really shouldn't do mean on data where -127 represents missing
    Ok(())
}

#[tokio::test]
async fn rusty_cloud_bed2() -> Result<(), Box<BedErrorPlus>> {
    use bed_reader::Bed;

    let file = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?;

    let mut bed = Bed::new(&file)?;
    let val = ReadOptions::builder()
        .iid_index(0)
        .sid_index(vec![1])
        .i8()
        .read(&mut bed)?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == 1.0); // really shouldn't do mean on data where -127 represents missing

    let url = path_to_url_string(file)?;
    let mut bed_cloud = BedCloud::new(&url, EMPTY_OPTIONS).await?;
    let val = ReadOptions::builder()
        .iid_index(0)
        .sid_index(vec![1])
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == 1.0); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[tokio::test]
async fn rusty_cloud_bed2_url() -> Result<(), Box<BedErrorPlus>> {
    let url = sample_bed_url("plink_sim_10s_100v_10pmiss.bed")?;
    let mut bed_cloud = BedCloud::new(&url, EMPTY_OPTIONS).await?;
    let val = ReadOptions::builder()
        .iid_index(0)
        .sid_index(vec![1])
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == 1.0); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[tokio::test]
async fn rusty_cloud_bed3() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("plink_sim_10s_100v_10pmiss.bed")?;

    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    let iid_bool: nd::Array1<bool> = (0..bed_cloud.iid_count().await?)
        .map(|elem| (elem % 2) != 0)
        .collect();
    let sid_bool: nd::Array1<bool> = (0..bed_cloud.sid_count().await?)
        .map(|elem| (elem % 8) != 0)
        .collect();
    let val = ReadOptions::builder()
        .missing_value(-127)
        .iid_index(iid_bool)
        .sid_index(sid_bool)
        .read_cloud(&mut bed_cloud)
        .await?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == -14.50344827586207); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[tokio::test]
async fn rusty_cloud_bed_allele() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("plink_sim_10s_100v_10pmiss.bed")?;

    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    let val = ReadOptions::builder()
        .count_a2()
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;

    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == -13.274); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[tokio::test]
async fn rusty_cloud_bed_order() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("plink_sim_10s_100v_10pmiss.bed")?;

    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    let val = ReadOptions::builder()
        .c()
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;

    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == -13.142); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[tokio::test]
async fn bad_header_cloud() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_object_path("badfile.bed")?;
    let bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .skip_early_check()
        .build()
        .await?;

    println!("{:?}", bed_cloud.object_path());

    // Attempt to create a new BedCloud instance and handle the error
    let result = BedCloud::from_object_path(&object_path).await;
    assert_error_variant!(result, BedErrorPlus::BedError(BedError::IllFormed(_)));

    Ok(())
}

#[tokio::test]
async fn bad_header_cloud_url() -> Result<(), Box<BedErrorPlus>> {
    println!("start");
    let url = sample_url("badfile.bed")?;
    println!("{:?}", url);
    let bed_cloud = BedCloud::builder(&url, EMPTY_OPTIONS)?
        .skip_early_check()
        .build()
        .await?;

    println!("{:?}", bed_cloud.object_path());

    // Attempt to create a new BedCloud instance and handle the error
    let result = BedCloud::new(&url, EMPTY_OPTIONS).await;
    println!("{:?}", result);
    assert_error_variant!(result, BedErrorPlus::BedError(BedError::IllFormed(_)));
    println!("done");

    Ok(())
}

#[tokio::test]
async fn doc_test_test_cloud() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("small.bed")?;

    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    let val = bed_cloud.read::<f64>().await?;
    assert_eq_nan(
        &val,
        &nd::array![
            [1.0, 0.0, f64::NAN, 0.0],
            [2.0, 0.0, f64::NAN, 2.0],
            [0.0, 1.0, 2.0, 0.0]
        ],
    );

    let object_path2 = sample_bed_object_path("some_missing.bed")?;
    let mut bed_cloud2 = BedCloud::from_object_path(&object_path2).await?;
    let val2 = ReadOptions::builder()
        .f64()
        .iid_index(s![..;2])
        .sid_index(20..30)
        .read_cloud(&mut bed_cloud2)
        .await?;
    assert!(val2.dim() == (50, 10));

    let mut bed_cloud3 = BedCloud::from_object_path(&object_path2).await?;
    println!("{:?}", bed_cloud3.iid().await?.slice(s![..5]));
    println!("{:?}", bed_cloud3.sid().await?.slice(s![..5]));
    println!(
        "{:?}",
        bed_cloud3
            .chromosome()
            .await?
            .iter()
            .collect::<HashSet<_>>()
    );
    let val3 = ReadOptions::builder()
        .sid_index(bed_cloud3.chromosome().await?.map(|elem| elem == "5"))
        .f64()
        .read_cloud(&mut bed_cloud3)
        .await?;
    assert!(val3.dim() == (100, 6));

    Ok(())
}

#[tokio::test]
async fn open_examples_cloud() -> Result<(), Box<BedErrorPlus>> {
    //     >>> from bed_reader import open_bed, sample_bed_file
    //     >>>
    //     >>> file_name = sample_bed_file("small.bed")
    //     >>> bed_cloud = open_bed(file_name)
    //     >>> print(bed_cloud.iid)
    //     ['iid1' 'iid2' 'iid3']
    //     >>> print(bed_cloud.sid)
    //     ['sid1' 'sid2' 'sid3' 'sid4']
    //     >>> print(bed_cloud.read())
    //     [[ 1.  0. nan  0.]
    //      [ 2.  0. nan  2.]
    //      [ 0.  1.  2.  0.]]
    //     >>> del bed_cloud  # optional: delete bed_cloud object

    let object_path = sample_bed_object_path("small.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;

    println!("{:?}", bed_cloud.iid().await?);
    println!("{:?}", bed_cloud.sid().await?);
    println!("{:?}", bed_cloud.read::<f64>().await?);

    // ["iid1", "iid2", "iid3"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1
    // ["sid1", "sid2", "sid3", "sid4"], shape=[4], strides=[1], layout=CFcf (0xf), const ndim=1
    // [[1.0, 0.0, NaN, 0.0],
    //  [2.0, 0.0, NaN, 2.0],
    //  [0.0, 1.0, 2.0, 0.0]], shape=[3, 4], strides=[1, 3], layout=Ff (0xa), const ndim=2

    // Open the file and read data for one SNP (variant)
    // at index position 2.

    // .. doctest::

    //     >>> import numpy as np
    //     >>> with open_bed(file_name) as bed_cloud:
    //     ...     print(bed_cloud.read(np.s_[:,2]))
    //     [[nan]
    //      [nan]
    //      [ 2.]]

    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    println!(
        "{:?}",
        ReadOptions::builder()
            .sid_index(2)
            .f64()
            .read_cloud(&mut bed_cloud)
            .await?
    );

    // [[NaN],
    //  [NaN],
    //  [2.0]], shape=[3, 1], strides=[1, 3], layout=CFcf (0xf), const ndim=2

    // Replace :attr:`iid`.

    //     >>> bed_cloud = open_bed(file_name, properties={"iid":["sample1","sample2","sample3"]})
    //     >>> print(bed_cloud.iid) # replaced
    //     ['sample1' 'sample2' 'sample3']
    //     >>> print(bed_cloud.sid) # same as before
    //     ['sid1' 'sid2' 'sid3' 'sid4']

    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .iid(["sample1", "sample2", "sample3"])
        .build()
        .await?;
    println!("{:?}", bed_cloud.iid().await?);
    println!("{:?}", bed_cloud.sid().await?);

    // ["sample1", "sample2", "sample3"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1
    // ["sid1", "sid2", "sid3", "sid4"], shape=[4], strides=[1], layout=CFcf (0xf), const ndim=

    // Do more testing of Rust
    let iid = nd::array!["sample1", "sample2", "sample3"];
    let mut _bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .iid(iid)
        .build()
        .await?;
    let iid = nd::array![
        "sample1".to_string(),
        "sample2".to_string(),
        "sample3".to_string()
    ];
    let mut _bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .iid(iid)
        .build()
        .await?;
    let iid = vec!["sample1", "sample2", "sample3"];
    let mut _bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .iid(iid)
        .build()
        .await?;
    let iid = vec![
        "sample1".to_string(),
        "sample2".to_string(),
        "sample3".to_string(),
    ];
    let mut _bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .iid(iid)
        .build()
        .await?;

    // Give the number of individuals (samples) and SNPs (variants) so that the .fam and
    // .bim files need never be opened.

    //     >>> with open_bed(file_name, iid_count=3, sid_count=4) as bed_cloud:
    //     ...     print(bed_cloud.read())
    //     [[ 1.  0. nan  0.]
    //      [ 2.  0. nan  2.]
    //      [ 0.  1.  2.  0.]]

    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .iid_count(3)
        .sid_count(4)
        .skip_early_check()
        .build()
        .await?;
    println!("{:?}", bed_cloud.read::<f64>().await?);

    //  [[1.0, 0.0, NaN, 0.0],
    //   [2.0, 0.0, NaN, 2.0],
    //   [0.0, 1.0, 2.0, 0.0]], shape=[3, 4], strides=[1, 3], layout=Ff (0xa), const ndim=2

    // Mark some properties as "don’t read or offer".

    //     >>> bed_cloud = open_bed(file_name, properties={
    //     ...    "father" : None, "mother" : None, "sex" : None, "pheno" : None,
    //     ...    "allele_1" : None, "allele_2":None })
    //     >>> print(bed_cloud.iid)        # read from file
    //     ['iid1' 'iid2' 'iid3']
    //     >>> print(bed_cloud.allele_2)   # not read and not offered
    //     None

    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .skip_allele_2()
        .build()
        .await?;
    println!("{:?}", bed_cloud.iid().await?);

    let result = bed_cloud.allele_2().await;
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::CannotUseSkippedMetadata(_))
    );

    Ok(())
}

#[tokio::test]
async fn metadata_etc_cloud() -> Result<(), Box<BedErrorPlus>> {
    // Initialize BedCloud with the sample file
    let url = sample_bed_url("small.bed")?;

    // Reading sex data
    let mut bed_cloud = BedCloud::new(&url, EMPTY_OPTIONS).await?;
    println!("{:?}", bed_cloud.sex().await?);
    // Expected output: [1, 2, 0], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1

    // Reading centiMorgan position
    let mut bed_cloud = BedCloud::new(&url, EMPTY_OPTIONS).await?;
    println!("{:?}", bed_cloud.cm_position().await?);
    // Expected output: [100.4, 2000.5, 4000.7, 7000.9], shape=[4], strides=[1], layout=CFcf (0xf), const ndim=1

    // Reading base pair position
    println!("{:?}", bed_cloud.bp_position().await?);
    // Expected output: [1, 100, 1000, 1004], shape=[4], strides=[1], layout=CFcf (0xf), const ndim=1

    // Reading family ID
    let mut bed_cloud = BedCloud::new(&url, EMPTY_OPTIONS).await?;
    println!("{:?}", bed_cloud.fid().await?);
    // Expected output: ["fid1", "fid1", "fid2"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1

    // Reading father ID
    println!("{:?}", bed_cloud.father().await?);
    // Expected output: ["iid23", "iid23", "iid22"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1

    // Reading mother ID
    println!("{:?}", bed_cloud.mother().await?);
    // Expected output: ["iid34", "iid34", "iid33"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1

    Ok(())
}

#[tokio::test]
async fn hello_father_cloud() -> Result<(), Box<BedErrorPlus>> {
    // Initialize BedCloud with the sample file and custom father metadata
    let object_path = sample_bed_object_path("small.bed")?;

    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .father(["f1", "f2", "f3"])
        .skip_mother()
        .build()
        .await?;

    println!("{:?}", bed_cloud.father().await?);
    // Expected output: ["f1", "f2", "f3"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1

    // Attempt to access mother data, expecting an error
    bed_cloud.mother().await.unwrap_err();

    Ok(())
}

#[tokio::test]
async fn max_concurrent_requests() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("plink_sim_10s_100v_10pmiss.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;

    // Read data with specified number of threads (or equivalent parallel processing setting)
    let val = ReadOptions::builder()
        .max_concurrent_requests(1)
        .max_chunk_size(1_000_000)
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == -13.142); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[tokio::test]
async fn fam_and_bim_cloud() -> Result<(), Box<BedErrorPlus>> {
    let deb_maf_mib = sample_object_paths(["small.deb", "small.maf", "small.mib"])?;

    // Build BedCloud with custom fam and bim paths
    let mut bed_cloud = BedCloud::builder_from_object_path(&deb_maf_mib[0])
        .fam_object_path(&deb_maf_mib[1])
        .bim_object_path(&deb_maf_mib[2])
        .build()
        .await?;

    // Read and process data
    println!("{:?}", bed_cloud.iid().await?);
    println!("{:?}", bed_cloud.sid().await?);
    let val: nd::Array2<i8> = bed_cloud.read().await?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == -20.5); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[tokio::test]
async fn fam_and_bim_cloud_url() -> Result<(), Box<BedErrorPlus>> {
    let mut deb_maf_mib = sample_urls(["small.deb", "small.maf", "small.mib"])?;

    // Build BedCloud with custom fam and bim paths
    let mut bed_cloud = BedCloud::builder(deb_maf_mib.remove(0), EMPTY_OPTIONS)?
        .fam_from_url(deb_maf_mib.remove(0), EMPTY_OPTIONS)? // Note: indexes shift
        .bim_from_url(deb_maf_mib.remove(0), EMPTY_OPTIONS)? // Note: indexes shift
        .build()
        .await?;

    // Read and process data
    println!("{:?}", bed_cloud.iid().await?);
    println!("{:?}", bed_cloud.sid().await?);
    let val: nd::Array2<i8> = bed_cloud.read().await?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == -20.5); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[tokio::test]
async fn readme_examples_cloud() -> Result<(), Box<BedErrorPlus>> {
    // Read genomic data from a .bed file.

    let object_path = sample_bed_object_path("small.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    let val = bed_cloud.read::<f64>().await?;
    println!("{val:?}");
    // [[1.0, 0.0, NaN, 0.0],
    // [2.0, 0.0, NaN, 2.0],
    // [0.0, 1.0, 2.0, 0.0]], shape=[3, 4], strides=[1, 3], layout=Ff (0xa), const ndim=2

    // Read every second individual and SNPs (variants) from 20 to 30.

    let object_path2 = sample_bed_object_path("some_missing.bed")?;
    let mut bed_cloud2 = BedCloud::from_object_path(&object_path2).await?;
    let val2 = ReadOptions::<f64>::builder()
        .iid_index(s![..;2])
        .sid_index(20..30)
        .read_cloud(&mut bed_cloud2)
        .await?;
    println!("{:?}", val2.dim()); // (50, 10)

    // List the first 5 individual (sample) ids, the first 5 SNP (variant) ids,
    // and every unique chromosome. Then, read every genomic value in chromosome 5.

    let mut bed_cloud3 = BedCloud::from_object_path(&object_path2).await?;
    let iid = bed_cloud3.iid().await?;
    let s = iid.slice(s![..5]);
    println!("{:?}", s);
    println!("{:?}", bed_cloud3.iid().await?.slice(s![..5]));
    println!("{:?}", bed_cloud3.sid().await?.slice(s![..5]));
    let unique = bed_cloud3
        .chromosome()
        .await?
        .iter()
        .collect::<HashSet<_>>();
    println!("{unique:?}");

    let is_5 = nd::Zip::from(bed_cloud3.chromosome().await?).par_map_collect(|elem| elem == "5");
    let val3 = ReadOptions::builder()
        .sid_index(is_5)
        .f64()
        .read_cloud(&mut bed_cloud3)
        .await?;
    println!("{:?}", val3.dim());
    // ["iid_0", "iid_1", "iid_2", "iid_3", "iid_4"], shape=[5], strides=[1], layout=CFcf (0xf), const ndim=1
    // ["sid_0", "sid_1", "sid_2", "sid_3", "sid_4"], shape=[5], strides=[1], layout=CFcf (0xf), const ndim=1
    // {"10", "11", "4", "21", "22", "14", "3", "12", "20", "15", "19", "8", "6", "18", "9", "2", "16", "13", "17", "1", "7", "5"}
    // (100, 6)
    Ok(())
}

#[tokio::test]
async fn range_cloud() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("small.bed")?;

    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    ReadOptions::builder()
        .iid_index(0..2)
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;
    ReadOptions::builder()
        .iid_index(0..=2)
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;
    ReadOptions::builder()
        .iid_index(..2)
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;
    ReadOptions::builder()
        .iid_index(..=2)
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;
    ReadOptions::builder()
        .iid_index(0..)
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;
    ReadOptions::builder()
        .iid_index(..)
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;

    Ok(())
}

#[tokio::test]
async fn nd_slice_cloud() -> Result<(), Box<BedErrorPlus>> {
    // ndarray operations (remain synchronous)
    let ndarray = nd::array![0, 1, 2, 3];
    println!("{:?}", ndarray.slice(nd::s![1..3])); // [1, 2]
    println!("{:?}", ndarray.slice(nd::s![1..3;-1])); // [2, 1]
    #[allow(clippy::reversed_empty_ranges)]
    let slice = nd::s![3..1;-1];
    println!("{:?}", ndarray.slice(slice)); // []

    // Reading BED file with various slice options
    let object_path = sample_bed_object_path("small.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;

    ReadOptions::builder()
        .iid_index(nd::s![0..2])
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;
    ReadOptions::builder()
        .iid_index(nd::s![..2])
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;
    ReadOptions::builder()
        .iid_index(nd::s![0..])
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;
    ReadOptions::builder()
        .iid_index(nd::s![0..2])
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;
    ReadOptions::builder()
        .iid_index(nd::s![-2..-1;-1])
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;

    Ok(())
}

#[tokio::test]
async fn skip_coverage_cloud() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("small.bed")?;

    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .skip_fid()
        .skip_iid()
        .skip_father()
        .skip_mother()
        .skip_sex()
        .skip_pheno()
        .skip_chromosome()
        .skip_sid()
        .skip_cm_position()
        .skip_bp_position()
        .skip_allele_1()
        .skip_allele_2()
        .build()
        .await?;

    // If the mother information is skipped, it should return an error.
    assert!(
        bed_cloud.mother().await.is_err(),
        "Mother data should not be available"
    );

    Ok(())
}

#[tokio::test]
async fn into_iter_cloud() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("small.bed")?;

    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .fid(["sample1", "sample2", "sample3"])
        .iid(["sample1", "sample2", "sample3"])
        .father(["sample1", "sample2", "sample3"])
        .mother(["sample1", "sample2", "sample3"])
        .sex([0, 0, 0])
        .pheno(["sample1", "sample2", "sample3"])
        .chromosome(["a", "b", "c", "d"])
        .sid(["a", "b", "c", "d"])
        .bp_position([0, 0, 0, 0])
        .cm_position([0.0, 0.0, 0.0, 0.0])
        .allele_1(["a", "b", "c", "d"])
        .allele_2(["a", "b", "c", "d"])
        .build()
        .await?;

    // Assuming pheno is an async method in BedCloud.
    let _ = bed_cloud.pheno().await?;
    Ok(())
}

fn rt1_cloud<R>(
    range_thing: R,
) -> Result<Result<nd::Array2<i8>, Box<BedErrorPlus>>, Box<BedErrorPlus>>
where
    R: std::ops::RangeBounds<usize>
        + std::fmt::Debug
        + Clone
        + std::slice::SliceIndex<[isize], Output = [isize]>
        + std::panic::RefUnwindSafe,
{
    println!("Running {:?}", &range_thing);

    let result1 = catch_unwind(|| {
        let rt = runtime::Runtime::new()?;

        rt.block_on(async {
            let object_path = sample_bed_object_path("toydata.5chrom.bed")?;
            let mut bed_cloud = BedCloud::from_object_path(&object_path).await.unwrap();
            let all: Vec<isize> = (0..(bed_cloud.iid_count().await.unwrap() as isize)).collect();
            let iid_index: &[isize] = &all[range_thing.clone()];

            ReadOptions::builder()
                .iid_index(iid_index)
                .i8()
                .read_cloud(&mut bed_cloud)
                .await
        })
    });
    match result1 {
        Err(_) => Err(BedError::PanickedThread().into()),
        Ok(bed_result) => Ok(bed_result),
    }
}

fn rt2_cloud(
    range_thing: bed_reader::Index,
) -> Result<Result<nd::Array2<i8>, Box<BedErrorPlus>>, Box<BedErrorPlus>> {
    println!("Running {:?}", &range_thing);
    let object_path = sample_bed_object_path("toydata.5chrom.bed")?;

    let result1 = catch_unwind(|| {
        let rt = runtime::Runtime::new()?;

        rt.block_on(async {
            let mut bed_cloud = BedCloud::from_object_path(&object_path).await.unwrap();
            ReadOptions::builder()
                .iid_index(range_thing)
                .i8()
                .read_cloud(&mut bed_cloud)
                .await
        })
    });
    match result1 {
        Err(_) => Err(BedError::PanickedThread().into()),
        Ok(bed_result) => Ok(bed_result),
    }
}

type RrArray2 = Result<Result<nd::Array2<i8>, Box<BedErrorPlus>>, Box<BedErrorPlus>>;
type RrUsize = Result<Result<usize, Box<BedErrorPlus>>, Box<BedErrorPlus>>;

fn rt23_cloud(range_thing: bed_reader::Index) -> (RrArray2, RrUsize) {
    (rt2_cloud(range_thing.clone()), rt3_cloud(range_thing))
}

fn rt3_cloud(
    range_thing: bed_reader::Index,
) -> Result<Result<usize, Box<BedErrorPlus>>, Box<BedErrorPlus>> {
    println!("Running {:?}", &range_thing);
    let object_path = sample_bed_object_path("toydata.5chrom.bed")?;

    let result1 = catch_unwind(|| {
        let rt = runtime::Runtime::new()?;

        rt.block_on(async {
            let mut bed_cloud = BedCloud::from_object_path(&object_path).await.unwrap();
            let count = bed_cloud.iid_count().await.unwrap();
            range_thing.len(count)
        })
    });
    match result1 {
        Err(_) => Err(BedError::PanickedThread().into()),
        Ok(bed_result) => Ok(bed_result),
    }
}

fn nds1_cloud(
    range_thing: SliceInfo1,
) -> Result<Result<nd::Array2<i8>, Box<BedErrorPlus>>, Box<BedErrorPlus>> {
    println!("Running {:?}", &range_thing);
    let object_path = sample_bed_object_path("toydata.5chrom.bed")?;

    let result1 = catch_unwind(|| {
        let rt = runtime::Runtime::new()?;

        rt.block_on(async {
            let mut bed_cloud = BedCloud::from_object_path(&object_path).await.unwrap();
            let all: nd::Array1<isize> =
                (0..(bed_cloud.iid_count().await.unwrap() as isize)).collect();
            let iid_index = &all.slice(&range_thing);

            ReadOptions::builder()
                .iid_index(iid_index)
                .i8()
                .read_cloud(&mut bed_cloud)
                .await
        })
    });
    match result1 {
        Err(_) => Err(BedError::PanickedThread().into()),
        Ok(bed_result) => Ok(bed_result),
    }
}

fn is_err2<T>(result_result: &Result<Result<T, Box<BedErrorPlus>>, Box<BedErrorPlus>>) -> bool {
    !matches!(result_result, Ok(Ok(_)))
}

fn assert_same_result(result1: RrArray2, result23: (RrArray2, RrUsize)) {
    let result2 = result23.0;
    let result3 = result23.1;
    let err1 = is_err2(&result1);
    let err2 = is_err2(&result2);
    let err3 = is_err2(&result3);

    if err1 || err2 || err3 {
        if !err1 || !err2 || !err3 {
            println!("{result1:?}");
            println!("{result2:?}");
            println!("{result3:?}");
            panic!("all should panic/error the same");
        }
        return;
    }

    let result1 = result1.unwrap().unwrap();
    let result2 = result2.unwrap().unwrap();
    let result3 = result3.unwrap().unwrap();
    println!("{result1:?}");
    println!("{result2:?}");
    println!("{result3:?}");
    assert!(
        allclose(&result1.view(), &result2.view(), 0, true),
        "not close"
    );
    assert!(result1.dim().0 == result3, "not same length");
}

#[test]
#[allow(clippy::reversed_empty_ranges)]
fn range_same_cloud() -> Result<(), Box<BedErrorPlus>> {
    let a = rt1_cloud(3..0);
    println!("{:?}", a);
    let b = rt23_cloud((3..0).into());
    println!("{:?}", b);

    assert_same_result(rt1_cloud(3..0), rt23_cloud((3..0).into()));
    assert_same_result(rt1_cloud(1000..), rt23_cloud((1000..).into()));

    assert_same_result(rt1_cloud(..), rt23_cloud((..).into()));
    assert_same_result(rt1_cloud(..3), rt23_cloud((..3).into()));
    assert_same_result(rt1_cloud(..=3), rt23_cloud((..=3).into()));
    assert_same_result(rt1_cloud(1..), rt23_cloud((1..).into()));
    assert_same_result(rt1_cloud(1..3), rt23_cloud((1..3).into()));
    assert_same_result(rt1_cloud(1..=3), rt23_cloud((1..=3).into()));
    assert_same_result(rt1_cloud(2..=2), rt23_cloud((2..=2).into()));
    Ok(())
}

#[test]
fn nd_slice_same_cloud() -> Result<(), Box<BedErrorPlus>> {
    assert_same_result(nds1_cloud(s![1000..]), rt23_cloud(s![1000..].into()));
    assert_same_result(nds1_cloud(s![..1000]), rt23_cloud(s![..1000].into()));
    assert_same_result(nds1_cloud(s![999..1000]), rt23_cloud(s![999..1000].into()));
    assert_same_result(nds1_cloud(s![-1000..]), rt23_cloud(s![-1000..].into()));
    assert_same_result(nds1_cloud(s![..-1000]), rt23_cloud(s![..-1000].into()));

    #[allow(clippy::reversed_empty_ranges)]
    assert_same_result(
        nds1_cloud(s![-999..-1000]),
        rt23_cloud(s![-999..-1000].into()),
    );
    #[allow(clippy::reversed_empty_ranges)]
    assert_same_result(nds1_cloud(s![3..0]), rt23_cloud(s![3..0].into()));
    #[allow(clippy::reversed_empty_ranges)]
    assert_same_result(nds1_cloud(s![-1..-2]), rt23_cloud(s![-1..-2].into()));

    assert_same_result(nds1_cloud(s![..-3]), rt23_cloud(s![..-3].into()));
    assert_same_result(nds1_cloud(s![..=-3]), rt23_cloud(s![..=-3].into()));
    assert_same_result(nds1_cloud(s![-1..]), rt23_cloud(s![-1..].into()));
    assert_same_result(nds1_cloud(s![-3..-1]), rt23_cloud(s![-3..-1].into()));
    assert_same_result(nds1_cloud(s![-3..=-1]), rt23_cloud(s![-3..=-1].into()));
    assert_same_result(nds1_cloud(s![-2..=-2]), rt23_cloud(s![-2..=-2].into()));

    #[allow(clippy::reversed_empty_ranges)]
    assert_same_result(nds1_cloud(s![1..-1]), rt23_cloud(s![1..-1].into()));

    assert_same_result(nds1_cloud(s![..]), rt23_cloud((s![..]).into()));
    assert_same_result(nds1_cloud(s![..3]), rt23_cloud((s![..3]).into()));
    assert_same_result(nds1_cloud(s![..=3]), rt23_cloud((s![..=3]).into()));
    assert_same_result(nds1_cloud(s![1..]), rt23_cloud((s![1..]).into()));
    assert_same_result(nds1_cloud(s![1..3]), rt23_cloud((s![1..3]).into()));
    assert_same_result(nds1_cloud(s![1..=3]), rt23_cloud((s![1..=3]).into()));
    assert_same_result(nds1_cloud(s![2..=2]), rt23_cloud(s![2..=2].into()));

    Ok(())
}

#[tokio::test]
async fn bool_read_cloud() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("small.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;

    let result = ReadOptions::builder()
        .iid_index([false, false, true, false])
        .i8()
        .read_cloud(&mut bed_cloud)
        .await;
    println!("{result:?}");
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::BoolArrayVectorWrongLength(_, _))
    );

    let _val = ReadOptions::builder()
        .iid_index([false, false, true])
        .i8()
        .read_cloud(&mut bed_cloud)
        .await?;

    Ok(())
}
#[tokio::test]
async fn i8_etc_cloud() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("small.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;

    let _val = ReadOptions::builder()
        .f()
        .i8()
        .iid_index([false, false, true])
        .read_cloud(&mut bed_cloud)
        .await?;

    Ok(())
}

#[tokio::test]
async fn fill_cloud() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("small.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;

    let read_options = ReadOptions::builder()
        .f()
        .i8()
        .iid_index([false, false, true])
        .build()?;

    let mut val = nd::Array2::<i8>::default((3, 4));
    let result = bed_cloud
        .read_and_fill_with_options(&mut val.view_mut(), &read_options)
        .await;
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InvalidShape(_, _, _, _))
    );

    let mut val = nd::Array2::<i8>::default((1, 4));
    bed_cloud
        .read_and_fill_with_options(&mut val.view_mut(), &read_options)
        .await?;

    assert_eq!(bed_cloud.dim().await?, (3, 4));

    Ok(())
}

#[tokio::test]
async fn read_options_builder_cloud() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("small.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;

    // Read the SNPs indexed by 2.
    let val = ReadOptions::builder()
        .sid_index(2)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);

    // Read the SNPs indexed by 2, 3, and 0.
    let val = ReadOptions::builder()
        .sid_index([2, 3, 0])
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert_eq_nan(
        &val,
        &nd::array![[f64::NAN, 0.0, 1.0], [f64::NAN, 2.0, 2.0], [2.0, 0.0, 0.0]],
    );

    // Read SNPs from 1 (inclusive) to 4 (exclusive).
    let val = ReadOptions::builder()
        .sid_index(1..4)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert_eq_nan(
        &val,
        &nd::array![[0.0, f64::NAN, 0.0], [0.0, f64::NAN, 2.0], [1.0, 2.0, 0.0]],
    );

    // Print unique chrom values. Then, read all SNPs in chrom 5.
    use std::collections::HashSet;

    println!(
        "{:?}",
        bed_cloud.chromosome().await?.iter().collect::<HashSet<_>>()
    );
    let val = ReadOptions::builder()
        .sid_index(bed_cloud.chromosome().await?.map(|elem| elem == "5"))
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);

    // Read 1st individual (across all SNPs).
    let val = ReadOptions::builder()
        .iid_index(0)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert_eq_nan(&val, &nd::array![[1.0, 0.0, f64::NAN, 0.0]]);

    // Read every 2nd individual.
    use ndarray::s;
    let val = ReadOptions::builder()
        .iid_index(s![..;2])
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert_eq_nan(
        &val,
        &nd::array![[1.0, 0.0, f64::NAN, 0.0], [0.0, 1.0, 2.0, 0.0]],
    );

    // Read last and 2nd-to-last individuals and the last SNPs
    let val = ReadOptions::builder()
        .iid_index([-1, -2])
        .sid_index(-1)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    println!("{:?}", &val);
    assert_eq_nan(&val, &nd::array![[0.0], [2.0]]);

    Ok(())
}

#[tokio::test]
async fn bed_builder_cloud() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("small.bed")?;

    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .build()
        .await?;
    println!("{:?}", bed_cloud.iid().await?);
    println!("{:?}", bed_cloud.sid().await?);
    let val = bed_cloud.read::<f64>().await?;

    assert_eq_nan(
        &val,
        &nd::array![
            [1.0, 0.0, f64::NAN, 0.0],
            [2.0, 0.0, f64::NAN, 2.0],
            [0.0, 1.0, 2.0, 0.0]
        ],
    );

    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .build()
        .await?;
    let val = ReadOptions::builder()
        .sid_index(2)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;

    assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);

    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .iid(["sample1", "sample2", "sample3"])
        .build()
        .await?;
    println!("{:?}", bed_cloud.iid().await?); // replaced
    println!("{:?}", bed_cloud.sid().await?); // same as before

    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .iid_count(3)
        .sid_count(4)
        .build()
        .await?;
    let val = bed_cloud.read::<f64>().await?;
    assert_eq_nan(
        &val,
        &nd::array![
            [1.0, 0.0, f64::NAN, 0.0],
            [2.0, 0.0, f64::NAN, 2.0],
            [0.0, 1.0, 2.0, 0.0]
        ],
    );

    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .skip_father()
        .skip_mother()
        .skip_sex()
        .skip_pheno()
        .skip_allele_1()
        .skip_allele_2()
        .build()
        .await?;
    println!("{:?}", bed_cloud.iid().await?);
    bed_cloud.allele_2().await.expect_err("Can't be read");

    Ok(())
}

#[tokio::test]
#[allow(clippy::needless_borrows_for_generic_args)]
async fn negative_indexing_cloud() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("small.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;

    // println!("{:?}", bed_cloud.read::<f64>().await?);
    // [[1.0, 0.0, NaN, 0.0],
    // [2.0, 0.0, NaN, 2.0],
    // [0.0, 1.0, 2.0, 0.0]], shape=[3, 4], strides=[1, 3], layout=Ff (0xa), const ndim=2
    //  iid range is -4ERROR -3 -2 -1 0 1 2 3ERROR
    //  sid range is -5ERROR -4 ... 3 4ERROR
    for index in [-4, 3] {
        match ReadOptions::builder()
            .iid_index(index)
            .i8()
            .read_cloud(&mut bed_cloud)
            .await
        {
            Err(ref boxed_error) => match **boxed_error {
                BedErrorPlus::BedError(BedError::IidIndexTooBig(x)) => {
                    assert_eq!(x, index);
                }
                _ => panic!("test failure"),
            },
            _ => panic!("test failure"),
        }
    }

    for index in [-3, 0] {
        let val = ReadOptions::builder()
            .iid_index(index)
            .i8()
            .read_cloud(&mut bed_cloud)
            .await?;
        // println!("{val:?}");
        assert!(val[[0, 0]] == 1);
    }

    for index in [-1, 2] {
        let val = ReadOptions::builder()
            .iid_index(index)
            .i8()
            .read_cloud(&mut bed_cloud)
            .await?;
        // println!("{val:?}");
        assert!(val[[0, 0]] == 0);
    }

    for index in [-5, 4] {
        match ReadOptions::builder()
            .sid_index(index)
            .i8()
            .read_cloud(&mut bed_cloud)
            .await
        {
            Err(ref boxed_error) => match **boxed_error {
                BedErrorPlus::BedError(BedError::SidIndexTooBig(x)) => {
                    assert_eq!(x, index);
                }
                _ => panic!("test failure"),
            },
            _ => panic!("test failure"),
        }
    }

    for index in [-4, 0] {
        let val = ReadOptions::builder()
            .sid_index(index)
            .i8()
            .read_cloud(&mut bed_cloud)
            .await?;
        // println!("{val:?}");
        assert!(val[[0, 0]] == 1);
    }

    for index in [-1, 3] {
        let val = ReadOptions::builder()
            .sid_index(index)
            .i8()
            .read_cloud(&mut bed_cloud)
            .await?;
        // println!("{val:?}");
        assert!(val[[0, 0]] == 0);
    }

    Ok(())
}

#[tokio::test]
async fn index_doc_cloud() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("some_missing.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;

    println!("{:?}", bed_cloud.dim().await?); // prints (100, 100)

    // Read all individuals and all SNPs
    let val = ReadOptions::builder()
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert!(val.dim() == (100, 100));

    // Read the individual at index position 10 and all SNPs
    let val = ReadOptions::builder()
        .iid_index(10)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert!(val.dim() == (1, 100));

    // Read the individuals at index positions 0, 5, 1st-from-the-end and the SNP at index position 3
    let val = ReadOptions::builder()
        .iid_index(vec![0, 5, -1])
        .sid_index(3)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert!(val.dim() == (3, 1));

    // Repeat, but with an ndarray
    let val = ReadOptions::builder()
        .iid_index(nd::array![0, 5, -1])
        .sid_index(3)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert!(val.dim() == (3, 1));

    // Repeat, but with a Rust array
    let val = ReadOptions::builder()
        .iid_index([0, 5, -1])
        .sid_index(3)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert!(val.dim() == (3, 1));

    // Create a boolean ndarray identifying SNPs in chromosome 5, then select those SNPs.
    let chrom_5 = bed_cloud.chromosome().await?.map(|elem| elem == "5");
    let val = ReadOptions::builder()
        .sid_index(chrom_5)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert!(val.dim() == (100, 6));

    // Use ndarray's slice macro, `s!`, to select every 2nd individual and every 3rd SNP.
    let val = ReadOptions::builder()
        .iid_index(s![..;2])
        .sid_index(s![..;3])
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert!(val.dim() == (50, 34));

    // Use ndarray's slice macro, [`s!`](https://docs.rs/ndarray/latest/ndarray/macro.s.html),
    // to select the 10th-from-last individual to the last, in reverse order,
    // and every 3rd SNP in reverse order.)
    let val = ReadOptions::builder()
        .iid_index(s![-10..;-1])
        .sid_index(s![..;-3])
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert!(val.dim() == (10, 34));
    Ok(())
}

#[tokio::test]
async fn cloud_index_options() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("some_missing.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;

    #[allow(clippy::let_unit_value)]
    let index: () = ();

    let all = ReadOptions::builder()
        .iid_index(index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert!(all.dim() == (100, 100));

    let mut index: [bool; 100] = [false; 100];
    index[0] = true;
    index[2] = true;
    let val = ReadOptions::builder()
        .iid_index(index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected = all
        .select(nd::Axis(0), [0, 2].as_slice())
        .select(nd::Axis(1), [0, 2].as_slice());
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let mut index: nd::Array1<bool> = nd::Array::from_elem(100, false);
    index[0] = true;
    index[2] = true;
    let val = ReadOptions::builder()
        .iid_index(&index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected = all
        .select(nd::Axis(0), [0, 2].as_slice())
        .select(nd::Axis(1), [0, 2].as_slice());
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let mut index: Vec<bool> = vec![false; 100];
    index[0] = true;
    index[2] = true;
    let val = ReadOptions::builder()
        .iid_index(&index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected = all
        .select(nd::Axis(0), [0, 2].as_slice())
        .select(nd::Axis(1), [0, 2].as_slice());
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let index: isize = 2;
    let val = ReadOptions::builder()
        .iid_index(index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected = all.slice(s![2isize..=2, 2isize..=2]);
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let index: isize = -1;
    let val = ReadOptions::builder()
        .iid_index(index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected = all.slice(s![99isize..=99, 99isize..=99]);
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let index: Vec<isize> = vec![0, 10, -2];
    let val = ReadOptions::builder()
        .iid_index(&index)
        .sid_index(&index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected_index = vec![0, 10, 98usize];
    let expected = all
        .select(nd::Axis(0), expected_index.as_slice())
        .select(nd::Axis(1), expected_index.as_slice());
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let index: &[isize] = &[0, 10, -2];
    let val = ReadOptions::builder()
        .iid_index(index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected_index = vec![0, 10, 98usize];
    let expected = all
        .select(nd::Axis(0), expected_index.as_slice())
        .select(nd::Axis(1), expected_index.as_slice());
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let val = ReadOptions::builder()
        .iid_index([0, 10, -2])
        .sid_index([0, 10, -2])
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected_index = vec![0, 10, 98usize];
    let expected = all
        .select(nd::Axis(0), expected_index.as_slice())
        .select(nd::Axis(1), expected_index.as_slice());
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let index: nd::Array1<isize> = nd::array![0, 10, -2];
    let val = ReadOptions::builder()
        .iid_index(&index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected_index = vec![0, 10, 98usize];
    let expected = all
        .select(nd::Axis(0), expected_index.as_slice())
        .select(nd::Axis(1), expected_index.as_slice());
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let index: std::ops::Range<usize> = 10..20;
    let val = ReadOptions::builder()
        .iid_index(&index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected = all.slice(s![10usize..20, 10usize..20]);
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let index: std::ops::RangeFrom<usize> = 50..;
    let val = ReadOptions::builder()
        .iid_index(&index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected = all.slice(s![50usize.., 50usize..]);
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let index: std::ops::RangeFull = ..;
    let val = ReadOptions::builder()
        .iid_index(index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected = all.slice(s![.., ..]);
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let index: std::ops::RangeTo<usize> = ..3;
    let val = ReadOptions::builder()
        .iid_index(index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected = all.slice(s![..3, ..3]);
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let index: std::ops::RangeToInclusive<usize> = ..=19;
    let val = ReadOptions::builder()
        .iid_index(index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected = all.slice(s![..=19, ..=19]);
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let index: std::ops::RangeInclusive<usize> = 1..=3;
    let val = ReadOptions::builder()
        .iid_index(&index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected = all.slice(s![1..=3, 1..=3]);
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    let index: SliceInfo1 = s![-20..-10;-2];
    let val = ReadOptions::builder()
        .iid_index(index)
        .sid_index(index)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    let expected = all.slice(s![-20..-10;-2,-20..-10;-2]);
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    Ok(())
}

#[tokio::test]
async fn cloud_set_metadata() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("small.bed")?;
    let metadata = Metadata::builder()
        .iid(["iid1", "iid2", "iid3"])
        .sid(["sid1", "sid2", "sid3", "sid4"])
        .build()?;
    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .metadata(&metadata)
        .build()
        .await?;
    let metadata2 = bed_cloud.metadata().await?;
    println!("{metadata2:?}");

    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    let metadata = bed_cloud.metadata().await?;
    println!("{metadata:?}");

    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .metadata(&metadata)
        .build()
        .await?;
    let metadata2 = bed_cloud.metadata().await?;
    println!("{metadata2:?}");

    Ok(())
}

#[tokio::test]
async fn cloud_metadata_print() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("small.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;

    let fid = bed_cloud.fid().await?;
    println!("{fid:?}"); // Outputs ndarray ["fid1", "fid1", "fid2"]
    let iid = bed_cloud.iid().await?;
    println!("{iid:?}"); // Outputs ndarray ["iid1", "iid2", "iid3"]
    let father = bed_cloud.father().await?;
    println!("{father:?}"); // Outputs ndarray ["iid23", "iid23", "iid22"]
    let mother = bed_cloud.mother().await?;
    println!("{mother:?}"); // Outputs ndarray ["iid34", "iid34", "iid33"]
    let sex = bed_cloud.sex().await?;
    println!("{sex:?}"); // Outputs ndarray [1, 2, 0]
    let pheno = bed_cloud.pheno().await?;
    println!("{pheno:?}"); // Outputs ndarray ["red", "red", "blue"]

    let chromosome = bed_cloud.chromosome().await?;
    println!("{chromosome:?}"); // Outputs ndarray ["1", "1", "5", "Y"
    let sid = bed_cloud.sid().await?;
    println!("{sid:?}"); // Outputs ndarray "sid1", "sid2", "sid3", "sid4"]
    let cm_position = bed_cloud.cm_position().await?;
    println!("{cm_position:?}"); // Outputs ndarray [100.4, 2000.5, 4000.7, 7000.9]
    let bp_position = bed_cloud.bp_position().await?;
    println!("{bp_position:?}"); // Outputs ndarray [1, 100, 1000, 1004]
    let allele_1 = bed_cloud.allele_1().await?;
    println!("{allele_1:?}"); // Outputs ndarray ["A", "T", "A", "T"]
    let allele_2 = bed_cloud.allele_2().await?;
    println!("{allele_2:?}"); // Outputs ndarray ["A", "C", "C", "G"]
    Ok(())
}

#[tokio::test]
async fn cloud_iid_index() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("some_missing.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;

    // Read the individual at index position 3

    let val = ReadOptions::builder()
        .iid_index(3)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert!(val.dim() == (1, 100));

    // Read the individuals at index positions 0, 5, and 1st-from-last.

    let val = ReadOptions::builder()
        .iid_index([0, 5, -1])
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;

    assert!(val.dim() == (3, 100));

    // Read the individuals at index positions 20 (inclusive) to 30 (exclusive).

    let val = ReadOptions::builder()
        .iid_index(20..30)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;

    assert!(val.dim() == (10, 100));

    // Read the individuals at every 2nd index position.

    let val = ReadOptions::builder()
        .iid_index(s![..;2])
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;

    assert!(val.dim() == (50, 100));

    // Read chromosome 5 of the the female individuals.

    let female = bed_cloud.sex().await?.map(|elem| *elem == 2);
    let chrom_5 = bed_cloud.chromosome().await?.map(|elem| elem == "5");
    let val = ReadOptions::builder()
        .iid_index(female)
        .sid_index(chrom_5)
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;

    println!("{:?}", val.dim());
    assert_eq!(val.dim(), (50, 6));

    Ok(())
}

#[tokio::test]
async fn cloud_struct_play() -> Result<(), Box<BedErrorPlus>> {
    // Bed
    // can't construct Bed directly because some fields are private
    // can't change pub properties because there are none

    // make ReadOptions directly or change? no, no pub fields
    // make WriteOptions directly or change? no, no pub fields
    // make Metadata directly or change? no, no pub fields

    // Can you change a value in a vector? No, because can't be borrowed as mutable
    let metadata = Metadata::builder().build()?.fill(100, 100)?;
    println!("{0:?}", metadata.iid());
    Ok(())
}

#[tokio::test]
async fn cloud_metadata_bed() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("small.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    let metadata = bed_cloud.metadata().await?;
    println!("{0:?}", metadata.iid()); // Outputs Some(["iid1", "iid2", "iid3"] ...)
    println!("{0:?}", metadata.sid()); // Outputs Some(["sid1", "sid2", "sid3", "sid4"] ...)
    Ok(())
}

#[tokio::test]
async fn cloud_read_and_fill_with_options() -> Result<(), Box<BedErrorPlus>> {
    use bed_reader::assert_eq_nan;
    use bed_reader::ReadOptions;
    use ndarray as nd;
    // Read the SNPs indexed by 2.
    let object_path = sample_bed_object_path("small.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;

    let read_options = ReadOptions::builder().sid_index(2).build()?;
    let mut val = nd::Array2::<f64>::default((3, 1));
    bed_cloud
        .read_and_fill_with_options(&mut val.view_mut(), &read_options)
        .await?;

    assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);

    let object_path = sample_bed_object_path("small.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    let mut val = nd::Array2::<f64>::default((3, 1));
    ReadOptions::builder()
        .sid_index(2)
        .read_and_fill_cloud(&mut bed_cloud, &mut val.view_mut())
        .await?;

    assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);

    Ok(())
}

#[tokio::test]
async fn cloud_bed_builder_metadata() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("small.bed")?;

    // show that can fine errors
    let metadata = Metadata::builder()
        .iid(["i1", "i2", "i3"])
        .sid(["s1", "s2", "s3", "s4"])
        .build()?;
    let result = BedCloud::builder_from_object_path(&object_path)
        .fid(["f1", "f2", "f3", "f4"])
        .metadata(&metadata)
        .build()
        .await;
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    let metadata = Metadata::builder()
        .iid(["i1", "i2", "i3"])
        .sid(["s1", "s2", "s3", "s4"])
        .build()?;
    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .fid(["f1", "f2", "f3"])
        .iid(["x1", "x2", "x3"])
        .metadata(&metadata)
        .build()
        .await?;
    println!("{0:?}", bed_cloud.fid().await?); // Outputs ndarray ["f1", "f2", "f3"]
    println!("{0:?}", bed_cloud.iid().await?); // Outputs ndarray ["i1", "i2", "i3"]
    println!("{0:?}", bed_cloud.sid().await?); // Outputs ndarray ["s1", "s2", "s3", "s4"]
    println!("{0:?}", bed_cloud.chromosome().await?); // Outputs ndarray ["1", "1", "5", "Y"]

    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .skip_fid()
        .fid(["f1", "f2", "f3"])
        .iid(["x1", "x2", "x3"])
        .metadata(&metadata)
        .build()
        .await?;
    bed_cloud.fid().await.expect_err("Should fail");
    println!("{0:?}", bed_cloud.iid().await?); // Outputs ndarray ["i1", "i2", "i3"]
    println!("{0:?}", bed_cloud.sid().await?); // Outputs ndarray ["s1", "s2", "s3", "s4"]
    println!("{0:?}", bed_cloud.chromosome().await?); // Outputs ndarray ["1", "1", "5", "Y"]

    Ok(())
}

#[tokio::test]
async fn cloud_metadata_read_fam_bim() -> Result<(), Box<BedErrorPlus>> {
    let skip_set = HashSet::<MetadataFields>::new();
    let metadata_empty = Metadata::new();
    let object_path = sample_object_path("small.fam")?;
    let (metadata_fam, iid_count) = metadata_empty
        .read_fam_cloud(&object_path, &skip_set)
        .await?;
    let object_path = sample_object_path("small.bim")?;
    let (metadata_bim, sid_count) = metadata_fam.read_bim_cloud(&object_path, &skip_set).await?;
    assert_eq!(iid_count, 3);
    assert_eq!(sid_count, 4);
    println!("{0:?}", metadata_bim.iid()); // Outputs optional ndarray Some(["iid1", "iid2", "iid3"]...)
    println!("{0:?}", metadata_bim.sid()); // Outputs optional ndarray Some(["sid1", "sid2", "sid3", "sid4"]...)
    println!("{0:?}", metadata_bim.chromosome()); // Outputs optional ndarray Some(["1", "1", "5", "Y"]...)

    Ok(())
}

#[tokio::test]
async fn cloud_read_options_properties() -> Result<(), Box<BedErrorPlus>> {
    let read_options = ReadOptions::builder().sid_index([2, 3, 0]).i8().build()?;
    assert_eq!(read_options.missing_value(), -127);
    println!("{0:?}", read_options.iid_index()); // Outputs 'All'
    println!("{0:?}", read_options.sid_index()); // Outputs 'Vec([2, 3, 0])'
    assert!(read_options.is_f());
    assert!(read_options.is_a1_counted());
    assert_eq!(read_options.num_threads(), None);

    let object_path = sample_bed_object_path("small.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    let val = bed_cloud.read_with_options(&read_options).await?;

    assert_eq_nan(&val, &nd::array![[-127, 0, 1], [-127, 2, 2], [2, 0, 0]]);

    Ok(())
}

// Bed, file set iid and fid and metadata and iid_count, and check if error is caught
#[tokio::test]
async fn cloud_bed_inconsistent_count() -> Result<(), Box<BedErrorPlus>> {
    // Bed: fid vs metadata
    let metadata = Metadata::builder().iid(["f1", "f2", "f3", "f4"]).build()?;
    let object_path = sample_bed_object_path("small.bed")?;
    let result = BedCloud::builder_from_object_path(&object_path)
        .fid(["f1", "f2", "f3"])
        .metadata(&metadata)
        .build()
        .await;
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    // Bed: file vs file:
    let bad_fam_object_path = sample_object_path("small.fam_bad")?;
    let mut bed_cloud = BedCloud::builder_from_object_path(&sample_bed_object_path("small.bed")?)
        .fam_object_path(&bad_fam_object_path)
        .skip_iid()
        .skip_father()
        .skip_mother()
        .skip_sex()
        .skip_pheno()
        .build()
        .await?;
    let result = bed_cloud.fid().await;
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::MetadataFieldCount(_, _, _))
    );

    // Bed: fid vs iid
    let result = BedCloud::builder_from_object_path(&object_path)
        .fid(["f1", "f2", "f3"])
        .iid(["i1", "i2", "i3", "i4"])
        .build()
        .await;
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    // Bed: iid vs file
    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .iid(["i1", "i2", "i3", "i4"])
        .build()
        .await?;
    let result = bed_cloud.fid().await;
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    Ok(())
}

// MetadataBuilders
#[tokio::test]
async fn cloud_metadata_inconsistent_count() -> Result<(), Box<BedErrorPlus>> {
    let skip_set = HashSet::<MetadataFields>::new();

    // Metadata: iid vs file
    let metadata = Metadata::builder().iid(["i1", "i2", "i3", "i4"]).build()?;
    let result = metadata
        .read_fam_cloud(&sample_object_path("small.fam")?, &skip_set)
        .await;
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    // Metadata: fid vs metadata
    let metadata = Metadata::builder().iid(["f1", "f2", "f3", "f4"]).build()?;
    let result = Metadata::builder()
        .fid(["f1", "f2", "f3"])
        .metadata(&metadata)
        .build();
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    // Metadata: file vs file:
    let metadata = Metadata::builder().build()?;
    let result = metadata
        .read_fam_cloud(&sample_object_path("small.fam_bad")?, &skip_set)
        .await;
    println!("{0:?}", result);
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::MetadataFieldCount(_, _, _))
    );

    // Metadata: fid vs iid
    let result = Metadata::builder()
        .fid(["f1", "f2", "f3"])
        .iid(["i1", "i2", "i3", "i4"])
        .build();
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    // Metadata: iid vs fill
    let metadata = Metadata::builder().iid(["i1", "i2", "i3", "i4"]).build()?;
    let result = metadata.fill(3, 4);
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    Ok(())
}

#[tokio::test]
async fn cloud_parsing_metadata() -> Result<(), Box<BedErrorPlus>> {
    let bed_fam_bim =
        sample_object_paths(["small.bed", "small.fam", "small.bim_bad_positions.bim"])?;
    let mut bed_cloud = BedCloud::builder_from_object_path(&bed_fam_bim[0])
        .bim_object_path(&bed_fam_bim[2])
        .build()
        .await?;
    let result = bed_cloud.cm_position().await;
    assert_error_variant!(result, BedErrorPlus::ParseIntError(_));
    Ok(())
}

#[tokio::test]
async fn cloud_read_fam() -> Result<(), Box<BedErrorPlus>> {
    let skip_set = HashSet::<MetadataFields>::new();
    let metadata_empty = Metadata::new();
    let fam_object_path = sample_object_path("small.fam")?;
    let (metadata_fam, _) = metadata_empty
        .read_fam_cloud(&fam_object_path, &skip_set)
        .await?;
    // metadata_empty.read_fam("bed_reader/tests/data/small.fam", &skip_set)?;
    println!("{:?}", metadata_fam.iid()); // Outputs optional ndarray Some(["iid1", "iid2", "iid3"]...)
    Ok(())
}

#[tokio::test]
async fn cloud_lib_intro() -> Result<(), Box<BedErrorPlus>> {
    let object_path = sample_bed_object_path("some_missing.bed")?;

    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    println!("{:?}", bed_cloud.iid().await?.slice(s![..5])); // Outputs ndarray: ["iid_0", "iid_1", "iid_2", "iid_3", "iid_4"]
    println!("{:?}", bed_cloud.sid().await?.slice(s![..5])); // Outputs ndarray: ["sid_0", "sid_1", "sid_2", "sid_3", "sid_4"]
    println!(
        "{:?}",
        bed_cloud.chromosome().await?.iter().collect::<HashSet<_>>()
    );
    // Outputs: {"12", "10", "4", "8", "19", "21", "9", "15", "6", "16", "13", "7", "17", "18", "1", "22", "11", "2", "20", "3", "5", "14"}
    let _ = ReadOptions::builder()
        .sid_index(bed_cloud.chromosome().await?.map(|elem| elem == "5"))
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;

    Ok(())
}

#[tokio::test]
async fn read_index_quick_cloud() -> Result<(), Box<BedErrorPlus>> {
    use bed_reader::Bed;

    let read_options = ReadOptions::builder()
        // .iid_index([-1, -2])
        .sid_index(-1)
        .f64()
        .build()?;

    let path = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(path)?;
    let val0 = bed.read_with_options(&read_options)?;
    println!("0: {:?}", &val0);

    let object_path = sample_bed_object_path("small.bed")?;
    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    let val1 = bed_cloud.read_with_options(&read_options).await?;
    println!("1: {:?}", &val1);

    assert_eq_nan(&val0, &val1);

    Ok(())
}

#[tokio::test]
async fn object_path_extension() -> Result<(), Box<BedErrorPlus>> {
    let mut object_path = sample_bed_object_path("plink_sim_10s_100v_10pmiss.bed")?;
    assert_eq!(object_path.size().await?, 303);
    object_path.set_extension("fam")?;
    assert_eq!(object_path.size().await?, 130);
    Ok(())
}

// The AWS tests are skipped to credentials are not available.
#[tokio::test]
async fn s3_cloud() -> Result<(), Box<BedErrorPlus>> {
    use rusoto_credential::{CredentialsError, ProfileProvider, ProvideAwsCredentials};
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

    let s3 = AmazonS3Builder::new()
        .with_region("us-west-2")
        .with_bucket_name("bedreader")
        .with_access_key_id(credentials.aws_access_key_id())
        .with_secret_access_key(credentials.aws_secret_access_key())
        .build()
        .unwrap();

    // LATER should we accept a string as a Path?
    let store_path = StorePath::parse("/v1/toydata.5chrom.bed").unwrap();
    let object_path = ObjectPath::new(Arc::new(s3), store_path);
    assert_eq!(object_path.size().await?, 1_250_003);
    Ok(())
}

#[tokio::test]
async fn dyn_cloud() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;
    let url = path_to_url_string(file_name)?;
    let url = Url::parse(&url).unwrap();
    let iid_count = 3;
    let sid_count = 4;
    let (object_store, store_path) = object_store::parse_url(&url).unwrap();
    let object_path: ObjectPath<Box<dyn ObjectStore>> =
        ObjectPath::new(Arc::new(object_store), store_path);

    let mut bed_cloud = BedCloud::builder_from_object_path(&object_path)
        .iid_count(iid_count)
        .sid_count(sid_count)
        .build()
        .await?;
    let val = bed_cloud.read::<i8>().await?;
    println!("{val:?}");

    assert_eq_nan(
        &val,
        &nd::array![[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]],
    );

    Ok(())
}

#[tokio::test]
async fn s3_url_cloud() -> Result<(), Box<BedErrorPlus>> {
    // Read my AWS credentials from file ~/.aws/credentials
    use rusoto_credential::{CredentialsError, ProfileProvider, ProvideAwsCredentials};
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
    Ok(())
}

#[tokio::test]
async fn s3_url_cloud2() -> Result<(), Box<BedErrorPlus>> {
    // Read my AWS credentials from file ~/.aws/credentials
    use rusoto_credential::{CredentialsError, ProfileProvider, ProvideAwsCredentials};
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
    let url = Url::parse(url).unwrap();
    let (object_store, store_path): (Box<dyn ObjectStore>, StorePath) =
        object_store::parse_url_opts(&url, options).unwrap();
    // print!("{:?}", object_store);
    // print!("{:?}", store_path);
    // // let store_path: StorePath = "/v1/toydata.5chrom.bed".into();
    let object_path: ObjectPath<Box<dyn ObjectStore>> =
        ObjectPath::new(Arc::new(object_store), store_path);
    // assert_eq!(object_path.size().await?, 1_250_003);

    let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
    let val = bed_cloud.read::<i8>().await?;
    assert_eq!(val.shape(), &[500, 10_000]);

    Ok(())
}

#[test]
fn object_path_2() -> Result<(), Box<BedErrorPlus>> {
    use bed_reader::{sample_bed_file, BedErrorPlus, ObjectPath};
    use object_store::{local::LocalFileSystem, path::Path as StorePath};
    use std::sync::Arc;
    use tokio::runtime::Runtime;

    Runtime::new().unwrap().block_on(async {
        let arc_object_store = Arc::new(LocalFileSystem::new()); // Arc-wrapped ObjectStore
        let file_path = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?; // regular Rust PathBuf
        let store_path = StorePath::from_filesystem_path(&file_path)?; // StorePath
        let object_path = ObjectPath::new(arc_object_store, store_path); // ObjectPath
        assert_eq!(object_path.size().await?, 303);
        Ok::<(), Box<BedErrorPlus>>(())
    })
}

// The AWS tests are skipped to credentials are not available.
#[tokio::test]
async fn s3_play_cloud() -> Result<(), Box<BedErrorPlus>> {
    use rusoto_credential::{CredentialsError, ProfileProvider, ProvideAwsCredentials};
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

    let s3 = AmazonS3Builder::new()
        .with_region("us-west-2")
        .with_url(url)
        .with_access_key_id(credentials.aws_access_key_id())
        .with_secret_access_key(credentials.aws_secret_access_key())
        .build()
        .unwrap();
    print!("{:?}", s3);

    // LATER should we accept a string as a Path?
    let store_path = StorePath::parse("/v1/toydata.5chrom.bed").unwrap();
    let object_path = ObjectPath::new(Arc::new(s3), store_path);
    assert_eq!(object_path.size().await?, 1_250_003);
    Ok(())
}

#[test]
/// Open the file and read data for one SNP (variant)
/// at index position 2.
fn read_me_cloud() -> Result<(), Box<BedErrorPlus>> {
    use bed_reader::{assert_eq_nan, BedCloud, ReadOptions, EMPTY_OPTIONS};
    use ndarray as nd;
    use {bed_reader::BedErrorPlus, tokio::runtime::Runtime}; // '#' needed for doctest
    Runtime::new().unwrap().block_on(async {
        let mut bed_cloud = BedCloud::new(
            "https://raw.githubusercontent.com/fastlmm/bed-sample-files/main/small.bed",
            EMPTY_OPTIONS,
        )
        .await?;
        let val = ReadOptions::builder()
            .sid_index(2)
            .f64()
            .read_cloud(&mut bed_cloud)
            .await?;
        assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
        Ok::<(), Box<BedErrorPlus>>(())
    })
}

#[test]
fn local_file_url_example() -> Result<(), Box<BedErrorPlus>> {
    use bed_reader::{sample_bed_file, BedCloud, ReadOptions, EMPTY_OPTIONS};
    use ndarray as nd;
    use {assert_eq_nan, bed_reader::BedErrorPlus, tokio::runtime::Runtime}; // '#' needed for doctest
    Runtime::new().unwrap().block_on(async {
        let file_name = sample_bed_file("small.bed")?.to_string_lossy().to_string();
        println!("{file_name:?}"); // For example, "C:\\Users\\carlk\\AppData\\Local\\fastlmm\\bed-reader\\cache\\small.bed"
        let url: String = path_to_url_string(file_name)?;
        println!("{url:?}"); // For example, "file:///C:/Users/carlk/AppData/Local/bed_reader/bed_reader/Cache/small.bed"
        let mut bed_cloud = BedCloud::new(url, EMPTY_OPTIONS).await?;
        let val = ReadOptions::builder()
            .sid_index(2)
            .f64()
            .read_cloud(&mut bed_cloud)
            .await?;
        assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
        Ok::<(), Box<BedErrorPlus>>(())
    })
}

#[test]
fn local_file_object_path_example() -> Result<(), Box<BedErrorPlus>> {
    use bed_reader::{sample_bed_file, BedCloud, ReadOptions};
    use ndarray as nd;
    use {assert_eq_nan, bed_reader::BedErrorPlus, tokio::runtime::Runtime}; // '#' needed for doctest
    Runtime::new().unwrap().block_on(async {
        let file_name = sample_bed_file("small.bed")?.to_string_lossy().to_string();
        println!("{file_name:?}"); // For example, "C:\\Users\\carlk\\AppData\\Local\\fastlmm\\bed-reader\\cache\\small.bed"

        let arc_object_store = Arc::new(LocalFileSystem::new());
        let path = StorePath::from_filesystem_path(&file_name)?;
        let object_path = ObjectPath::new(arc_object_store, path);

        let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
        let val = ReadOptions::builder()
            .sid_index(2)
            .f64()
            .read_cloud(&mut bed_cloud)
            .await?;
        assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);
        Ok::<(), Box<BedErrorPlus>>(())
    })
}

#[test]
fn aws_object_path_example() -> Result<(), Box<BedErrorPlus>> {
    use bed_reader::BedCloud;
    use rusoto_credential::{CredentialsError, ProfileProvider, ProvideAwsCredentials};
    use {bed_reader::BedErrorPlus, tokio::runtime::Runtime}; // '#' needed for doctest
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
        Ok::<(), Box<BedErrorPlus>>(())
    })
}

#[tokio::test]
async fn s3_article() -> Result<(), Box<BedErrorPlus>> {
    // Somehow, get AWS credentials
    use rusoto_credential::{CredentialsError, ProfileProvider, ProvideAwsCredentials};
    let credentials = if let Ok(provider) = ProfileProvider::new() {
        provider.credentials().await
    } else {
        Err(CredentialsError::new("No credentials found"))
    };
    let Ok(credentials) = credentials else {
        eprintln!("Skipping test because no AWS credentials found");
        return Ok(());
    };

    // Create a dictionary with your AWS region and credentials and the AWS region.
    let options = [
        ("aws_region", "us-west-2"),
        ("aws_access_key_id", credentials.aws_access_key_id()),
        ("aws_secret_access_key", credentials.aws_secret_access_key()),
    ];

    // Open the bed file with a URL and any needed cloud options, then use as before.
    let mut bed_cloud = BedCloud::new("s3://bedreader/v1/some_missing.bed", options).await?;
    println!("{:?}", bed_cloud.iid().await?.slice(s![..5]));
    println!("{:?}", bed_cloud.sid().await?.slice(s![..5]));
    println!(
        "{:?}",
        bed_cloud.chromosome().await?.iter().collect::<HashSet<_>>()
    );
    let val = ReadOptions::builder()
        .sid_index(bed_cloud.chromosome().await?.map(|elem| elem == "5"))
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert!(val.dim() == (100, 6));
    Ok(())
}

#[tokio::test]
async fn http_one() -> Result<(), Box<BedErrorPlus>> {
    // Open the bed file with a URL and any needed cloud options, then use as before.
    let mut bed_cloud = BedCloud::new("https://raw.githubusercontent.com/fastlmm/bed-reader/rustybed/bed_reader/tests/data/some_missing.bed", EMPTY_OPTIONS).await?;
    println!("{:?}", bed_cloud.iid().await?.slice(s![..5]));
    println!("{:?}", bed_cloud.sid().await?.slice(s![..5]));
    println!(
        "{:?}",
        bed_cloud.chromosome().await?.iter().collect::<HashSet<_>>()
    );
    let val = ReadOptions::builder()
        .sid_index(bed_cloud.chromosome().await?.map(|elem| elem == "5"))
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    assert!(val.dim() == (100, 6));
    Ok(())
}

#[tokio::test]
async fn http_two() -> Result<(), Box<BedErrorPlus>> {
    let local_fam_file = sample_file("synthetic_v1_chr-10.fam")?;
    let local_bim_file = sample_file("synthetic_v1_chr-10.bim")?;
    let empty_skip_set = HashSet::<MetadataFields>::new();
    let metadata = Metadata::new()
        .read_fam(local_fam_file, &empty_skip_set)?
        .0
        .read_bim(local_bim_file, &empty_skip_set)?
        .0;

    // Open the bed file with a URL and any needed cloud options, then use as before.
    let mut bed_cloud = BedCloud::builder(
        "https://www.ebi.ac.uk/biostudies/files/S-BSST936/genotypes/synthetic_v1_chr-10.bed",
        [("timeout", "100s")],
    )?
    .metadata(&metadata)
    .skip_early_check()
    .build()
    .await?;
    println!(
        "iid_count={}, sid_count={}",
        bed_cloud.iid_count().await?.separate_with_underscores(),
        bed_cloud.sid_count().await?.separate_with_underscores()
    );
    println!("{:?}", bed_cloud.iid().await?.slice(s![..5]));
    println!("{:?}", bed_cloud.sid().await?.slice(s![..5]));
    println!(
        "{:?}",
        bed_cloud.chromosome().await?.iter().collect::<HashSet<_>>()
    );
    let val = ReadOptions::builder()
        .iid_index(..10)
        .sid_index(s![..;bed_cloud.sid_count().await?/10])
        .f64()
        .read_cloud(&mut bed_cloud)
        .await?;
    println!("{:?}", val);
    assert!(val.dim() == (10, 10) || val.dim() == (10, 11));
    Ok(())
}

#[test]
fn http_cloud_urls_md_1() -> Result<(), Box<BedErrorPlus>> {
    use bed_reader::{BedCloud, EMPTY_OPTIONS};
    use ndarray as nd;
    use tokio::runtime::Runtime; // '#' needed for doctest
    Runtime::new()
        .unwrap()
        .block_on(async {
            let mut bed_cloud = BedCloud::new(
                "https://raw.githubusercontent.com/fastlmm/bed-sample-files/main/small.bed",
                EMPTY_OPTIONS,
            )
            .await?;
            let val: nd::Array2<f32> = bed_cloud.read().await?;
            let missing_count = val.iter().filter(|x| x.is_nan()).count();
            let missing_fraction = missing_count as f32 / val.len() as f32;
            println!("{missing_fraction:.2}"); // Outputs 0.17
            assert_eq!(missing_count, 2);
            Ok::<(), Box<dyn std::error::Error>>(())
        })
        .unwrap();
    Ok::<(), Box<BedErrorPlus>>(())
}

#[test]
fn http_cloud_urls_md_2() -> Result<(), Box<BedErrorPlus>> {
    use bed_reader::BedCloud;
    use std::collections::BTreeSet;
    use tokio::runtime::Runtime; // '#' needed for doctest
    Runtime::new()
        .unwrap()
        .block_on(async {
            let mut bed_cloud = BedCloud::builder(
                "https://raw.githubusercontent.com/fastlmm/bed-sample-files/main/toydata.5chrom.bed",
                [("timeout", "100s")],
            )?.skip_early_check().build().await?;
            println!("{:?}", bed_cloud.iid().await?.slice(s![..5])); // Outputs ndarray: ["per0", "per1", "per2", "per3", "per4"]
            println!("{:?}", bed_cloud.sid().await?.slice(s![..5])); // Outputs ndarray: ["null_0", "null_1", "null_2", "null_3", "null_4"]
            println!(
                "{:?}",
                bed_cloud
                    .chromosome()
                    .await?
                    .iter()
                    .collect::<BTreeSet<_>>()
            ); // Outputs: {"1", "2", "3", "4", "5"}
            let val = ReadOptions::builder()
                .sid_index(bed_cloud.chromosome().await?.map(|elem| elem == "5"))
                .f32()
                .read_cloud(&mut bed_cloud)
                .await?;
            assert_eq!(val.dim(), (500, 440));

            Ok::<(), Box<dyn std::error::Error>>(())
        })
        .unwrap();
    Ok::<(), Box<BedErrorPlus>>(())
}

#[test]
fn http_cloud_urls_md_3() -> Result<(), Box<BedErrorPlus>> {
    use bed_reader::BedCloud;
    use tokio::runtime::Runtime; // '#' needed for doctest
    Runtime::new()
        .unwrap()
        .block_on(async {
            let mut bed_cloud = BedCloud::builder(
                "https://www.ebi.ac.uk/biostudies/files/S-BSST936/genotypes/synthetic_v1_chr-10.bed",
                [("timeout", "100s")],
            )?
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
            Ok::<(), Box<dyn std::error::Error>>(())
        })
        .unwrap();
    Ok::<(), Box<BedErrorPlus>>(())
}

#[test]
fn http_cloud_urls_md_4() -> Result<(), Box<BedErrorPlus>> {
    use bed_reader::{assert_eq_nan, BedCloud, ObjectPath};
    use object_store::{http::HttpBuilder, ClientOptions};
    use std::time::Duration;
    use tokio::runtime::Runtime; // '#' needed for doctest
    Runtime::new()
        .unwrap()
        .block_on(async {
            let client_options = ClientOptions::new().with_timeout(Duration::from_secs(1000));
            let http_store = HttpBuilder::new()
                .with_url("https://raw.githubusercontent.com/")
                .with_client_options(client_options)
                .build()?;
            let arc_object_store = Arc::new(http_store);
            let path = StorePath::from_url_path("fastlmm/bed-sample-files/main/small.bed")?;
            let object_path = ObjectPath::new(arc_object_store, path);
            let mut bed_cloud = BedCloud::from_object_path(&object_path).await?;
            let val = ReadOptions::builder()
                .sid_index(bed_cloud.chromosome().await?.map(|elem| elem == "5"))
                .f64()
                .read_cloud(&mut bed_cloud)
                .await?;
            assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);

            Ok::<(), Box<dyn std::error::Error>>(())
        })
        .unwrap();
    Ok::<(), Box<BedErrorPlus>>(())
}
