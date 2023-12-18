use bed_reader::allclose;
use bed_reader::assert_eq_nan;
use bed_reader::assert_error_variant;
use bed_reader::sample_bed_file;
use bed_reader::sample_file;
use bed_reader::sample_files;
use bed_reader::Bed;
use bed_reader::BedError;
use bed_reader::BedErrorPlus;
use bed_reader::Metadata;
use bed_reader::MetadataFields;
use bed_reader::ReadOptions;
use bed_reader::SliceInfo1;
use bed_reader::WriteOptions;
use ndarray as nd;
use ndarray::s;
use ndarray_rand::{rand::prelude::StdRng, rand::SeedableRng, rand_distr::Uniform, RandomExt};
use std::collections::HashSet;
use std::panic::catch_unwind;
use temp_testdir::TempDir;

#[test]
fn rusty_bed1() -> Result<(), Box<BedErrorPlus>> {
    let file = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?;
    let mut bed = Bed::new(&file)?;
    let val = bed.read::<i8>()?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    assert!(mean == -13.142); // really shouldn't do mean on data where -127 represents missing

    let mut bed = Bed::new(&file)?;
    let val = ReadOptions::builder().count_a2().i8().read(&mut bed)?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    assert!(mean == -13.274); // really shouldn't do mean on data where -127 represents missing
    Ok(())
}

#[test]
fn rusty_bed2() -> Result<(), Box<BedErrorPlus>> {
    let file = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?;
    let mut bed = Bed::new(file)?;

    let val = ReadOptions::builder()
        .iid_index(0)
        .sid_index(vec![1])
        .i8()
        .read(&mut bed)?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == 1.0); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[test]
fn rusty_bed3() -> Result<(), Box<BedErrorPlus>> {
    let file = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?;
    let mut bed = Bed::new(file)?;
    let iid_bool: nd::Array1<bool> = (0..bed.iid_count()?).map(|elem| (elem % 2) != 0).collect();
    let sid_bool: nd::Array1<bool> = (0..bed.sid_count()?).map(|elem| (elem % 8) != 0).collect();
    let val = ReadOptions::builder()
        .missing_value(-127)
        .iid_index(iid_bool)
        .sid_index(sid_bool)
        .read(&mut bed)?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == -14.50344827586207); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[test]
fn rusty_bed_allele() -> Result<(), Box<BedErrorPlus>> {
    let file = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?;
    let mut bed = Bed::new(file)?;
    let val = ReadOptions::builder().count_a2().i8().read(&mut bed)?;

    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == -13.274); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[test]
fn rusty_bed_order() -> Result<(), Box<BedErrorPlus>> {
    let file = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?;
    let mut bed = Bed::new(file)?;
    let val = ReadOptions::builder().c().i8().read(&mut bed)?;

    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == -13.142); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[test]
fn bad_header() -> Result<(), Box<BedErrorPlus>> {
    let filename = sample_file("badfile.bed")?;
    let bed = Bed::builder(&filename).skip_early_check().build()?;
    println!("{:?}", bed.path());

    let result = Bed::new(&filename);
    assert_error_variant!(result, BedErrorPlus::BedError(BedError::IllFormed(_)));

    Ok(())
}

#[test]
fn doc_test_test() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(file_name)?;
    let val = bed.read::<f64>()?;
    assert_eq_nan(
        &val,
        &nd::array![
            [1.0, 0.0, f64::NAN, 0.0],
            [2.0, 0.0, f64::NAN, 2.0],
            [0.0, 1.0, 2.0, 0.0]
        ],
    );

    let file_name2 = sample_bed_file("some_missing.bed")?;
    let mut bed2 = Bed::new(&file_name2)?;
    let val2 = ReadOptions::builder()
        .f64()
        .iid_index(s![..;2])
        .sid_index(20..30)
        .read(&mut bed2)?;
    assert!(val2.dim() == (50, 10));

    let mut bed3 = Bed::new(&file_name2)?;
    println!("{:?}", bed3.iid()?.slice(s![..5]));
    println!("{:?}", bed3.sid()?.slice(s![..5]));
    println!("{:?}", bed3.chromosome()?.iter().collect::<HashSet<_>>());
    let val3 = ReadOptions::builder()
        .sid_index(bed3.chromosome()?.map(|elem| elem == "5"))
        .f64()
        .read(&mut bed3)?;
    assert!(val3.dim() == (100, 6));

    Ok(())
}

#[test]
fn open_examples() -> Result<(), Box<BedErrorPlus>> {
    //     >>> from bed_reader import open_bed, sample_bed_file
    //     >>>
    //     >>> file_name = sample_bed_file("small.bed")
    //     >>> bed = open_bed(file_name)
    //     >>> print(bed.iid)
    //     ['iid1' 'iid2' 'iid3']
    //     >>> print(bed.sid)
    //     ['sid1' 'sid2' 'sid3' 'sid4']
    //     >>> print(bed.read())
    //     [[ 1.  0. nan  0.]
    //      [ 2.  0. nan  2.]
    //      [ 0.  1.  2.  0.]]
    //     >>> del bed  # optional: delete bed object

    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(&file_name)?;
    println!("{:?}", bed.iid()?);
    println!("{:?}", bed.sid()?);
    println!("{:?}", bed.read::<f64>()?);

    // ["iid1", "iid2", "iid3"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1
    // ["sid1", "sid2", "sid3", "sid4"], shape=[4], strides=[1], layout=CFcf (0xf), const ndim=1
    // [[1.0, 0.0, NaN, 0.0],
    //  [2.0, 0.0, NaN, 2.0],
    //  [0.0, 1.0, 2.0, 0.0]], shape=[3, 4], strides=[1, 3], layout=Ff (0xa), const ndim=2

    // Open the file and read data for one SNP (variant)
    // at index position 2.

    // .. doctest::

    //     >>> import numpy as np
    //     >>> with open_bed(file_name) as bed:
    //     ...     print(bed.read(np.s_[:,2]))
    //     [[nan]
    //      [nan]
    //      [ 2.]]

    let mut bed = Bed::new(&file_name)?;
    println!(
        "{:?}",
        ReadOptions::builder().sid_index(2).f64().read(&mut bed)?
    );

    // [[NaN],
    //  [NaN],
    //  [2.0]], shape=[3, 1], strides=[1, 3], layout=CFcf (0xf), const ndim=2

    // Replace :attr:`iid`.

    //     >>> bed = open_bed(file_name, properties={"iid":["sample1","sample2","sample3"]})
    //     >>> print(bed.iid) # replaced
    //     ['sample1' 'sample2' 'sample3']
    //     >>> print(bed.sid) # same as before
    //     ['sid1' 'sid2' 'sid3' 'sid4']

    let mut bed = Bed::builder(&file_name)
        .iid(["sample1", "sample2", "sample3"])
        .build()?;
    println!("{:?}", bed.iid()?);
    println!("{:?}", bed.sid()?);

    // ["sample1", "sample2", "sample3"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1
    // ["sid1", "sid2", "sid3", "sid4"], shape=[4], strides=[1], layout=CFcf (0xf), const ndim=

    // Do more testing of Rust
    let iid = nd::array!["sample1", "sample2", "sample3"];
    let mut _bed = Bed::builder(&file_name).iid(iid).build()?;
    let iid = nd::array![
        "sample1".to_string(),
        "sample2".to_string(),
        "sample3".to_string()
    ];
    let mut _bed = Bed::builder(&file_name).iid(iid).build()?;
    let iid = vec!["sample1", "sample2", "sample3"];
    let mut _bed = Bed::builder(&file_name).iid(iid).build()?;
    let iid = vec![
        "sample1".to_string(),
        "sample2".to_string(),
        "sample3".to_string(),
    ];
    let mut _bed = Bed::builder(&file_name).iid(iid).build()?;

    // Give the number of individuals (samples) and SNPs (variants) so that the .fam and
    // .bim files need never be opened.

    //     >>> with open_bed(file_name, iid_count=3, sid_count=4) as bed:
    //     ...     print(bed.read())
    //     [[ 1.  0. nan  0.]
    //      [ 2.  0. nan  2.]
    //      [ 0.  1.  2.  0.]]

    let mut bed = Bed::builder(&file_name).iid_count(3).sid_count(4).build()?;
    println!("{:?}", bed.read::<f64>()?);

    //  [[1.0, 0.0, NaN, 0.0],
    //   [2.0, 0.0, NaN, 2.0],
    //   [0.0, 1.0, 2.0, 0.0]], shape=[3, 4], strides=[1, 3], layout=Ff (0xa), const ndim=2

    // Mark some properties as "don’t read or offer".

    //     >>> bed = open_bed(file_name, properties={
    //     ...    "father" : None, "mother" : None, "sex" : None, "pheno" : None,
    //     ...    "allele_1" : None, "allele_2":None })
    //     >>> print(bed.iid)        # read from file
    //     ['iid1' 'iid2' 'iid3']
    //     >>> print(bed.allele_2)   # not read and not offered
    //     None

    let mut bed = Bed::builder(&file_name).skip_allele_2().build()?;
    println!("{:?}", bed.iid()?);

    let result = bed.allele_2();
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::CannotUseSkippedMetadata(_))
    );

    Ok(())
}

#[test]
fn metadata_etc() -> Result<(), Box<BedErrorPlus>> {
    // >>> file_name = sample_bed_file("small.bed")
    // >>> with open_bed(file_name) as bed:
    // ...     print(bed.sex)
    // [1 2 0]
    let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
    println!("{:?}", bed.sex()?);
    // [1, 2, 0], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1

    let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
    println!("{:?}", bed.cm_position()?);
    // [100.4, 2000.5, 4000.7, 7000.9], shape=[4], strides=[1], layout=CFcf (0xf), const ndim=1

    println!("{:?}", bed.bp_position()?);
    // [1, 100, 1000, 1004], shape=[4], strides=[1], layout=CFcf (0xf), const ndim=1

    let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
    println!("{:?}", bed.fid()?);
    // ["fid1", "fid1", "fid2"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1

    println!("{:?}", bed.father()?);
    // ["iid23", "iid23", "iid22"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1

    println!("{:?}", bed.mother()?);
    // ["iid34", "iid34", "iid33"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1

    Ok(())
}

#[test]
fn hello_father() -> Result<(), Box<BedErrorPlus>> {
    let mut bed = Bed::builder(sample_bed_file("small.bed")?)
        .father(["f1", "f2", "f3"])
        .skip_mother()
        .build()?;
    println!("{:?}", bed.father()?);
    // ["f1", "f2", "f3"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1
    bed.mother().unwrap_err();

    Ok(())
}

#[test]
fn num_threads() -> Result<(), Box<BedErrorPlus>> {
    let file = sample_bed_file("plink_sim_10s_100v_10pmiss.bed")?;
    let mut bed = Bed::new(file)?;

    let val = ReadOptions::builder().num_threads(4).i8().read(&mut bed)?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == -13.142); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[test]
fn fam_and_bim() -> Result<(), Box<BedErrorPlus>> {
    let deb_maf_mib = sample_files(["small.deb", "small.maf", "small.mib"])?;
    let mut bed = Bed::builder(&deb_maf_mib[0])
        .fam_path(&deb_maf_mib[1])
        .bim_path(&deb_maf_mib[2])
        .build()?;

    println!("{:?}", bed.iid()?);
    println!("{:?}", bed.sid()?);
    let val: nd::Array2<i8> = bed.read()?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{mean:?}");
    assert!(mean == -20.5); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[test]
fn readme_examples() -> Result<(), Box<BedErrorPlus>> {
    // Read genomic data from a .bed file.

    // >>> import numpy as np
    // >>> from bed_reader import open_bed, sample_bed_file
    // >>>
    // >>> file_name = sample_bed_file("small.bed")
    // >>> bed = open_bed(file_name)
    // >>> val = bed.read()
    // >>> print(val)
    // [[ 1.  0. nan  0.]
    //  [ 2.  0. nan  2.]
    //  [ 0.  1.  2.  0.]]
    // >>> del bed

    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(file_name)?;
    let val = bed.read::<f64>()?;
    println!("{val:?}");
    // [[1.0, 0.0, NaN, 0.0],
    // [2.0, 0.0, NaN, 2.0],
    // [0.0, 1.0, 2.0, 0.0]], shape=[3, 4], strides=[1, 3], layout=Ff (0xa), const ndim=2

    // Read every second individual and SNPs (variants) from 20 to 30.

    // >>> file_name2 = sample_bed_file("some_missing.bed")
    // >>> bed2 = open_bed(file_name2)
    // >>> val2 = bed2.read(index=np.s_[::2,20:30])
    // >>> print(val2.shape)
    // (50, 10)
    // >>> del bed2

    let file_name2 = sample_bed_file("some_missing.bed")?;
    let mut bed2 = Bed::new(&file_name2)?;
    let val2 = ReadOptions::<f64>::builder()
        .iid_index(s![..;2])
        .sid_index(20..30)
        .read(&mut bed2)?;
    println!("{:?}", val2.dim()); // (50, 10)

    // List the first 5 individual (sample) ids, the first 5 SNP (variant) ids,
    // and every unique chromosome. Then, read every genomic value in chromosome 5.

    // >>> with open_bed(file_name2) as bed3:
    // ...     print(bed3.iid[:5])
    // ...     print(bed3.sid[:5])
    // ...     print(np.unique(bed3.chromosome))
    // ...     val3 = bed3.read(index=np.s_[:,bed3.chromosome=='5'])
    // ...     print(val3.shape)
    // ['iid_0' 'iid_1' 'iid_2' 'iid_3' 'iid_4']
    // ['sid_0' 'sid_1' 'sid_2' 'sid_3' 'sid_4']
    // ['1' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '2' '20' '21' '22'
    //  '3' '4' '5' '6' '7' '8' '9']
    // (100, 6)

    let mut bed3 = Bed::new(&file_name2)?;
    let iid = bed3.iid()?;
    let s = iid.slice(s![..5]);
    println!("{:?}", s);
    println!("{:?}", bed3.iid()?.slice(s![..5]));
    println!("{:?}", bed3.sid()?.slice(s![..5]));
    let unique = bed3.chromosome()?.iter().collect::<HashSet<_>>();
    println!("{unique:?}");
    // let is_5 = bed3.chromosome()?.map(|elem| elem == "5");
    let is_5 = nd::Zip::from(bed3.chromosome()?).par_map_collect(|elem| elem == "5");
    let val3 = ReadOptions::builder()
        .sid_index(is_5)
        .f64()
        .read(&mut bed3)?;
    println!("{:?}", val3.dim());
    // ["iid_0", "iid_1", "iid_2", "iid_3", "iid_4"], shape=[5], strides=[1], layout=CFcf (0xf), const ndim=1
    // ["sid_0", "sid_1", "sid_2", "sid_3", "sid_4"], shape=[5], strides=[1], layout=CFcf (0xf), const ndim=1
    // {"10", "11", "4", "21", "22", "14", "3", "12", "20", "15", "19", "8", "6", "18", "9", "2", "16", "13", "17", "1", "7", "5"}
    // (100, 6)
    Ok(())
}

#[test]
fn write_docs() -> Result<(), Box<BedErrorPlus>> {
    // In this example, all properties are given.

    //     >>> from bed_reader import to_bed, tmp_path
    //     >>>
    //     >>> output_file = tmp_path() / "small.bed"
    //     >>> val = [[1.0, 0.0, np.nan, 0.0],
    //     ...        [2.0, 0.0, np.nan, 2.0],
    //     ...        [0.0, 1.0, 2.0, 0.0]]
    //     >>> properties = {
    //     ...    "fid": ["fid1", "fid1", "fid2"],
    //     ...    "iid": ["iid1", "iid2", "iid3"],
    //     ...    "father": ["iid23", "iid23", "iid22"],
    //     ...    "mother": ["iid34", "iid34", "iid33"],
    //     ...    "sex": [1, 2, 0],
    //     ...    "pheno": ["red", "red", "blue"],
    //     ...    "chromosome": ["1", "1", "5", "Y"],
    //     ...    "sid": ["sid1", "sid2", "sid3", "sid4"],
    //     ...    "cm_position": [100.4, 2000.5, 4000.7, 7000.9],
    //     ...    "bp_position": [1, 100, 1000, 1004],
    //     ...    "allele_1": ["A", "T", "A", "T"],
    //     ...    "allele_2": ["A", "C", "C", "G"],
    //     ... }
    //     >>> to_bed(output_file, val, properties=properties)

    let output_folder = TempDir::default();

    let output_file = output_folder.join("small.bed");
    let val = nd::array![
        [1.0, 0.0, f64::NAN, 0.0],
        [2.0, 0.0, f64::NAN, 2.0],
        [0.0, 1.0, 2.0, 0.0]
    ];
    WriteOptions::builder(output_file)
        .fid(["fid1", "fid1", "fid2"])
        .iid(["iid1", "iid2", "iid3"])
        .father(["iid23", "iid23", "iid22"])
        .mother(["iid34", "iid34", "iid33"])
        .sex([1, 2, 0])
        .pheno(["red", "red", "blue"])
        .chromosome(["1", "1", "5", "Y"])
        .sid(["sid1", "sid2", "sid3", "sid4"])
        .cm_position([100.4, 2000.5, 4000.7, 7000.9])
        .bp_position([1, 100, 1000, 1004])
        .allele_1(["A", "T", "A", "T"])
        .allele_2(["A", "C", "C", "G"])
        .write(&val)?;

    // Here, no properties are given, so default values are assigned.
    // If we then read the new file and list the chromosome property,
    // it is an array of '0's, the default chromosome value.

    //     >>> output_file2 = tmp_path() / "small2.bed"
    //     >>> val = [[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]]
    //     >>> to_bed(output_file2, val)
    //     >>>
    //     >>> from bed_reader import open_bed
    //     >>> with open_bed(output_file2) as bed2:
    //     ...     print(bed2.chromosome)
    //     ['0' '0' '0' '0']
    let output_file2 = output_folder.join("small2.bed");
    let val = nd::array![[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]];
    Bed::write(&val, &output_file2)?;
    let mut bed2 = Bed::new(&output_file2)?;
    println!("{:?}", bed2.chromosome()?);
    // ["0", "0", "0", "0"], shape=[4], strides=[1], layout=CFcf (0xf), const ndim=1

    Ok(())
}

#[test]
fn read_write() -> Result<(), Box<BedErrorPlus>> {
    // with open_bed(shared_datadir / "small.bed") as bed:
    // val = bed.read()
    // properties = bed.properties

    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(file_name)?;
    let val = bed.read::<f64>()?;
    let metadata = bed.metadata()?;
    println!("{metadata:?}");

    // output_file = tmp_path / "small.deb"
    // fam_file = tmp_path / "small.maf"
    // bim_file = tmp_path / "small.mib"

    let temp_out = TempDir::default();
    let output_file = temp_out.join("small.deb");
    let fam_file = temp_out.join("small.maf");
    let bim_file = temp_out.join("small.mib");

    // to_bed(
    // output_file,
    // val,
    // properties=properties,
    // fam_filepath=fam_file,
    // bim_filepath=bim_file,
    // )

    WriteOptions::builder(&output_file)
        .metadata(&metadata)
        .fam_path(&fam_file)
        .bim_path(&bim_file)
        .write(&val)?;

    // // assert output_file.exists() and fam_file.exists() and bim_file.exists()
    assert!(
        output_file.exists() & fam_file.exists() & bim_file.exists(),
        "don't exist"
    );

    // with open_bed(output_file, fam_filepath=fam_file, bim_filepath=bim_file) as deb:
    // val2 = deb.read()
    // properties2 = deb.properties
    let mut deb = Bed::builder(&output_file)
        .fam_path(&fam_file)
        .bim_path(&bim_file)
        .build()?;
    let val2 = deb.read::<f64>()?;
    let metadata2 = deb.metadata()?;

    // assert np.allclose(val, val2, equal_nan=true)
    assert!(
        allclose(&val.view(), &val2.view(), 1e-08, true),
        "not close"
    );
    println!("{metadata:?}");
    println!("{metadata2:?}");
    assert!(metadata == metadata2, "meta not equal");

    Ok(())
}

#[test]
fn range() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(file_name)?;
    ReadOptions::builder().iid_index(0..2).i8().read(&mut bed)?;
    ReadOptions::builder()
        .iid_index(0..=2)
        .i8()
        .read(&mut bed)?;
    ReadOptions::builder().iid_index(..2).i8().read(&mut bed)?;
    ReadOptions::builder().iid_index(..=2).i8().read(&mut bed)?;
    ReadOptions::builder().iid_index(0..).i8().read(&mut bed)?;
    ReadOptions::builder().iid_index(..).i8().read(&mut bed)?;

    Ok(())
}

#[test]
fn nd_slice() -> Result<(), Box<BedErrorPlus>> {
    let ndarray = nd::array![0, 1, 2, 3];
    println!("{:?}", ndarray.slice(nd::s![1..3])); // [1, 2]
    println!("{:?}", ndarray.slice(nd::s![1..3;-1])); // [2, 1]
    #[allow(clippy::reversed_empty_ranges)]
    let slice = nd::s![3..1;-1];
    println!("{:?}", ndarray.slice(slice)); // []

    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(file_name)?;
    ReadOptions::builder()
        .iid_index(nd::s![0..2])
        .i8()
        .read(&mut bed)?;
    ReadOptions::builder()
        .iid_index(nd::s![..2])
        .i8()
        .read(&mut bed)?;
    ReadOptions::builder()
        .iid_index(nd::s![0..])
        .i8()
        .read(&mut bed)?;
    ReadOptions::builder()
        .iid_index(nd::s![0..2])
        .i8()
        .read(&mut bed)?;
    ReadOptions::builder()
        .iid_index(nd::s![-2..-1;-1])
        .i8()
        .read(&mut bed)?;

    Ok(())
}

#[test]
fn skip_coverage() -> Result<(), Box<BedErrorPlus>> {
    let mut bed = Bed::builder(sample_bed_file("small.bed")?)
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
        .build()?;
    bed.mother().unwrap_err();

    Ok(())
}

#[test]
fn into_iter() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::builder(file_name)
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
        .build()?;

    let _ = bed.pheno()?;
    Ok(())
}

#[allow(clippy::reversed_empty_ranges)]
#[test]
fn range_same() -> Result<(), Box<BedErrorPlus>> {
    let a = rt1(3..0);
    println!("{:?}", a);
    let b = rt23((3..0).into());
    println!("{:?}", b);

    assert_same_result(rt1(3..0), rt23((3..0).into()));
    assert_same_result(rt1(1000..), rt23((1000..).into()));

    assert_same_result(rt1(..), rt23((..).into()));
    assert_same_result(rt1(..3), rt23((..3).into()));
    assert_same_result(rt1(..=3), rt23((..=3).into()));
    assert_same_result(rt1(1..), rt23((1..).into()));
    assert_same_result(rt1(1..3), rt23((1..3).into()));
    assert_same_result(rt1(1..=3), rt23((1..=3).into()));
    assert_same_result(rt1(2..=2), rt23((2..=2).into()));
    Ok(())
}

#[test]
fn nd_slice_same() -> Result<(), Box<BedErrorPlus>> {
    assert_same_result(nds1(s![1000..]), rt23(s![1000..].into()));
    assert_same_result(nds1(s![..1000]), rt23(s![..1000].into()));
    assert_same_result(nds1(s![999..1000]), rt23(s![999..1000].into()));
    assert_same_result(nds1(s![-1000..]), rt23(s![-1000..].into()));
    assert_same_result(nds1(s![..-1000]), rt23(s![..-1000].into()));
    #[allow(clippy::reversed_empty_ranges)]
    assert_same_result(nds1(s![-999..-1000]), rt23(s![-999..-1000].into()));
    #[allow(clippy::reversed_empty_ranges)]
    assert_same_result(nds1(s![3..0]), rt23(s![3..0].into()));
    #[allow(clippy::reversed_empty_ranges)]
    assert_same_result(nds1(s![-1..-2]), rt23(s![-1..-2].into()));

    assert_same_result(nds1(s![..-3]), rt23(s![..-3].into()));
    assert_same_result(nds1(s![..=-3]), rt23(s![..=-3].into()));
    assert_same_result(nds1(s![-1..]), rt23(s![-1..].into()));
    assert_same_result(nds1(s![-3..-1]), rt23(s![-3..-1].into()));
    assert_same_result(nds1(s![-3..=-1]), rt23(s![-3..=-1].into()));
    assert_same_result(nds1(s![-3..=-1]), rt23(s![-3..=-1].into()));
    assert_same_result(nds1(s![-2..=-2]), rt23(s![-2..=-2].into()));
    #[allow(clippy::reversed_empty_ranges)]
    assert_same_result(nds1(s![1..-1]), rt23(s![1..-1].into()));

    assert_same_result(nds1(s![..]), rt23((s![..]).into()));
    assert_same_result(nds1(s![..3]), rt23((s![..3]).into()));
    assert_same_result(nds1(s![..=3]), rt23((s![..=3]).into()));
    assert_same_result(nds1(s![1..]), rt23((s![1..]).into()));
    assert_same_result(nds1(s![1..3]), rt23((s![1..3]).into()));
    assert_same_result(nds1(s![1..=3]), rt23((s![1..=3]).into()));
    assert_same_result(nds1(s![2..=2]), rt23(s![2..=2].into()));

    Ok(())
}

#[test]
fn counts_and_files() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;

    // We give the wrong number for iid_count and then expect an error
    let mut bed = Bed::builder(&file_name)
        .iid(["i1", "i2", "i3", "i4"])
        .build()?;
    let iid_count = bed.iid_count()?;
    assert_eq!(iid_count, 4);
    let fid_result = bed.fid();
    assert_error_variant!(
        fid_result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    assert_error_variant!(
        Bed::builder(&file_name)
            .fid(["f1", "f1", "f1"])
            .iid(["i1", "i2", "i3", "i4"])
            .build(),
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    let mut bed = Bed::builder(&file_name)
        .bim_path(sample_file("small.bad_bim")?)
        .build()?;

    assert_error_variant!(
        bed.sid(),
        BedErrorPlus::BedError(BedError::MetadataFieldCount(_, _, _))
    );

    // We give the wrong number for iid_count and then expect an error
    let mut bed = Bed::builder(&file_name).iid_count(4).build()?;
    assert_eq!(bed.iid_count()?, 4);
    assert_error_variant!(
        bed.iid(),
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    let mut bed = Bed::builder(&file_name)
        .bim_path(sample_file("small.bim")?)
        .build()?;
    assert_eq!(bed.iid_count()?, 3);
    assert_eq!(bed.sid_count()?, 4);
    let mut bed = Bed::new(&file_name)?;
    assert_eq!(bed.iid_count()?, 3);
    assert_eq!(bed.sid_count()?, 4);

    let mut bed = Bed::builder(&file_name).build()?;
    let _ = bed.iid()?;
    let _ = bed.iid_count()?;

    Ok(())
}

#[test]
fn bool_read() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(file_name)?;
    let result = ReadOptions::builder()
        .iid_index([false, false, true, false])
        .i8()
        .read(&mut bed);
    println!("{result:?}");
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::BoolArrayVectorWrongLength(_, _))
    );

    let _val = ReadOptions::builder()
        .iid_index([false, false, true])
        .i8()
        .read(&mut bed)?;

    Ok(())
}

#[test]
fn i8_etc() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(file_name)?;
    let _val = ReadOptions::builder()
        .f()
        .i8()
        .iid_index([false, false, true])
        .read(&mut bed)?;

    Ok(())
}

#[test]
fn fill() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(file_name)?;
    let read_options = ReadOptions::builder()
        .f()
        .i8()
        .iid_index([false, false, true])
        .build()?;

    let mut val = nd::Array2::<i8>::default((3, 4));
    let result = bed.read_and_fill_with_options(&mut val.view_mut(), &read_options);
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InvalidShape(_, _, _, _))
    );

    let mut val = nd::Array2::<i8>::default((1, 4));
    bed.read_and_fill_with_options(&mut val.view_mut(), &read_options)?;

    assert_eq!(bed.dim()?, (3, 4));

    Ok(())
}

#[test]
fn read_options_builder() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;

    let mut bed = Bed::new(file_name)?;
    // Read the SNPs indexed by 2.
    let val = ReadOptions::builder().sid_index(2).f64().read(&mut bed)?;

    assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);

    // Read the SNPs indexed by 2, 3, and 0.
    let val = ReadOptions::builder()
        .sid_index([2, 3, 0])
        .f64()
        .read(&mut bed)?;

    assert_eq_nan(
        &val,
        &nd::array![[f64::NAN, 0.0, 1.0], [f64::NAN, 2.0, 2.0], [2.0, 0.0, 0.0]],
    );

    //  Read SNPs from 1 (inclusive) to 4 (exclusive).
    let val = ReadOptions::builder()
        .sid_index(1..4)
        .f64()
        .read(&mut bed)?;

    assert_eq_nan(
        &val,
        &nd::array![[0.0, f64::NAN, 0.0], [0.0, f64::NAN, 2.0], [1.0, 2.0, 0.0]],
    );

    // Print unique chrom values. Then, read all SNPs in chrom 5.
    use std::collections::HashSet;

    println!("{:?}", bed.chromosome()?.iter().collect::<HashSet<_>>());
    // This outputs: {"5", "1", "Y"}.
    let val = ReadOptions::builder()
        .sid_index(bed.chromosome()?.map(|elem| elem == "5"))
        .f64()
        .read(&mut bed)?;

    assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);

    // Read 1st individual (across all SNPs).
    let val = ReadOptions::builder().iid_index(0).f64().read(&mut bed)?;

    assert_eq_nan(&val, &nd::array![[1.0, 0.0, f64::NAN, 0.0]]);

    // Read every 2nd individual.
    use ndarray::s;
    let val = ReadOptions::builder()
        .iid_index(s![..;2])
        .f64()
        .read(&mut bed)?;

    assert_eq_nan(
        &val,
        &nd::array![[1.0, 0.0, f64::NAN, 0.0], [0.0, 1.0, 2.0, 0.0]],
    );

    // Read last and 2nd-to-last individuals and the last SNPs
    let val = ReadOptions::builder()
        .iid_index([-1, -2])
        .sid_index(-1)
        .f64()
        .read(&mut bed)?;

    println!("{:?}", &val);
    assert_eq_nan(&val, &nd::array![[0.0], [2.0]]);
    Ok(())
}

#[test]
fn bed_builder() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::builder(&file_name).build()?;
    println!("{:?}", bed.iid()?);
    println!("{:?}", bed.sid()?);
    let val = bed.read::<f64>()?;

    assert_eq_nan(
        &val,
        &nd::array![
            [1.0, 0.0, f64::NAN, 0.0],
            [2.0, 0.0, f64::NAN, 2.0],
            [0.0, 1.0, 2.0, 0.0]
        ],
    );

    let mut bed = Bed::builder(&file_name).build()?;
    let val = ReadOptions::builder().sid_index(2).f64().read(&mut bed)?;

    assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);

    let mut bed = Bed::builder(&file_name)
        .iid(["sample1", "sample2", "sample3"])
        .build()?;
    println!("{:?}", bed.iid()?); // replaced
    println!("{:?}", bed.sid()?); // same as before

    let mut bed = Bed::builder(&file_name).iid_count(3).sid_count(4).build()?;
    let val = bed.read::<f64>()?;
    assert_eq_nan(
        &val,
        &nd::array![
            [1.0, 0.0, f64::NAN, 0.0],
            [2.0, 0.0, f64::NAN, 2.0],
            [0.0, 1.0, 2.0, 0.0]
        ],
    );

    let mut bed = Bed::builder(&file_name)
        .skip_father()
        .skip_mother()
        .skip_sex()
        .skip_pheno()
        .skip_allele_1()
        .skip_allele_2()
        .build()?;
    println!("{:?}", bed.iid()?);
    bed.allele_2().expect_err("Can't be read");

    Ok::<(), Box<BedErrorPlus>>(())
}

#[test]
#[allow(clippy::needless_borrows_for_generic_args)]
fn negative_indexing() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(&file_name)?;
    // println!("{:?}", bed.read::<f64>()?);
    // [[1.0, 0.0, NaN, 0.0],
    // [2.0, 0.0, NaN, 2.0],
    // [0.0, 1.0, 2.0, 0.0]], shape=[3, 4], strides=[1, 3], layout=Ff (0xa), const ndim=2
    //  iid range is -4ERROR -3 -2 -1 0 1 2 3ERROR
    //  sid range is -5ERROR -4 ... 3 4ERROR
    for index in [-4, 3] {
        match ReadOptions::builder().iid_index(index).i8().read(&mut bed) {
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
            .read(&mut bed)?;
        // println!("{val:?}");
        assert!(val[[0, 0]] == 1,);
    }
    for index in [-1, 2] {
        let val = ReadOptions::builder()
            .iid_index(index)
            .i8()
            .read(&mut bed)?;
        // println!("{val:?}");
        assert!(val[[0, 0]] == 0,);
    }

    for index in [-5, 4] {
        match ReadOptions::builder().sid_index(index).i8().read(&mut bed) {
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
            .read(&mut bed)?;
        // println!("{val:?}");
        assert!(val[[0, 0]] == 1,);
    }
    for index in [-1, 3] {
        let val = ReadOptions::builder()
            .sid_index(index)
            .i8()
            .read(&mut bed)?;
        // println!("{val:?}");
        assert!(val[[0, 0]] == 0,);
    }

    Ok(())
}

#[test]
fn index_doc() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("some_missing.bed")?;
    let mut bed = Bed::new(file_name)?;
    println!("{:?}", bed.dim()?); // prints (100, 100)

    // Read all individuals and all SNPs
    let val = ReadOptions::builder().f64().read(&mut bed)?;
    assert!(val.dim() == (100, 100));

    // Read the individual at index position 10 and all SNPs
    let val = ReadOptions::builder().iid_index(10).f64().read(&mut bed)?;
    assert!(val.dim() == (1, 100));

    // Read the individuals at index positions 0,5, 1st-from-the-end and
    // the SNP at index position 3
    let val = ReadOptions::builder()
        .iid_index(vec![0, 5, -1])
        .sid_index(3)
        .f64()
        .read(&mut bed)?;
    assert!(val.dim() == (3, 1));
    // Repeat, but with an ndarray
    let val = ReadOptions::builder()
        .iid_index(nd::array![0, 5, -1])
        .sid_index(3)
        .f64()
        .read(&mut bed)?;
    assert!(val.dim() == (3, 1));
    // Repeat, but with an Rust array
    let val = ReadOptions::builder()
        .iid_index([0, 5, -1])
        .sid_index(3)
        .f64()
        .read(&mut bed)?;
    assert!(val.dim() == (3, 1));

    // Create a boolean ndarray identifying SNPs in chromosome 5,
    // then select those SNPs.
    let chrom_5 = bed.chromosome()?.map(|elem| elem == "5");
    let val = ReadOptions::builder()
        .sid_index(chrom_5)
        .f64()
        .read(&mut bed)?;
    assert!(val.dim() == (100, 6));

    // Use ndarray's slice macro, [`s!`](https://docs.rs/ndarray/latest/ndarray/macro.s.html),
    // to select every 2nd individual and every 3rd SNP.
    let val = ReadOptions::builder()
        .iid_index(s![..;2])
        .sid_index(s![..;3])
        .f64()
        .read(&mut bed)?;
    assert!(val.dim() == (50, 34));
    // Use ndarray's slice macro, [`s!`](https://docs.rs/ndarray/latest/ndarray/macro.s.html),
    // to select the 10th-from-last individual to the last, in reverse order,
    // and every 3rd SNP in reverse order.)
    let val = ReadOptions::builder()
        .iid_index(s![-10..;-1])
        .sid_index(s![..;-3])
        .f64()
        .read(&mut bed)?;
    assert!(val.dim() == (10, 34));
    Ok(())
}

#[test]
fn index_options() -> Result<(), Box<BedErrorPlus>> {
    let mut bed = Bed::new(sample_bed_file("some_missing.bed")?)?;
    #[allow(clippy::let_unit_value)]
    let index: () = ();

    let all = ReadOptions::builder()
        .iid_index(index)
        .sid_index(index)
        .f64()
        .read(&mut bed)?;
    assert!(all.dim() == (100, 100));

    let mut index: [bool; 100] = [false; 100];
    index[0] = true;
    index[2] = true;
    let val = ReadOptions::builder()
        .iid_index(index)
        .sid_index(index)
        .f64()
        .read(&mut bed)?;
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
        .read(&mut bed)?;
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
        .read(&mut bed)?;
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
        .read(&mut bed)?;
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
        .read(&mut bed)?;
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
        .read(&mut bed)?;
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
        .read(&mut bed)?;
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
        .read(&mut bed)?;
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
        .read(&mut bed)?;
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
        .read(&mut bed)?;
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
        .read(&mut bed)?;
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
        .read(&mut bed)?;
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
        .read(&mut bed)?;
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
        .read(&mut bed)?;
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
        .read(&mut bed)?;
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
        .read(&mut bed)?;
    let expected = all.slice(s![-20..-10;-2,-20..-10;-2]);
    assert!(
        allclose(&val.view(), &expected.view(), 1e-08, true),
        "not close"
    );

    Ok(())
}

#[test]
fn set_metadata() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;
    let metadata = Metadata::builder()
        .iid(["iid1", "iid2", "iid3"])
        .sid(["sid1", "sid2", "sid3", "sid4"])
        .build()?;
    let mut bed = Bed::builder(&file_name).metadata(&metadata).build()?;
    let metadata2 = bed.metadata()?;
    println!("{metadata2:?}");

    let mut bed = Bed::new(&file_name)?;
    let metadata = bed.metadata()?;
    println!("{metadata:?}");

    let mut bed = Bed::builder(&file_name).metadata(&metadata).build()?;
    let metadata2 = bed.metadata()?;
    println!("{metadata2:?}");

    Ok(())
}

#[test]
fn metadata_print() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(file_name)?;

    let fid = bed.fid()?;
    println!("{fid:?}"); // Outputs ndarray ["fid1", "fid1", "fid2"]
    let iid = bed.iid()?;
    println!("{iid:?}"); // Outputs ndarray ["iid1", "iid2", "iid3"]
    let father = bed.father()?;
    println!("{father:?}"); // Outputs ndarray ["iid23", "iid23", "iid22"]
    let mother = bed.mother()?;
    println!("{mother:?}"); // Outputs ndarray ["iid34", "iid34", "iid33"]
    let sex = bed.sex()?;
    println!("{sex:?}"); // Outputs ndarray [1, 2, 0]
    let pheno = bed.pheno()?;
    println!("{pheno:?}"); // Outputs ndarray ["red", "red", "blue"]

    let chromosome = bed.chromosome()?;
    println!("{chromosome:?}"); // Outputs ndarray ["1", "1", "5", "Y"
    let sid = bed.sid()?;
    println!("{sid:?}"); // Outputs ndarray "sid1", "sid2", "sid3", "sid4"]
    let cm_position = bed.cm_position()?;
    println!("{cm_position:?}"); // Outputs ndarray [100.4, 2000.5, 4000.7, 7000.9]
    let bp_position = bed.bp_position()?;
    println!("{bp_position:?}"); // Outputs ndarray [1, 100, 1000, 1004]
    let allele_1 = bed.allele_1()?;
    println!("{allele_1:?}"); // Outputs ndarray ["A", "T", "A", "T"]
    let allele_2 = bed.allele_2()?;
    println!("{allele_2:?}"); // Outputs ndarray ["A", "C", "C", "G"]
    Ok(())
}

#[test]
fn iid_index() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("some_missing.bed")?;
    let mut bed = Bed::new(file_name)?;

    // Read the individual at index position 3

    let val = ReadOptions::builder().iid_index(3).f64().read(&mut bed)?;
    assert!(val.dim() == (1, 100));

    // Read the individuals at index positions 0, 5, and 1st-from-last.

    let val = ReadOptions::builder()
        .iid_index([0, 5, -1])
        .f64()
        .read(&mut bed)?;

    assert!(val.dim() == (3, 100));

    // Read the individuals at index positions 20 (inclusive) to 30 (exclusive).

    let val = ReadOptions::builder()
        .iid_index(20..30)
        .f64()
        .read(&mut bed)?;

    assert!(val.dim() == (10, 100));

    // Read the individuals at every 2nd index position.

    let val = ReadOptions::builder()
        .iid_index(s![..;2])
        .f64()
        .read(&mut bed)?;

    assert!(val.dim() == (50, 100));

    // Read chromosome 5 of the the female individuals.

    let female = bed.sex()?.map(|elem| *elem == 2);
    let chrom_5 = bed.chromosome()?.map(|elem| elem == "5");
    let val = ReadOptions::builder()
        .iid_index(female)
        .sid_index(chrom_5)
        .f64()
        .read(&mut bed)?;

    println!("{:?}", val.dim());
    assert_eq!(val.dim(), (50, 6));

    Ok(())
}

#[test]
fn write_options_metadata() -> Result<(), Box<BedErrorPlus>> {
    // This test check interesting combinations of these options:
    // A: setting metadata (twice?)
    // B: giving val (early, late)
    // C: setting iid_count

    let output_folder = TempDir::default();
    let output_file = output_folder.join("small.bed");

    // <none>
    let write_options_result = WriteOptions::<f32>::builder(&output_file).build(0, 0);
    println!("{write_options_result:?}"); // Outputs fields of None
    write_options_result?;

    // A
    let write_options_result = WriteOptions::<f32>::builder(&output_file)
        .iid(["iid1", "iid2", "iid3"])
        .build(3, 4);
    println!("{write_options_result:?}"); // Outputs field of iid
    let iid_count = write_options_result?.iid_count();
    assert_eq!(iid_count, 3);

    // AA
    let write_options_result = WriteOptions::<f32>::builder(&output_file)
        .chromosome(["1", "1", "1"])
        .sid(["sid1", "sid2", "sid3", "sid4"])
        .build(3, 4);
    println!("{write_options_result:?}"); // Outputs field of iid
    assert_error_variant!(
        write_options_result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    // B early
    let val = nd::array![
        [1.0, 0.0, f64::NAN, 0.0],
        [2.0, 0.0, f64::NAN, 2.0],
        [0.0, 1.0, 2.0, 0.0]
    ];
    WriteOptions::builder(&output_file).write(&val)?;
    assert_eq!(3, Bed::new(&output_file)?.iid_count()?);

    // B late
    let write_options = WriteOptions::<f64>::builder(&output_file).build(3, 4)?;
    Bed::write_with_options(&val, &write_options)?;
    let sid_count = write_options.sid_count();
    println!("{sid_count:?}");
    assert_eq!(4, sid_count);

    // BB inconsistent
    let val2 = nd::array![[1.0, 0.0, f64::NAN, 0.0], [0.0, 1.0, 2.0, 0.0]];
    let write_options = WriteOptions::<f64>::builder(&output_file).build(3, 4)?;
    Bed::write_with_options(&val, &write_options)?;
    let result = Bed::write_with_options(&val2, &write_options);
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    // C
    let write_options_result = WriteOptions::<f32>::builder(&output_file).build(3, 4);
    println!("{write_options_result:?}"); // Outputs field of iid
    let sid_count = write_options_result?.sid_count();
    assert_eq!(sid_count, 4);

    // AC inconsistent
    let write_options_result = WriteOptions::<f32>::builder(&output_file)
        .iid(["iid1", "iid2", "iid3"])
        .build(4, 4);
    assert_error_variant!(
        write_options_result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    // AB early inconsistent
    let result = WriteOptions::builder(&output_file)
        .iid(["iid1", "iid2", "iid3", "iid4"])
        .write(&val);
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    // AB late inconsistent
    let write_options = WriteOptions::builder(&output_file)
        .iid(["iid1", "iid2", "iid3", "iid4"])
        .build(4, 4)?;
    let result = Bed::write_with_options(&val, &write_options);
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    // BC late inconsistent
    let write_options = WriteOptions::builder(&output_file).build(4, 4)?;
    let result = Bed::write_with_options(&val, &write_options);
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    // ABC early consistent
    WriteOptions::builder(&output_file)
        .sid(["sid1", "sid2", "sid3", "sid4"])
        .write(&val)?;

    // ABC late consistent
    let write_options = WriteOptions::builder(&output_file)
        .sid(["sid1", "sid2", "sid3", "sid4"])
        .build(3, 4)?;
    Bed::write_with_options(&val, &write_options)?;

    // Here Allele_1 and Allele_2 get their default values of A1,A1,... and A2,A2,...
    let write_options = WriteOptions::builder(output_file)
        .fid(["fid1", "fid1", "fid2"])
        .iid(["iid1", "iid2", "iid3"])
        .father(["iid23", "iid23", "iid22"])
        .mother(["iid34", "iid34", "iid33"])
        .sex([1, 2, 0])
        .pheno(["red", "red", "blue"])
        .chromosome(["1", "1", "5", "Y"])
        .sid(["sid1", "sid2", "sid3", "sid4"])
        .cm_position([100.4, 2000.5, 4000.7, 7000.9])
        .bp_position([1, 100, 1000, 1004])
        .f32()
        .build(3, 4)?;

    let metadata = write_options.metadata();
    println!("{metadata:?}");

    Ok(())
}

#[test]
fn metadata_use() -> Result<(), Box<BedErrorPlus>> {
    // Extract metadata from a file.
    // Create a random file with the same metadata.

    let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
    let metadata = bed.metadata()?;
    let shape = bed.dim()?;

    let mut rng = StdRng::seed_from_u64(0);
    let val = nd::Array::random_using(shape, Uniform::from(-1..3), &mut rng);

    let temp_out = TempDir::default();
    let output_file = temp_out.join("random.bed");
    WriteOptions::builder(output_file)
        .metadata(&metadata)
        .missing_value(-1)
        .write(&val)?;

    Ok(())
}

#[test]
fn metadata_same() -> Result<(), Box<BedErrorPlus>> {
    let iid_count = 1_000;
    let sid_count = 5_000;
    let file_count = 10;

    let metadata = Metadata::builder()
        .iid((0..iid_count).map(|iid_index| format!("iid_{iid_index}")))
        .sid((0..sid_count).map(|sid_index| format!("sid_{sid_index}")))
        .build()?;

    let temp_out = TempDir::default();

    for file_index in 0..file_count {
        let output_file = temp_out.join(format!("random{file_index}.bed"));

        let mut rng = StdRng::seed_from_u64(0);
        let val = nd::Array::random_using((iid_count, sid_count), Uniform::from(-1..3), &mut rng);
        WriteOptions::builder(output_file)
            .metadata(&metadata)
            .missing_value(-1)
            .write(&val)?;
    }
    Ok(())
}

// A - apply to reading
// B - extract from reading
// C - apply to writing
// D - extra from options when writing

// CD
// create 10 files with the same metadata
// structs: Metadata, Bed, ReadOptions, WriteOptions

#[test]
fn struct_play() -> Result<(), Box<BedErrorPlus>> {
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

#[test]
fn metadata_bed() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(file_name)?;
    let metadata = bed.metadata()?;
    println!("{0:?}", metadata.iid()); // Outputs Some(["iid1", "iid2", "iid3"] ...)
    println!("{0:?}", metadata.sid()); // Outputs Some(["sid1", "sid2", "sid3", "sid4"] ...)
    Ok(())
}

fn rt1<R>(range_thing: R) -> Result<Result<nd::Array2<i8>, Box<BedErrorPlus>>, Box<BedErrorPlus>>
where
    R: std::ops::RangeBounds<usize>
        + std::fmt::Debug
        + Clone
        + std::slice::SliceIndex<[isize], Output = [isize]>
        + std::panic::RefUnwindSafe,
{
    println!("Running {:?}", &range_thing);
    let file_name = sample_bed_file("toydata.5chrom.bed")?;

    let result1 = catch_unwind(|| {
        let mut bed = Bed::new(&file_name).unwrap();
        let all: Vec<isize> = (0..(bed.iid_count().unwrap() as isize)).collect();
        let mut bed = Bed::new(&file_name).unwrap();
        let iid_index: &[isize] = &all[range_thing.clone()];
        ReadOptions::builder()
            .iid_index(iid_index)
            .i8()
            .read(&mut bed)
    });
    match result1 {
        Err(_) => Err(BedError::PanickedThread().into()),
        Ok(bed_result) => Ok(bed_result),
    }
}

fn rt2(
    range_thing: bed_reader::Index,
) -> Result<Result<nd::Array2<i8>, Box<BedErrorPlus>>, Box<BedErrorPlus>> {
    println!("Running {:?}", &range_thing);
    let file_name = sample_bed_file("toydata.5chrom.bed")?;

    let result2 = catch_unwind(|| {
        let mut bed = Bed::new(file_name).unwrap();
        ReadOptions::builder()
            .iid_index(range_thing)
            .i8()
            .read(&mut bed)
    });
    match result2 {
        Err(_) => Err(BedError::PanickedThread().into()),
        Ok(bed_result) => Ok(bed_result),
    }
}

fn rt23(range_thing: bed_reader::Index) -> (RrArray2, RrUsize) {
    (rt2(range_thing.clone()), rt3(range_thing))
}

fn rt3(
    range_thing: bed_reader::Index,
) -> Result<Result<usize, Box<BedErrorPlus>>, Box<BedErrorPlus>> {
    println!("Running {:?}", &range_thing);
    let file_name = sample_bed_file("toydata.5chrom.bed").unwrap();

    let result3 = catch_unwind(|| {
        let mut bed = Bed::new(file_name).unwrap();
        range_thing.len(bed.iid_count().unwrap()).unwrap()
    });
    match result3 {
        Err(_) => Err(BedError::PanickedThread().into()),
        Ok(bed_result) => Ok(Ok(bed_result)),
    }
}

fn nds1(
    range_thing: SliceInfo1,
) -> Result<Result<nd::Array2<i8>, Box<BedErrorPlus>>, Box<BedErrorPlus>> {
    println!("Running {:?}", &range_thing);
    let file_name = sample_bed_file("toydata.5chrom.bed")?;

    let result1 = catch_unwind(|| {
        let mut bed = Bed::new(&file_name).unwrap();
        let all: nd::Array1<isize> = (0..(bed.iid_count().unwrap() as isize)).collect();
        let mut bed = Bed::new(&file_name).unwrap();
        let iid_index = &all.slice(&range_thing);
        ReadOptions::builder()
            .iid_index(iid_index)
            .i8()
            .read(&mut bed)
    });
    match result1 {
        Err(_) => Err(BedError::PanickedThread().into()),
        Ok(bed_result) => Ok(bed_result),
    }
}

type RrArray2 = Result<Result<nd::Array2<i8>, Box<BedErrorPlus>>, Box<BedErrorPlus>>;
type RrUsize = Result<Result<usize, Box<BedErrorPlus>>, Box<BedErrorPlus>>;

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

fn is_err2<T>(result_result: &Result<Result<T, Box<BedErrorPlus>>, Box<BedErrorPlus>>) -> bool {
    !matches!(result_result, Ok(Ok(_)))
}

#[test]
fn read_and_fill_with_options() -> Result<(), Box<BedErrorPlus>> {
    use bed_reader::assert_eq_nan;
    use bed_reader::{Bed, ReadOptions};
    use ndarray as nd;
    // Read the SNPs indexed by 2.
    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(file_name)?;
    let read_options = ReadOptions::builder().sid_index(2).build()?;
    let mut val = nd::Array2::<f64>::default((3, 1));
    bed.read_and_fill_with_options(&mut val.view_mut(), &read_options)?;

    assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);

    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(file_name)?;
    let mut val = nd::Array2::<f64>::default((3, 1));
    ReadOptions::builder()
        .sid_index(2)
        .read_and_fill(&mut bed, &mut val.view_mut())?;

    assert_eq_nan(&val, &nd::array![[f64::NAN], [f64::NAN], [2.0]]);

    Ok(())
}

#[test]
fn bed_builder_metadata() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("small.bed")?;

    // show that can fine errors
    let metadata = Metadata::builder()
        .iid(["i1", "i2", "i3"])
        .sid(["s1", "s2", "s3", "s4"])
        .build()?;
    let result = Bed::builder(&file_name)
        .fid(["f1", "f2", "f3", "f4"])
        .metadata(&metadata)
        .build();
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    let metadata = Metadata::builder()
        .iid(["i1", "i2", "i3"])
        .sid(["s1", "s2", "s3", "s4"])
        .build()?;
    let mut bed = Bed::builder(&file_name)
        .fid(["f1", "f2", "f3"])
        .iid(["x1", "x2", "x3"])
        .metadata(&metadata)
        .build()?;
    println!("{0:?}", bed.fid()?); // Outputs ndarray ["f1", "f2", "f3"]
    println!("{0:?}", bed.iid()?); // Outputs ndarray ["i1", "i2", "i3"]
    println!("{0:?}", bed.sid()?); // Outputs ndarray ["s1", "s2", "s3", "s4"]
    println!("{0:?}", bed.chromosome()?); // Outputs ndarray ["1", "1", "5", "Y"]

    let mut bed = Bed::builder(file_name)
        .skip_fid()
        .fid(["f1", "f2", "f3"])
        .iid(["x1", "x2", "x3"])
        .metadata(&metadata)
        .build()?;
    bed.fid().expect_err("Should fail");
    println!("{0:?}", bed.iid()?); // Outputs ndarray ["i1", "i2", "i3"]
    println!("{0:?}", bed.sid()?); // Outputs ndarray ["s1", "s2", "s3", "s4"]
    println!("{0:?}", bed.chromosome()?); // Outputs ndarray ["1", "1", "5", "Y"]

    Ok(())
}

#[test]
fn write_options_builder_metadata() -> Result<(), Box<BedErrorPlus>> {
    let output_folder = TempDir::default();

    let output_file = output_folder.join("small_m.bed");
    let metadata = Metadata::builder()
        .iid(["i1", "i2", "i3"])
        .sid(["s1", "s2", "s3", "s4"])
        .build()?;
    let write_options = WriteOptions::builder(output_file)
        .f32()
        .fid(["f1", "f2", "f3"])
        .iid(["x1", "x2", "x3"])
        .metadata(&metadata)
        .build(3, 4)?;
    println!("{0:?}", write_options.fid()); // Outputs ndarray ["f1", "f2", "f3"]
    println!("{0:?}", write_options.iid()); // Outputs ndarray ["i1", "i2", "i3"]
    println!("{0:?}", write_options.sid()); // Outputs ndarray ["s1", "s2", "s3", "s4"]
    println!("{0:?}", write_options.chromosome()); // Outputs ndarray ["0", "0", "0", "0"]
    Ok(())
}

#[test]
fn metadata_builder_metadata() -> Result<(), Box<BedErrorPlus>> {
    let metadata1 = Metadata::builder()
        .iid(["i1", "i2", "i3"])
        .sid(["s1", "s2", "s3", "s4"])
        .build()?;
    let metadata2 = Metadata::builder()
        .fid(["f1", "f2", "f3"])
        .iid(["x1", "x2", "x3"])
        .metadata(&metadata1)
        .build()?;
    println!("{0:?}", metadata2.fid()); // Outputs optional ndarray Some(["f1", "f2", "f3"]...)
    println!("{0:?}", metadata2.iid()); // Outputs optional ndarray Some(["i1", "i2", "i3"]...)
    println!("{0:?}", metadata2.sid()); // Outputs optional ndarray Some(["s1", "s2", "s3", "s4"]...)
    println!("{0:?}", metadata2.chromosome()); // Outputs None

    Ok(())
}
#[test]
fn metadata_builder() -> Result<(), Box<BedErrorPlus>> {
    let metadata = Metadata::builder()
        .iid(["i1", "i2", "i3"])
        .sid(["s1", "s2", "s3", "s4"])
        .build()?;
    let mut rng = StdRng::seed_from_u64(0);
    let val = nd::Array::random_using((3, 4), Uniform::from(-1..3), &mut rng);

    let temp_out = TempDir::default();
    let output_file = temp_out.join("random.bed");
    WriteOptions::builder(output_file)
        .metadata(&metadata)
        .missing_value(-1)
        .write(&val)?;

    Ok(())
}

#[test]
fn metadata_read_fam_bim() -> Result<(), Box<BedErrorPlus>> {
    let skip_set = HashSet::<MetadataFields>::new();
    let metadata_empty = Metadata::new();
    let (metadata_fam, iid_count) =
        metadata_empty.read_fam(sample_file("small.fam")?, &skip_set)?;
    let (metadata_bim, sid_count) = metadata_fam.read_bim(sample_file("small.bim")?, &skip_set)?;
    assert_eq!(iid_count, 3);
    assert_eq!(sid_count, 4);
    println!("{0:?}", metadata_bim.iid()); // Outputs optional ndarray Some(["iid1", "iid2", "iid3"]...)
    println!("{0:?}", metadata_bim.sid()); // Outputs optional ndarray Some(["sid1", "sid2", "sid3", "sid4"]...)
    println!("{0:?}", metadata_bim.chromosome()); // Outputs optional ndarray Some(["1", "1", "5", "Y"]...)

    Ok(())
}

#[test]
fn metadata_fill() -> Result<(), Box<BedErrorPlus>> {
    let metadata0 = Metadata::builder()
        .iid(["i1", "i2", "i3"])
        .sid(["s1", "s2", "s3", "s4"])
        .build()?;
    let metadata_filled = metadata0.fill(3, 4)?;
    println!("{0:?}", metadata_filled.iid()); // Outputs optional ndarray Some(["i1", "i2", "i3"]...)
    println!("{0:?}", metadata_filled.sid()); // Outputs optional ndarray Some(["s1", "s2", "s3", "s4"]...)
    println!("{0:?}", metadata_filled.chromosome()); // Outputs optional ndarray Some(["0", "0", "0", "0"]...)

    Ok(())
}

#[test]
fn metadata_fam_bim_write() -> Result<(), Box<BedErrorPlus>> {
    let metadata0 = Metadata::builder()
        .iid(["i1", "i2", "i3"])
        .sid(["s1", "s2", "s3", "s4"])
        .build()?;
    let metadata_filled = metadata0.fill(3, 4)?;

    let temp_out = TempDir::default();
    let output_file = temp_out.join("no_bed.fam");
    metadata_filled.write_fam(output_file)?;

    let temp_out = TempDir::default();
    let output_file = temp_out.join("no_bed.bim");
    metadata_filled.write_bim(output_file)?;

    Ok(())
}

#[test]
fn metadata_iid_sid() -> Result<(), Box<BedErrorPlus>> {
    let metadata = Metadata::builder().iid(["i1", "i2", "i3"]).build()?;
    println!("{0:?}", metadata.iid()); // Outputs optional ndarray Some(["i1", "i2", "i3"]...)
    println!("{0:?}", metadata.sid()); // Outputs None
    Ok(())
}

#[test]
fn read_options_properties() -> Result<(), Box<BedErrorPlus>> {
    let read_options = ReadOptions::builder().sid_index([2, 3, 0]).i8().build()?;
    assert_eq!(read_options.missing_value(), -127);
    println!("{0:?}", read_options.iid_index()); // Outputs 'All'
    println!("{0:?}", read_options.sid_index()); // Outputs 'Vec([2, 3, 0])'
    assert!(read_options.is_f());
    assert!(read_options.is_a1_counted());
    assert_eq!(read_options.num_threads(), None);

    let file_name = sample_bed_file("small.bed")?;
    let mut bed = Bed::new(file_name)?;
    let val = bed.read_with_options(&read_options)?;

    assert_eq_nan(&val, &nd::array![[-127, 0, 1], [-127, 2, 2], [2, 0, 0]]);

    Ok(())
}

#[test]
fn write_options_properties() -> Result<(), Box<BedErrorPlus>> {
    let output_folder = TempDir::default();
    let output_file = output_folder.join("small.bed");
    let write_options = WriteOptions::builder(output_file)
        .f64()
        .iid(["i1", "i2", "i3"])
        .sid(["s1", "s2", "s3", "s4"])
        .build(3, 4)?;

    println!("{0:?}", write_options.path()); // Outputs "...small.bed"
    println!("{0:?}", write_options.fam_path()); // Outputs "...small.fam"
    println!("{0:?}", write_options.bim_path()); // Outputs "...small.bim"

    println!("{0:?}", write_options.fid()); // Outputs ndarray ["0", "0", "0"]
    println!("{0:?}", write_options.iid()); // Outputs ndarray ["i1", "i2", "i3"]
    println!("{0:?}", write_options.father()); // Outputs ndarray ["0", "0", "0"]
    println!("{0:?}", write_options.mother()); // Outputs ndarray ["0", "0", "0"]
    println!("{0:?}", write_options.sex()); // Outputs ndarray [0, 0, 0]
    println!("{0:?}", write_options.pheno()); // Outputs ndarray ["0", "0", "0"]

    println!("{0:?}", write_options.chromosome()); // Outputs ndarray ["0", "0", "0", "0"]
    println!("{0:?}", write_options.sid()); // Outputs ndarray ["s1", "s2", "s3", "s4"]
    println!("{0:?}", write_options.cm_position()); // Outputs ndarray [0.0, 0.0, 0.0, 0.0]
    println!("{0:?}", write_options.bp_position()); // Outputs ndarray [0, 0, 0, 0]
    println!("{0:?}", write_options.allele_1()); // Outputs ndarray ["A1", "A1", "A1", "A1"]
    println!("{0:?}", write_options.allele_2()); // Outputs ndarray ["A2", "A2", "A2", "A2"]

    let metadata = write_options.metadata();
    println!("{0:?}", metadata.iid()); // Outputs optional ndarray Some(["i1", "i2", "i3"])

    assert_eq!(write_options.iid_count(), 3);
    assert_eq!(write_options.sid_count(), 4);
    assert_eq!(write_options.dim(), (3, 4));

    let val = nd::array![
        [1.0, 0.0, f64::NAN, 0.0],
        [2.0, 0.0, f64::NAN, 2.0],
        [0.0, 1.0, 2.0, 0.0]
    ];

    Bed::write_with_options(&val, &write_options)?;

    Ok(())
}

#[test]
fn write_options_builder() -> Result<(), Box<BedErrorPlus>> {
    let output_folder = TempDir::default();
    let output_file = output_folder.join("small.bed");
    let val = nd::array![[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]];

    WriteOptions::builder(output_file)
        .num_threads(1)
        .write(&val)?;

    let output_folder = TempDir::default();
    let output_file = output_folder.join("small.deb");
    let val = nd::array![[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]];

    WriteOptions::builder(output_file)
        .fam_path(output_folder.join("small.maf"))
        .bim_path(output_folder.join("small.mib"))
        .write(&val)?;
    Ok(())
}

// Bed, file set iid and fid and metadata and iid_count, and check if error is caught
#[test]
fn bed_inconsistent_count() -> Result<(), Box<BedErrorPlus>> {
    // Bed: fid vs metadata
    let metadata = Metadata::builder().iid(["f1", "f2", "f3", "f4"]).build()?;
    let result = Bed::builder(sample_bed_file("small.bed")?)
        .fid(["f1", "f2", "f3"])
        .metadata(&metadata)
        .build();
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    // Bed: file vs file:
    let mut bed = Bed::builder(sample_bed_file("small.bed")?)
        .fam_path(sample_file("small.fam_bad")?)
        .skip_iid()
        .skip_father()
        .skip_mother()
        .skip_sex()
        .skip_pheno()
        .build()?;
    let result = bed.fid();
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::MetadataFieldCount(_, _, _))
    );

    // Bed: fid vs iid
    let result = Bed::builder(sample_bed_file("small.bed")?)
        .fid(["f1", "f2", "f3"])
        .iid(["i1", "i2", "i3", "i4"])
        .build();
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    // Bed: iid vs file
    let mut bed = Bed::builder(sample_bed_file("small.bed")?)
        .iid(["i1", "i2", "i3", "i4"])
        .build()?;
    let result = bed.fid();
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    Ok(())
}

// WriteOptions, file set iid and fid and metadata and iid_count, and check if error is caught
#[test]
fn write_options_inconsistent_count() -> Result<(), Box<BedErrorPlus>> {
    let temp_out = TempDir::default();
    let output_file = temp_out.join("small.bed");

    // WriteOptions: fid vs build
    let metadata = Metadata::builder().build()?;
    let result = WriteOptions::builder(&output_file)
        .i8()
        .fid(["f1", "f2", "f3", "f4"])
        .metadata(&metadata)
        .build(3, 4);
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    // WriteOptions: metadata vs build
    let metadata = Metadata::builder().iid(["f1", "f2", "f3", "f4"]).build()?;
    let result = WriteOptions::builder(&output_file)
        .i8()
        .metadata(&metadata)
        .build(3, 4);
    assert_error_variant!(
        result,
        BedErrorPlus::BedError(BedError::InconsistentCount(_, _, _))
    );

    Ok(())
}

// MetadataBuilders
#[test]
fn metadata_inconsistent_count() -> Result<(), Box<BedErrorPlus>> {
    let skip_set = HashSet::<MetadataFields>::new();

    // Metadata: iid vs file
    let metadata = Metadata::builder().iid(["i1", "i2", "i3", "i4"]).build()?;
    let result = metadata.read_fam(sample_file("small.fam")?, &skip_set);
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
    let result = metadata.read_fam(sample_file("small.fam_bad")?, &skip_set);
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

#[test]
fn write_options_examples() -> Result<(), Box<BedErrorPlus>> {
    let output_folder = TempDir::default();
    let output_file = output_folder.join("small.bed");
    let write_options = WriteOptions::builder(output_file)
        .i8()
        .iid(["i1", "i2", "i3"])
        .sid(["s1", "s2", "s3", "s4"])
        .build(3, 4)?;
    assert!(write_options.is_a1_counted());
    assert!(write_options.num_threads().is_none());
    assert!(write_options.missing_value() == -127);

    Ok(())
}

#[test]
fn parsing_metadata() -> Result<(), Box<BedErrorPlus>> {
    let bed_fam_bim = sample_files(["small.bed", "small.fam", "small.bim_bad_positions.bim"])?;
    let mut bed = Bed::builder(&bed_fam_bim[0])
        .bim_path(&bed_fam_bim[2])
        .build()?;
    let result = bed.cm_position();
    assert_error_variant!(result, BedErrorPlus::ParseIntError(_));
    Ok(())
}

#[test]
fn read_fam() -> Result<(), Box<BedErrorPlus>> {
    let skip_set = HashSet::<MetadataFields>::new();
    let metadata_empty = Metadata::new();
    let (metadata_fam, _) = metadata_empty.read_fam(sample_file("small.fam")?, &skip_set)?;
    // metadata_empty.read_fam("bed_reader/tests/data/small.fam", &skip_set)?;
    println!("{:?}", metadata_fam.iid()); // Outputs optional ndarray Some(["iid1", "iid2", "iid3"]...)
    Ok(())
}

#[test]
fn lib_intro() -> Result<(), Box<BedErrorPlus>> {
    let file_name = sample_bed_file("some_missing.bed")?;

    let mut bed = Bed::new(file_name)?;
    println!("{:?}", bed.iid()?.slice(s![..5])); // Outputs ndarray: ["iid_0", "iid_1", "iid_2", "iid_3", "iid_4"]
    println!("{:?}", bed.sid()?.slice(s![..5])); // Outputs ndarray: ["sid_0", "sid_1", "sid_2", "sid_3", "sid_4"]
    println!("{:?}", bed.chromosome()?.iter().collect::<HashSet<_>>());
    // Outputs: {"12", "10", "4", "8", "19", "21", "9", "15", "6", "16", "13", "7", "17", "18", "1", "22", "11", "2", "20", "3", "5", "14"}
    let _ = ReadOptions::builder()
        .sid_index(bed.chromosome()?.map(|elem| elem == "5"))
        .f64()
        .read(&mut bed)?;

    Ok(())
}

#[test]
fn metadata_iid_sid2() -> Result<(), Box<BedErrorPlus>> {
    let metadata = Metadata::builder()
        .iid(["sample1", "sample2", "sample3"])
        .build()?;
    println!("{:?}", metadata.iid()); // Outputs ndarray Some(["sample1", "sample2", "sample3"])

    let metadata = Metadata::builder()
        .sid(["SNP1", "SNP2", "SNP3", "SNP4"])
        .build()?;
    println!("{:?}", metadata.sid()); // Outputs ndarray Some(["SNP1", "SNP2", "SNP3", "SNP4"])

    Ok(())
}
