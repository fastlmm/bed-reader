#[cfg(test)]
use crate::api::write;
#[cfg(test)]
use crate::api::Bed;
#[cfg(test)]
use crate::api::ReadOptions;
#[cfg(test)]
use crate::api::WriteOptions;
#[cfg(test)]
use crate::BedError;
#[cfg(test)]
use crate::BedErrorPlus;
// !!!cmk later use read_all or new macros to make reading all easier.
// !!!cmk later is there a way to set default value based on the result type (if given)
#[cfg(test)]
use crate::tests::allclose;
#[cfg(test)]
use ndarray as nd;
#[cfg(test)]
use ndarray::s;
#[cfg(test)]
use temp_testdir::TempDir;

#[test]
fn rusty_bed1() -> Result<(), BedErrorPlus> {
    let file = "bed_reader/tests/data/plink_sim_10s_100v_10pmiss.bed";
    let mut bed = Bed::new(file)?;
    let val = bed.read::<i8>()?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    assert!(mean == -13.142); // really shouldn't do mean on data where -127 represents missing

    let mut bed = Bed::builder(file).count_a2().build()?;
    let val = bed.read::<i8>()?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    assert!(mean == -13.274); // really shouldn't do mean on data where -127 represents missing
    Ok(())
}

#[test]
fn rusty_bed2() -> Result<(), BedErrorPlus> {
    // !!!cmk later reading one iid is very common. Make it easy.
    let file = "bed_reader/tests/data/plink_sim_10s_100v_10pmiss.bed";
    let mut bed = Bed::new(file)?;

    let val = ReadOptions::builder()
        .iid_index(0)
        .sid_index(vec![1])
        .i8()
        .read(&mut bed)?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{:?}", mean);
    assert!(mean == 1.0); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

// !!!cmk later ask reddit help (mention builder library creator)
// macro_rules! read {
//     ($b:expr,$x:expr) => {
//         $b.read(ReadOptionsBuilder::default().$x.build())
//     }; // () => {
//        //      ReadOptionsBuilder::default().build()
//        //};
// }

#[cfg(test)]
use std::collections::HashSet;
#[cfg(test)]
use std::fs;
#[cfg(test)]
use std::path::PathBuf;

#[test]
fn rusty_bed3() -> Result<(), BedErrorPlus> {
    // !!!cmk later also show mixing bool and full and none
    let file = "bed_reader/tests/data/plink_sim_10s_100v_10pmiss.bed";
    let mut bed = Bed::new(file)?;
    let iid_bool: nd::Array1<bool> = (0..bed.iid_count()?).map(|elem| (elem % 2) != 0).collect();
    let sid_bool: nd::Array1<bool> = (0..bed.sid_count()?).map(|elem| (elem % 8) != 0).collect();
    let val = ReadOptions::builder()
        .missing_value(-127)
        .iid_index(iid_bool)
        .sid_index(sid_bool)
        .read(&mut bed)?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{:?}", mean);
    assert!(mean == -14.50344827586207); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[test]
fn rusty_bed_allele() -> Result<(), BedErrorPlus> {
    let file = "bed_reader/tests/data/plink_sim_10s_100v_10pmiss.bed";
    let mut bed = Bed::builder(file).count_a2().build()?;
    let val: nd::Array2<i8> = bed.read()?;

    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{:?}", mean);
    assert!(mean == -13.274); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[test]
fn rusty_bed_order() -> Result<(), BedErrorPlus> {
    let file = "bed_reader/tests/data/plink_sim_10s_100v_10pmiss.bed";
    let mut bed = Bed::new(file)?;
    let val = ReadOptions::builder().c().i8().read(&mut bed)?;

    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{:?}", mean);
    assert!(mean == -13.142); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[test]
fn bad_header() -> Result<(), BedErrorPlus> {
    let filename = "bed_reader/tests/data/badfile.bed";

    let _ = Bed::builder(filename).skip_early_check().build()?;

    let result = Bed::new(filename);

    match result {
        Err(BedErrorPlus::BedError(BedError::IllFormed(_))) => (),
        _ => panic!("test failure"),
    };

    Ok(())
}

#[test]
fn open_examples() -> Result<(), BedErrorPlus> {
    //     >>> from bed_reader import open_bed, sample_file
    //     >>>
    //     >>> file_name = sample_file("small.bed")
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

    let file_name = "bed_reader/tests/data/small.bed";
    let mut bed = Bed::new(file_name)?;
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

    let mut bed = Bed::new(file_name)?;
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

    let mut bed = Bed::builder(file_name)
        .iid(["sample1", "sample2", "sample3"])
        .build()?;
    println!("{:?}", bed.iid()?);
    println!("{:?}", bed.sid()?);

    // ["sample1", "sample2", "sample3"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1
    // ["sid1", "sid2", "sid3", "sid4"], shape=[4], strides=[1], layout=CFcf (0xf), const ndim=

    // Do more testing of Rust
    let iid = nd::array!["sample1", "sample2", "sample3"];
    let mut _bed = Bed::builder(file_name).iid(iid).build()?;
    let iid = nd::array![
        "sample1".to_string(),
        "sample2".to_string(),
        "sample3".to_string()
    ];
    let mut _bed = Bed::builder(file_name).iid(iid).build()?;
    let iid = vec!["sample1", "sample2", "sample3"];
    let mut _bed = Bed::builder(file_name).iid(iid).build()?;
    let iid = vec![
        "sample1".to_string(),
        "sample2".to_string(),
        "sample3".to_string(),
    ];
    let mut _bed = Bed::builder(file_name).iid(iid).build()?;

    // Give the number of individuals (samples) and SNPs (variants) so that the .fam and
    // .bim files need never be opened.

    //     >>> with open_bed(file_name, iid_count=3, sid_count=4) as bed:
    //     ...     print(bed.read())
    //     [[ 1.  0. nan  0.]
    //      [ 2.  0. nan  2.]
    //      [ 0.  1.  2.  0.]]

    let mut bed = Bed::builder(file_name).iid_count(3).sid_count(4).build()?;
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

    // !!! cmk later document that if you skip and then give default value, its the last that matters.
    // !!! cmk later test that sid_count/iid_count will raise an error if any metadata gives a different count
    let mut bed = Bed::builder(file_name).skip_allele_2().build()?;
    println!("{:?}", bed.iid()?);

    let result = bed.allele_2();
    match result {
        Err(BedErrorPlus::BedError(BedError::CannotUseSkippedMetadata(_))) => (),
        _ => panic!("test failure"),
    };

    Ok(())
}

#[test]
fn metadata_etc() -> Result<(), BedErrorPlus> {
    // >>> file_name = sample_file("small.bed")
    // >>> with open_bed(file_name) as bed:
    // ...     print(bed.sex)
    // [1 2 0]
    let mut bed = Bed::new("bed_reader/tests/data/small.bed")?;
    println!("{:?}", bed.sex()?);
    // [1, 2, 0], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1

    println!("{:?}", bed.cm_position()?);
    // [100.4, 2000.5, 4000.7, 7000.9], shape=[4], strides=[1], layout=CFcf (0xf), const ndim=1

    println!("{:?}", bed.bp_position()?);
    // [1, 100, 1000, 1004], shape=[4], strides=[1], layout=CFcf (0xf), const ndim=1

    println!("{:?}", bed.fid()?);
    // ["fid1", "fid1", "fid2"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1

    println!("{:?}", bed.father()?);
    // ["iid23", "iid23", "iid22"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1

    println!("{:?}", bed.mother()?);
    // ["iid34", "iid34", "iid33"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1

    Ok(())
}

#[test]
fn hello_father() -> Result<(), BedErrorPlus> {
    let mut bed = Bed::builder("bed_reader/tests/data/small.bed")
        .father(["f1", "f2", "f3"])
        .skip_mother()
        .build()?;
    println!("{:?}", bed.father()?);
    // ["f1", "f2", "f3"], shape=[3], strides=[1], layout=CFcf (0xf), const ndim=1
    bed.mother().unwrap_err();

    Ok(())
}

#[test]
fn num_threads() -> Result<(), BedErrorPlus> {
    let file = "bed_reader/tests/data/plink_sim_10s_100v_10pmiss.bed";
    let mut bed = Bed::new(file)?;

    let val = ReadOptions::builder().num_threads(4).i8().read(&mut bed)?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{:?}", mean);
    assert!(mean == -13.142); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[test]
fn fam_and_bim() -> Result<(), BedErrorPlus> {
    let mut bed = Bed::builder("bed_reader/tests/data/small.deb")
        .fam_path("bed_reader/tests/data/small.maf")
        .bim_path("bed_reader/tests/data/small.mib")
        .build()?;

    let val: nd::Array2<i8> = bed.read()?;
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{:?}", mean);
    assert!(mean == -20.5); // really shouldn't do mean on data where -127 represents missing

    Ok(())
}

#[test]
fn readme_examples() -> Result<(), BedErrorPlus> {
    // Read genomic data from a .bed file.

    // >>> import numpy as np
    // >>> from bed_reader import open_bed, sample_file
    // >>>
    // >>> file_name = sample_file("small.bed")
    // >>> bed = open_bed(file_name)
    // >>> val = bed.read()
    // >>> print(val)
    // [[ 1.  0. nan  0.]
    //  [ 2.  0. nan  2.]
    //  [ 0.  1.  2.  0.]]
    // >>> del bed

    // !!!cmk later document use statements
    // !!!cmk ask is there a rust crate for pulling down files if needed (using hash to check if file correct), like Python's Pooch
    let file_name = "bed_reader/tests/data/small.bed";
    let mut bed = Bed::new(file_name)?;
    let val = bed.read::<f64>()?;
    println!("{:?}", val);
    // [[1.0, 0.0, NaN, 0.0],
    // [2.0, 0.0, NaN, 2.0],
    // [0.0, 1.0, 2.0, 0.0]], shape=[3, 4], strides=[1, 3], layout=Ff (0xa), const ndim=2

    // Read every second individual and SNPs (variants) from 20 to 30.

    // >>> file_name2 = sample_file("some_missing.bed")
    // >>> bed2 = open_bed(file_name2)
    // >>> val2 = bed2.read(index=np.s_[::2,20:30])
    // >>> print(val2.shape)
    // (50, 10)
    // >>> del bed2

    let file_name2 = "bed_reader/tests/data/some_missing.bed";
    let mut bed2 = Bed::new(file_name2)?;
    let val2 = ReadOptions::<f64>::builder()
        .iid_index(s![..;2])
        .sid_index(20..30)
        .read(&mut bed2)?;
    println!("{:?}", val2.shape());
    // [50, 10]

    // List the first 5 individual (sample) ids, the first 5 SNP (variant) ids, and every unique chromosome. Then, read every value in chromosome 5.

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

    let mut bed3 = Bed::new(file_name2)?;
    println!("{:?}", bed3.iid()?.slice(s![..5]));
    println!("{:?}", bed3.sid()?.slice(s![..5]));
    let unique = bed3.chromosome()?.iter().collect::<HashSet<_>>();
    println!("{:?}", unique);
    // let is_5 = bed3.get_chromosome()?.map(|elem| elem == "5");
    let is_5 = nd::Zip::from(bed3.chromosome()?).par_map_collect(|elem| elem == "5");
    let val3 = ReadOptions::builder()
        .sid_index(is_5)
        .f64()
        .read(&mut bed3)?;
    println!("{:?}", val3.shape());
    // ["iid_0", "iid_1", "iid_2", "iid_3", "iid_4"], shape=[5], strides=[1], layout=CFcf (0xf), const ndim=1
    // ["sid_0", "sid_1", "sid_2", "sid_3", "sid_4"], shape=[5], strides=[1], layout=CFcf (0xf), const ndim=1
    // {"10", "11", "4", "21", "22", "14", "3", "12", "20", "15", "19", "8", "6", "18", "9", "2", "16", "13", "17", "1", "7", "5"}
    // [100, 6]
    Ok(())
}

// !!!cmk later re-write python_module.rs to use the new Rust API (may need .fill() and .fill_with_defaults())
// !!!cmk 0 more options
// .count_a2()
// .missing_value(f64::NAN)
// num_threads=None,

#[test]
fn write_docs() -> Result<(), BedErrorPlus> {
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

    // !!!cmk 0 make this a function
    let output_folder = tmp_path()?;

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
    write(&val, &output_file2)?;
    let mut bed2 = Bed::new(&output_file2)?;
    // println!("{:?}", bed2.chromosome()?);
    println!("{:?}", bed2.sid()?);
    // ['0' '0' '0' '0']

    Ok(())
}

#[cfg(test)]
fn tmp_path() -> Result<PathBuf, BedErrorPlus> {
    let output_path = PathBuf::from(TempDir::default().as_ref());
    fs::create_dir(&output_path)?;
    Ok(output_path)
}

#[test]
fn read_write() -> Result<(), BedErrorPlus> {
    // with open_bed(shared_datadir / "small.bed") as bed:
    // val = bed.read()
    // properties = bed.properties

    let file_name = "bed_reader/tests/data/small.bed";
    let mut bed = Bed::new(file_name)?;
    let val = bed.read::<f64>()?;
    let metadata = bed.metadata()?;
    println!("{:?}", metadata);

    // output_file = tmp_path / "small.deb"
    // fam_file = tmp_path / "small.maf"
    // bim_file = tmp_path / "small.mib"

    let temp_out = tmp_path()?;
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

    // assert np.allclose(val, val2, equal_nan=True)
    assert!(
        allclose(&val.view(), &val2.view(), 1e-08, true),
        "not close"
    );
    println!("{:?}", metadata);
    println!("{:?}", metadata2);
    assert!(metadata == metadata2, "meta not equal");

    Ok(())
}
