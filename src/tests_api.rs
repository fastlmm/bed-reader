#[cfg(test)]
use crate::BedError;
#[cfg(test)]
use crate::BedErrorPlus;
// !!!cmk later use read_all or new macros to make reading all easier.
// !!!cmk later is there a way to set default value based on the result type (if given)
#[cfg(test)]
use crate::api::BedBuilder;
#[cfg(test)]
use crate::api::ReadArgBuilder;
#[cfg(test)]
use ndarray as nd;

#[test]
fn rusty_bed1() {
    let file = "bed_reader/tests/data/plink_sim_10s_100v_10pmiss.bed";
    let bed = BedBuilder::default()
        .filename(file.to_string())
        .build()
        .unwrap();
    let val = bed
        .read(
            ReadArgBuilder::default()
                .missing_value(-127)
                .build()
                .unwrap(),
        )
        .unwrap();
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    assert!(mean == -13.142); // really shouldn't do mean on data where -127 represents missing

    let bed = BedBuilder::default()
        .filename(file.to_string())
        .count_a1(false)
        .build()
        .unwrap();
    let val = bed
        .read(
            ReadArgBuilder::default()
                .missing_value(-127)
                .build()
                .unwrap(),
        )
        .unwrap();
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    assert!(mean == -13.274); // really shouldn't do mean on data where -127 represents missing
}

#[test]
fn rusty_bed2() {
    // !!!cmk later reading one iid is very common. Make it easy.
    let file = "bed_reader/tests/data/plink_sim_10s_100v_10pmiss.bed";
    let bed = BedBuilder::default()
        .filename(file.to_string())
        .build()
        .unwrap();

    let val = bed
        .read(
            ReadArgBuilder::default()
                .missing_value(-127)
                // !!!cmk 0 could it be any slice of usize?
                // !!!cmk 0 or work with everything that works with ndarray indexing
                .iid_index([0].to_vec().into())
                .sid_index([1].to_vec().into())
                .build()
                .unwrap(),
        )
        .unwrap();
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{:?}", mean);
    assert!(mean == 1.0); // really shouldn't do mean on data where -127 represents missing
}

// !!!cmk later ask reddit help (mention builder library creator)
// macro_rules! read {
//     ($b:expr,$x:expr) => {
//         $b.read(ReadArgBuilder::default().$x.build())
//     }; // () => {
//        //      ReadArgBuilder::default().build()
//        //};
// }

#[cfg(test)]
use std::{collections::HashSet, f64::NAN};

#[test]
fn rusty_bed3() {
    // !!!cmk later also show mixing bool and full and none
    let file = "bed_reader/tests/data/plink_sim_10s_100v_10pmiss.bed";
    let mut bed = BedBuilder::default()
        .filename(file.to_string())
        .build()
        .unwrap();
    let iid_bool: nd::Array1<bool> = (0..bed.get_iid_count())
        .map(|elem| (elem % 2) != 0)
        .collect();
    let sid_bool: nd::Array1<bool> = (0..bed.get_sid_count())
        .map(|elem| (elem % 8) != 0)
        .collect();
    let val = bed
        .read(
            ReadArgBuilder::default()
                .missing_value(-127)
                .iid_index(iid_bool.into())
                .sid_index(sid_bool.into())
                .build()
                .unwrap(),
        )
        .unwrap();
    let mean = val.mapv(|elem| elem as f64).mean().unwrap();
    println!("{:?}", mean);
    assert!(mean == -14.50344827586207); // really shouldn't do mean on data where -127 represents missing
}

#[test]
fn bad_header() {
    let filename = "bed_reader/tests/data/badfile.bed";
    let result = BedBuilder::default().filename(filename.to_string()).build();

    match result {
        Err(BedErrorPlus::BedError(BedError::IllFormed(_))) => (),
        _ => panic!("test failure"),
    };
}

#[test]
fn readme_examples() {
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
    // !!!cmk later pull down sample file
    let file_name = "bed_reader/tests/data/small.bed";
    let bed = BedBuilder::default()
        .filename(file_name.to_string())
        .build()
        .unwrap();
    //nd let bed = Bed::new(filename)?;
    let val = bed
        .read(
            ReadArgBuilder::default()
                .missing_value(NAN)
                .build()
                .unwrap(),
        )
        .unwrap();
    //nd val = bed.read(NAN).unwrap();
    //nd val = bed.read!().unwrap();
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
    let bed2 = BedBuilder::default()
        .filename(file_name2.to_string())
        .build()
        .unwrap();
    let val2 = bed2
        .read(
            ReadArgBuilder::default()
                .iid_index(nd::s![..;2].into())
                .sid_index(nd::s![20..30].into())
                .missing_value(NAN)
                .build()
                .unwrap(),
        )
        .unwrap();
    //nd val2 = bed2.read!(index(s![..,20..30]).missing_value(NAN)).unwrap();
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

    let mut bed3 = BedBuilder::default()
        .filename(file_name2.to_string())
        .build()
        .unwrap();
    println!("{:?}", bed3.get_iid().slice(nd::s![5..]));
    println!("{:?}", bed3.get_sid().slice(nd::s![5..]));
    let unique: HashSet<_> = bed3.get_chromosome().iter().collect();
    println!("{:?}", unique);
    let is_5 = bed3.get_chromosome().mapv(|elem| elem == "5");
    let val3 = bed3
        .read(
            ReadArgBuilder::default()
                .sid_index(is_5.into())
                .missing_value(NAN)
                .build()
                .unwrap(),
        )
        .unwrap();
    //nd val3 = bed3.read!(sid_index(is_5).missing_value(NAN)).unwrap();
    println!("{:?}", val3.shape());
    // [100, 6]
}
