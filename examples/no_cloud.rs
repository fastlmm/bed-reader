use bed_reader::{sample_bed_file, Bed, BedErrorPlus, ReadOptions};
use ndarray::s;
use std::collections::HashSet;

fn main() -> Result<(), Box<BedErrorPlus>> {
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
