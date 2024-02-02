#[cfg(feature = "tokio")]
use bed_reader::{BedCloud, BedErrorPlus, ReadOptions};
#[cfg(feature = "tokio")]
use ndarray::s;
#[cfg(feature = "tokio")]
use std::collections::HashSet;

#[cfg(feature = "tokio")]
#[tokio::main]
async fn main() -> Result<(), Box<BedErrorPlus>> {
    let url = "https://raw.githubusercontent.com/fastlmm/bed-sample-files/main/some_missing.bed";
    let cloud_options = [("timeout", "10s")];

    let mut bed_cloud = BedCloud::new_with_options(url, cloud_options).await?;
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

#[cfg(not(feature = "tokio"))]
pub fn main() {
    println!("This example requires the 'tokio/full' feature to be enabled");
}
