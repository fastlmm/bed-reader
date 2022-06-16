bed-reader
==========

[<img alt="github" src="https://img.shields.io/badge/github-bed--reader-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/fastlmm/bed-reader)
[<img alt="crates.io" src="https://img.shields.io/crates/v/bed-reader.svg?style=for-the-badge&color=fc8d62&logo=rust" height="20">](https://crates.io/crates/bed-reader)
[<img alt="docs.rs" src="https://img.shields.io/badge/docs.rs-bed--reader-66c2a5?style=for-the-badge&labelColor=555555&logoColor=white&logo=data:image/svg+xml;base64,PHN2ZyByb2xlPSJpbWciIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgdmlld0JveD0iMCAwIDUxMiA1MTIiPjxwYXRoIGZpbGw9IiNmNWY1ZjUiIGQ9Ik00ODguNiAyNTAuMkwzOTIgMjE0VjEwNS41YzAtMTUtOS4zLTI4LjQtMjMuNC0zMy43bC0xMDAtMzcuNWMtOC4xLTMuMS0xNy4xLTMuMS0yNS4zIDBsLTEwMCAzNy41Yy0xNC4xIDUuMy0yMy40IDE4LjctMjMuNCAzMy43VjIxNGwtOTYuNiAzNi4yQzkuMyAyNTUuNSAwIDI2OC45IDAgMjgzLjlWMzk0YzAgMTMuNiA3LjcgMjYuMSAxOS45IDMyLjJsMTAwIDUwYzEwLjEgNS4xIDIyLjEgNS4xIDMyLjIgMGwxMDMuOS01MiAxMDMuOSA1MmMxMC4xIDUuMSAyMi4xIDUuMSAzMi4yIDBsMTAwLTUwYzEyLjItNi4xIDE5LjktMTguNiAxOS45LTMyLjJWMjgzLjljMC0xNS05LjMtMjguNC0yMy40LTMzLjd6TTM1OCAyMTQuOGwtODUgMzEuOXYtNjguMmw4NS0zN3Y3My4zek0xNTQgMTA0LjFsMTAyLTM4LjIgMTAyIDM4LjJ2LjZsLTEwMiA0MS40LTEwMi00MS40di0uNnptODQgMjkxLjFsLTg1IDQyLjV2LTc5LjFsODUtMzguOHY3NS40em0wLTExMmwtMTAyIDQxLjQtMTAyLTQxLjR2LS42bDEwMi0zOC4yIDEwMiAzOC4ydi42em0yNDAgMTEybC04NSA0Mi41di03OS4xbDg1LTM4Ljh2NzUuNHptMC0xMTJsLTEwMiA0MS40LTEwMi00MS40di0uNmwxMDItMzguMiAxMDIgMzguMnYuNnoiPjwvcGF0aD48L3N2Zz4K" height="20">](https://docs.rs/bed-reader)
[<img alt="build status" src="https://img.shields.io/github/workflow/status/fastlmm/bed-reader/CI/master?style=for-the-badge" height="20">](https://github.com/fastlmm/bed-reader/actions?query=branch%3Amaster)

Read and write the PLINK BED format, simply and efficiently.

Features
--------

* Fast and multi-threaded
* Supports many indexing methods. Slice data by individuals (samples) and/or SNPs (variants).
* The Python-facing APIs for this library is used by [PySnpTools](https://github.com/fastlmm/PySnpTools), [FaST-LMM](https://github.com/fastlmm/FaST-LMM), and [PyStatGen](https://github.com/pystatgen).
* Supports [PLINK 1.9](https://www.cog-genomics.org/plink2/formats).

Examples
--------

Read all genotype data from a .bed file.

```rust
use ndarray as nd;
use bed_reader::{Bed, ReadOptions, assert_eq_nan, sample_bed_file};

let file_name = sample_bed_file("small.bed")?;
let mut bed = Bed::new(file_name)?;
let val = ReadOptions::builder().f64().read(&mut bed)?;

assert_eq_nan(
    &val,
    &nd::array![
        [1.0, 0.0, f64::NAN, 0.0],
        [2.0, 0.0, f64::NAN, 2.0],
        [0.0, 1.0, 2.0, 0.0]
    ],
);
# use bed_reader::BedErrorPlus; // '#' needed for doctest
# Ok::<(), BedErrorPlus>(())
```

Read every second individual (samples) and SNPs (variants) 20 to 30.

```rust
# // '#' needed for doctest
# use bed_reader::{Bed, ReadOptions, assert_eq_nan, sample_bed_file};
use ndarray::s;

let file_name = sample_bed_file("some_missing.bed")?;
let mut bed = Bed::new(file_name)?;
let val = ReadOptions::builder()
    .iid_index(s![..;2])
    .sid_index(20..30)
    .f64()
    .read(&mut bed)?;

assert!(val.dim() == (50, 10));
# use bed_reader::BedErrorPlus; // '#' needed for doctest
# Ok::<(), BedErrorPlus>(())
```

List the first 5 individual (sample) ids, the first 5 SNP (variant) ids,
and every unique chromosome. Then, read every genomic value in chromosome 5.

```rust
# use ndarray::s; // '#' needed for doctest
# use bed_reader::{Bed, ReadOptions, assert_eq_nan, sample_bed_file};
# let file_name = sample_bed_file("some_missing.bed")?;
use std::collections::HashSet;

let mut bed = Bed::new(file_name)?;
println!("{:?}", bed.iid()?.slice(s![..5])); // Outputs ndarray: ["iid_0", "iid_1", "iid_2", "iid_3", "iid_4"]
println!("{:?}", bed.sid()?.slice(s![..5])); // Outputs ndarray: ["sid_0", "sid_1", "sid_2", "sid_3", "sid_4"]
println!("{:?}", bed.chromosome()?.iter().collect::<HashSet<_>>());
// Outputs: {"12", "10", "4", "8", "19", "21", "9", "15", "6", "16", "13", "7", "17", "18", "1", "22", "11", "2", "20", "3", "5", "14"}
let val = ReadOptions::builder()
    .sid_index(bed.chromosome()?.map(|elem| elem == "5"))
    .f64()
    .read(&mut bed)?;

assert!(val.dim() == (100, 6));
# use bed_reader::BedErrorPlus; // '#' needed for doctest
# Ok::<(), BedErrorPlus>(())
```

Project Links
-----

* [**Installation**](https://crates.io/crates/bed-reader)
* [**Documentation**](https://docs.rs/bed-reader/)
* [**Questions via email**](mailto://fastlmm-dev@python.org)
* [**Source code**](https://github.com/fastlmm/bed-reader)
* [**Discussion**](https://github.com/fastlmm/bed-reader/discussions/)
* [**Bug Reports**](https://github.com/fastlmm/bed-reader/issues)
* [**Project Website**](https://fastlmm.github.io/)

