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
* The [Python-facing API](https://pypi.org/project/bed-reader/) for this library is used by [PySnpTools](https://github.com/fastlmm/PySnpTools), [FaST-LMM](https://github.com/fastlmm/FaST-LMM), and [PyStatGen](https://github.com/pystatgen).
* Supports [PLINK 1.9](https://www.cog-genomics.org/plink2/formats).

Examples
--------

Read genomic data from a .bed file.

```rust
use bed_reader;

fn main() {
    let file_name = sample_file("small.bed");
    mut bed = bed_reader::open_bed(file_name)?;
    let val = bed.read()?;
    println!("{val:?}");
}
```

Read every second individual and SNPs (variants) from 20 to 30.

```rust
use bed_reader;

fn main() {
    let file_name2 = sample_file("some_missing.bed");
    mut bed2 = bed_reader::open_bed(file_name2)?;
    let val2 = bed2.read(index=(::2,20:30))?;
    println!("{val2:?}");
}
```

List the first 5 individual (sample) ids, the
first 5 SNP (variant) ids, and every unique
chromosome. Then, read every value in chromosome 5.

use bed_reader;

cmk See slicing macro s! https://docs.rs/ndarray/latest/ndarray/macro.s.html

```rust
fn main() {
    let file_name2 = sample_file("some_missing.bed");
    mut bed3 = bed_reader::open_bed(file_name2)?;
    println!("{:?}", bed3.iid[:5]);
    println!("{:?}", bed3.sid[:5]);
    println!("{:?}", bed3.chromosome.unique());
    let val3 = bed.read(index=np.s_[:,bed3.chromosome=='5'])?;
    println!("{:?}", val3);
}
```
cmk how do you show output?

Links
-----

cmk update
- **Questions to**: [fastlmm-dev@python.org](mailto:fastlmm-dev@python.org)
- [**Bug reports**](https://github.com/fastlmm/bed-reader/issues)
- [**Mailing list**](https://mail.python.org/mailman3/lists/fastlmm-user.python.org)
- [**Project Website**](https://fastlmm.github.io/)