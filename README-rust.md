cmk badges go here

[![crates.io](https://img.shields.io/crates/v/rustybed.svg)](https://crates.io/crates/rustbed)

Read and write the PLINK BED format, simply and efficiently. 

Features:

* Fast and multi-threaded
* Supports many indexing methods. Slice data by individuals (samples) and/or SNPs (variants).
* Used by [PySnpTools](https://github.com/fastlmm/PySnpTools), [FaST-LMM](https://github.com/fastlmm/FaST-LMM), and [PyStatGen](https://github.com/pystatgen). cmk update
* Supports [PLINK 1.9](https://www.cog-genomics.org/plink2/formats).

Setup
=======

Add this to your Cargo.toml:

```
[dependencies]
bed_reader = "0.1.0"
```

Usage
========

Read genomic data from a .bed file.

```rust
use bed_reader;

fn main() {
    let file_name = sample_file("small.bed");
    mut bed = bed_reader::open_bed(file_name)?;
    let val = bed.read()?;
    println!("{:?}", val);
}
```

Read every second individual and SNPs (variants) from 20 to 30.

```rust
use bed_reader;

fn main() {
    let file_name2 = sample_file("some_missing.bed");
    mut bed2 = bed_reader::open_bed(file_name2)?;
    let val2 = bed2.read(index=(::2,20:30))?;
    println!("{:?}", val2);
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

Project Links
==============

cmk update
- [**Documentation**](http://fastlmm.github.io/bed-reader)
- **Questions to**: [fastlmm-dev@python.org](mailto:fastlmm-dev@python.org)
- [**Source code**](https://github.com/fastlmm/bed-reader)
- [**PyPI**](https://pypi.org/project/bed-reader)
- [**Bug reports**](https://github.com/fastlmm/bed-reader/issues)
- [**Mailing list**](https://mail.python.org/mailman3/lists/fastlmm-user.python.org)
- [**Project Website**](https://fastlmm.github.io/)