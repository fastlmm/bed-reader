<!-- This README.md must be keep in sync with doc/source/index.rst manually -->

<!-- markdownlint-disable MD046 -->
<!-- MD046 is about 'code block style' -->
<!-- markdownlint-disable MD041 -->
<!-- MD041 is about 'first line in a file should be a top-level heading'. -->
[![PyPI version](https://badge.fury.io/py/bed-reader.svg)](https://badge.fury.io/py/bed-reader)
[![Build Status](https://github.com/fastlmm/bed-reader/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/fastlmm/bed-reader/actions/workflows/ci.yml)
![PyPI](https://img.shields.io/pypi/v/bed-reader)
<!-- markdownlint-enable MD041 -->

Read and write the PLINK BED format, simply and efficiently.

*This is the Python README. For Rust, see [README-rust.md](https://crates.io/crates/bed-reader).*

Highlights
====================

* Fast multi-threaded Rust engine.
* Supports all Python indexing methods. Slice data by individuals (samples) and/or SNPs (variants).
* Used by [PySnpTools](https://github.com/fastlmm/PySnpTools), [FaST-LMM](https://github.com/fastlmm/FaST-LMM), and [PyStatGen](https://github.com/pystatgen).
* Supports [PLINK 1.9](https://www.cog-genomics.org/plink2/formats).
* Read data locally or from the cloud, efficiently and directly.

Install
====================

**Full version**: With all optional dependencies:

    pip install bed-reader[samples,sparse]

**Minimal version**: Depends only on `numpy`:

    pip install bed-reader

Usage
========

Read genomic data from a .bed file.

```python
>>> import numpy as np
>>> from bed_reader import open_bed, sample_file
>>>
>>> file_name = sample_file("small.bed")
>>> bed = open_bed(file_name)
>>> val = bed.read()
>>> print(val)
[[ 1.  0. nan  0.]
 [ 2.  0. nan  2.]
 [ 0.  1.  2.  0.]]
>>> del bed

```

Read every second individual and SNPs (variants) from 20 to 30.

```python
>>> file_name2 = sample_file("some_missing.bed")
>>> bed2 = open_bed(file_name2)
>>> val2 = bed2.read(index=np.s_[::2,20:30])
>>> print(val2.shape)
(50, 10)
>>> del bed2

```

List the first 5 individual (sample) ids, the
first 5 SNP (variant) ids, and every unique
chromosome. Then, read every genomic value in chromosome 5.

```python
>>> with open_bed(file_name2) as bed3:
...     print(bed3.iid[:5])
...     print(bed3.sid[:5])
...     print(np.unique(bed3.chromosome))
...     val3 = bed3.read(index=np.s_[:,bed3.chromosome=='5'])
...     print(val3.shape)
['iid_0' 'iid_1' 'iid_2' 'iid_3' 'iid_4']
['sid_0' 'sid_1' 'sid_2' 'sid_3' 'sid_4']
['1' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '2' '20' '21' '22'
 '3' '4' '5' '6' '7' '8' '9']
(100, 6)

```

From the cloud: open a file and read data for one SNP (variant)
at index position 2.

```python
>>> with open_bed("https://raw.githubusercontent.com/fastlmm/bed-sample-files/main/small.bed") as bed:
...     val = bed.read(index=np.s_[:,2], dtype="float64")
...     print(val)
[[nan]
 [nan]
 [ 2.]]

```

Project Links
==============

* [**Documentation**](http://fastlmm.github.io/bed-reader)
* **Questions to**: [fastlmm-dev@python.org](mailto:fastlmm-dev@python.org)
* [**Source code**](https://github.com/fastlmm/bed-reader)
* [**PyPI**](https://pypi.org/project/bed-reader)
* [**Bug reports**](https://github.com/fastlmm/bed-reader/issues)
* [**Mailing list**](https://mail.python.org/mailman3/lists/fastlmm-user.python.org)
* [**Project Website**](https://fastlmm.github.io/)
* [**Change Log**](https://github.com/fastlmm/bed-reader/blob/master/CHANGELOG.md)
