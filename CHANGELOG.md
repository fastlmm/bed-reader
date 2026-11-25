# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.0] - in beta (changes since 1.2.1)

### Highlights

* **Faster local reads, especially with many threads.** Each thread reads its own run
  of SNPs with positional file reads, in place of one serialized reader. Results are
  unchanged. `num_threads` is now a maximum: a small read may use fewer threads.
* **Precision default.** The `ARRAY_DTYPE` environment variable (`float32` or
  `float64`) sets the default dtype for `open_bed.read` and `open_bed.read_sparse`; an
  explicit dtype wins. Without it, reads still default to float32.

### Behavior changes

* An explicit `dtype=None` in `open_bed.read` and `open_bed.read_sparse` now means the
  default (`ARRAY_DTYPE`, then float32), not NumPy's implicit float64.
* Unsupported explicit dtypes (anything but `int8`, `float32`, or `float64`) raise
  `ValueError` on every path. Before, `force_python_only=True` and files with no
  individuals or SNPs returned arrays of the unsupported dtype.
* When several indexes are out of range, the error now always names the first one.

### Fixes

* Constant fractional SNPs (for example, dosages that are all 0.1) now standardize to
  exactly zero with `standardize_f32` and `standardize_f64` (used by PySnpTools' Unit
  and Beta), and training stores an infinite standard deviation, so applying the
  trained statistics to new values also gives zero. Before, a rounded mean left a tiny
  variance: 0.1 repeated 1,000 times standardized to values up to 1.3e-7 in float64
  (1.7 in float32: up to 2.8e-3), and applying the trained statistics gave results up
  to 3.8e8.
* Float32 standardization accumulates per-SNP sums in float64. With float32 sums, a SNP
  whose values were mostly 2 (such as a rare variant counted by its common allele) had
  its standard deviation off by up to 3.4e-4 relative; it is now within float32
  rounding. Float32 results change slightly; float64 results are unchanged bit for bit.
* On mode 0 (individual-major) .bed files, an out-of-range individual index was
  reported as a SNP index, and the reverse. Local and cloud reads now share one index
  check.

### Other

* About 4 times faster float32 and float64 standardization of C-order (row-major)
  arrays (14,925 x 10,000 at 8 threads).
* Sample-file downloads retry up to three times on transient network errors.
* Rust: the crate uses edition 2024; the minimum supported Rust version remains 1.85.
  `numpy` is enabled only by the `extension-module` feature, and `ndarray-npy` is a
  test-only dependency, so non-Python builds compile fewer dependencies.
* The experimental GPU features of the 1.3.0a1 alpha were withdrawn; bed-reader reads
  into CPU arrays, and PySnpTools and FaST-LMM move the data to the GPU.

## [1.2.1] - 2026-11-13

* Fix GitHub history.

## [1.2.0] - 2026-11-13

- Require Python 3.11 or later and prepare Python 3.15 wheels with PyO3 0.30.1,
  Rust NumPy 0.30, and maturin 1.16.
- Require Rust 1.85 or later. Upgrade cloud-file to 0.3.1 and use its u64 byte
  ranges without narrowing file offsets. Rust callers supplying a CloudFile
  must use cloud-file 0.3.
- Refresh compatible dependencies and remove the vulnerable quick-xml chain.
  Retain shared ndarray 0.16 and compiler-compatible psm/crc-fast versions in
  the lockfile.
- Track the uv lockfile, pin CI actions, and test installed wheels and source
  distributions on boundary Python versions, including Linux and Windows ARM64.
- Preserve caller arrays when replacing missing integer metadata, including
  read-only arrays, and avoid NumPy's invalid NaN-to-integer conversion warning.
- Prepare a tag-triggered PyPI release workflow that reuses the full CI suite,
  validates metadata, and publishes the qualified artifacts with trusted publishing.
  Compare complete direct and source-derived wheels without ignoring native code.

## [1.1.1] - 2026-06-05

- Fix Rust README CI badge.

## [1.1.0] - 2026-06-04

- Add support for Python 3.14.

## [1.06] - 2024-10-17

- Add support for Python 3.13. Fix minor reported bug.

## [1.05] - 2024-7-1

- Add support for Python's Numpy 2

## [1.0.4] - 2024-4-24

- Fix a minor bug in helper functions used by PySnpTools.

## [1.0.3] - 2024-4-17

- On the Python file, add `create-bed`, a method for writing large
  bed files SNP-by-SNP (or individual-by-individual).
  Also, add full support for the less common individual-major files.

## [1.0.2] - 2024-3-16

- Add support for cloud files to both Rust and Python.

## [1.0.0] - 2023-11-5

### Changed

- (Python) Core project now depends only on `numpy`.
  If you want to download samples files or create
  sparse matrices, say `pip install bed-reader[samples,sparse]`.

- (Python) Added new `read_sparse()` method for creating
  `sympy.sparse.csc_matrix` and `sympy.sparse.csr_matrix`. This
  method can save memory when a *.bed file contains mostly zeros.

- (Rust) Rust methods that returned `Result` now return
  `Result<_,Box<BedErrorPlus>>`. Before, they returned
  `Result<_, BedErrorPlus>`. This saves memory.

## [0.2.27] - 2023-10-29

- (Python) Add support for Python 3.12.
