# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0] - Unreleased

- Require Python 3.11 or later and prepare Python 3.15 wheels with PyO3 0.30.1,
  Rust NumPy 0.30, and maturin 1.16.
- Require Rust 1.85 or later. Upgrade cloud-file to 0.3.1 and use its u64 byte
  ranges without narrowing file offsets. Rust callers supplying a CloudFile
  must use cloud-file 0.3.
- Refresh compatible dependencies and remove the vulnerable quick-xml chain.
  Retain shared ndarray 0.16 and compiler-compatible psm/crc-fast versions in
  the lockfile.
- Track the uv lockfile, pin CI actions, and test installed wheels and source
  distributions on boundary Python versions, including Linux ARM64 CI.

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
