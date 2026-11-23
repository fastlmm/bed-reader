# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.0b1] - 2026-11-22

* Withdraw the experimental GPU features introduced in 1.3.0a1: CuPy reads,
  GPU allocation and standardization helpers, and Rust CubeCL backends.
  bed-reader reads into CPU arrays. PySnpTools and FaST-LMM retain GPU
  computation by transferring CPU reads to CuPy.

* Behavior change from 1.2.1: an explicit `dtype=None` in `open_bed.read` and
  `open_bed.read_sparse` now returns float32, not NumPy's implicit float64.
  An omitted or explicit `dtype=None` uses the `ARRAY_DTYPE` environment
  variable (`float32` or `float64`) when set, otherwise float32. Explicit
  dtypes override the environment.

* Reject unsupported explicit dtypes (anything but `int8`, `float32`, or
  `float64`) with `ValueError` in `open_bed.read` and `open_bed.read_sparse` on
  every path. Previously `force_python_only=True` and files with no
  individuals or SNPs returned arrays of the unsupported dtype.

* Rust: `numpy` is now enabled only by the `extension-module` feature and
  `ndarray-npy` is a test-only dependency, so non-Python builds compile fewer
  dependencies.

* Fix index errors on mode 0 (individual-major) .bed files: an out-of-range
  individual index was reported as a SNP index, and the reverse. Local and
  cloud reads now resolve indexes with one shared function.

## [1.3.0a1] - 2026-11-17

* Add experimental NVIDIA GPU reads for Python through CuPy and optional
  CubeCL CUDA/HIP backends for Rust. CPU execution remains the default.
  Python GPU wheels target Windows and Linux x86-64; WSL2 uses the Linux wheel.
  The `cupy` extra installs CuPy's CUDA 13 distribution and follows CuPy's
  NVIDIA driver/toolkit prerequisites.

* Support GPU Unit/Beta standardization of fractional float32/float64 inputs,
  with centered variance calculation, missing values, constant columns,
  trained statistics, and both memory orders. Int8 standardization remains
  unsupported. Warn about float64 GPU arithmetic costs without changing dtype.

* Add experimental Rust backend selection while retaining CPU-only builds.
  Optional GPU dependencies require Rust 1.95 or later.

* Add `ARRAY_DTYPE` precision defaults for Python dense/sparse reads and GPU
  allocation. Explicit dtype wins; without the setting, existing defaults are
  unchanged. The optional dtype parameter now defaults to `None`.
  Explicit `dtype=None` now selects the same default as omission (float32 for
  reads without the environment), rather than NumPy's implicit float64 dtype.

* Add `is_cubecl_backed` to identify CubeCL-backed CuPy arrays and their views
  without copying values or synchronizing the GPU.

* Add experimental `empty_gpu` allocation of CubeCL-backed CuPy matrices, with
  view lifetime tracking and safe native owner transfer between Python threads.
  Expose `gpu_memory_cleanup` to release unused native allocations without
  invalidating live arrays.

* Retry sample-file downloads up to three times on transient network errors.

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
