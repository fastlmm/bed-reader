# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

* Return a GPU read's memory to CUDA as soon as the array and its views are
  dropped, so CuPy and other libraries can use it without a manual cleanup.
  Read results no longer stay in bed-reader's memory pool.

* Make `empty_gpu`, `is_cubecl_backed`, `gpu_memory_cleanup`, and
  `array_dtype` internal; they support PySnpTools rather than end users.
  `ARRAY_DTYPE` itself is unchanged.

* Report CUDA errors while copying data to the GPU as Python exceptions
  instead of panics.

* Rust: share `ReadOptions<TVal>` between CPU, cloud, and GPU reads. GPU reads
  end in `ReadOptionsBuilder::read_cubecl`; `Bed::read_cubecl` reads everything and
  `Bed::read_cubecl_with_options` takes built options. GPU items live in `gpu`
  (`CubeClMatrix`, `client`, `memory_cleanup`, `Error`), present only with `cubecl-cuda`;
  GPU calls return `Box<gpu::Error>`, wrapping ordinary read errors. `BedErrorPlus`
  and existing CPU option signatures match 1.2. GPU support remains Beta.
  CPU-only builds support Rust 1.85; GPU support needs Rust 1.95. Use an MSRV-aware
  lockfile with CPU-compatible ICU dependencies. Keep the `cubecl-cuda` feature
  name and remove the untested AMD `cubecl-hip` Cargo feature. Internal
  `_cubecl-cpu-test` runs kernel tests on CubeCL's CPU runtime, including when
  `--all-features` is used.

* Rust: Rust 1.85 users no longer need `cargo update -p psm --precise 0.1.28`.
  `ndarray-npy` is now a test-only dependency and `numpy` is enabled only by
  `extension-module`, so neither reaches non-Python builds.

* Python: GPU reads raise `ValueError` when CuPy's current device is not
  device 0, as GPU allocation already did, instead of returning an array on
  device 0. To use another GPU, set `CUDA_VISIBLE_DEVICES`.

* Fix GPU standardization of constant fractional SNPs. A rounded mean (three
  0.1s average to 0.10000000000000002) gave them a tiny variance, so every value
  standardized to -1 (or a Beta-weighted equivalent) instead of 0.

* Fix index errors on mode 0 (individual-major) .bed files: an out-of-range
  individual index was reported as a SNP index, and the reverse. CPU, cloud,
  and GPU reads now resolve indexes with one shared function.

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
