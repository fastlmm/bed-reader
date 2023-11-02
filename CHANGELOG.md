# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2023-11-02

### Changed

- (Python) Added an installation option `pip install bed-reader[base]` that depends only on `numpy` and does not need `pandas` or `pooch`.
- (Rust) Rust methods that returned `Result` used to return
  `Result<_, BedErrorPlus>`. To save memory, they now
  return `Result<_,Box<BedErrorPlus>>`.

## [0.2.27] - 2023-10-29

- (Python) Add support for Python 3.12.
