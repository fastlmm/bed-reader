# Useful

## Python

### install packages, compile Rust for Python

```bash
pip install -r requirements.txt
pip install -r requirements-dev.txt
maturin develop
```

### test

```bash
pytest bed_reader
pytest --doctest-modules bed_reader\_open_bed.py
pytest --doctest-modules README.md
pytest bed_reader/tests/test_open_bed_cloud.py::test_http_two -s
pytest bed_reader/tests/test_open_bed_cloud.py::test_http_cloud_urls_rst_1 -s
pytest --doctest-modules doc/index.rst

pytest --collect-only   # check test discovery
```

### generate doc

```bash
doc> make html & build\html\index.html
```

### create a local *.whl

```bash
maturin build --release
```

### Release

```bash
cd doc
make.bat html
build\html\index.html
xcopy /c /e /s /h build\html ..\docs
# push changes





## Rust

```bash
cargo doc --no-deps --open
cargo check

cargo test --doc supplemental_document_cloud_urls

cargo test --no-default-features
cargo test --all-features
cargo check --no-default-features
maturin develop
cargo publish --dry-run
```

## CMD

set PROMPT=$P$G
