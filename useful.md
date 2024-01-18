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

## Rust

```bash
cargo doc --no-deps --open
cargo check

cargo test --no-default-features
cargo check --no-default-features
maturin develop
```
