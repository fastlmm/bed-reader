# Useful

## Note

To get some Rust PyO3 debugging to work I had to set this before starting vscode:

```cmd
set PYO3_PYTHON=C:/Users/carlk/OneDrive/programs/bed-reader/.venv\Scripts\python.exe
```

## Python

### test

```bash
uvx ruff check .
uv sync --extra min_dev
pytest bed_reader/tests/test_opt_dep.py
uv sync --extra samples --extra sparse --extra dev
pytest .

pytest --doctest-modules bed_reader\_open_bed.py
pytest --doctest-modules README.md
pytest bed_reader/tests/test_open_bed_cloud.py::test_http_two -s
pytest bed_reader/tests/test_open_bed_cloud.py::test_http_cloud_urls_rst_1 -s
pytest --doctest-modules doc/source/index.rst

pytest --collect-only   # check test discovery
```

### Benchmarking

```cmd
uv sync --extra samples --extra sparse --extra dev
uv pip install matplotlib
cd C:\Users\carlk\OneDrive\programs\bed-reader\bed_reader\tests\benchmark
python benchmark.py
```

Look in `M:\deldir\bench` for results.

### install packages, compile Rust for Python

```bash
uv sync --extra samples --extra sparse --extra dev
```

### generate doc

```bash
doc> make html & build\html\index.html
```

### create a local *.whl

```bash
maturin build --release
```

### Release to PyPi

```bash
cd doc
make.bat html
build\html\index.html
xcopy /c /e /s /h build\html ..\docs
# push changes
# download artifacts to .e.g. C:\Users\carlk\Downloads\wheels (25)
twine upload --repository testpypi bed_reader-1*
pip install --index-url https://test.pypi.org/simple/ "bed-reader"
pip uninstall bed-reader
# test with install tests, look at README and web pages
twine upload bed_reader-1*

## Rust

### New Clippy Rules
```bash
rustup update
cargo --version
# Update ci.yml to this new version
cargo clippy --all-targets --all-features -- -D warnings
```

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
