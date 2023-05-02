# test Rust
cargo test

# install packages, compile Rust for Python, and test Python
pip install -r requirements.txt
pip install -r requirements-dev.txt
maturin develop
pytest bed_reader
maturin develop --target-dir c:\deldir\