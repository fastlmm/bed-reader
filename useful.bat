pip# test Rust
cargo test

# install packages, compile Rust for Python, and test Python
pip install -r requirements.txt
pip install -r requirements-dev.txt
maturin develop
pytest bed_reader
maturin develop --target-dir c:\deldir\

pytest --doctest-modules bed_reader\_open_bed.py

# generate doc
\doc>make html & build\html\index.html

# create a local *.whl
maturin build --release