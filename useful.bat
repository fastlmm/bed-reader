# test async in python:
pytest bed_reader\tests\test_open_bed_cloud.py -k read1
pytest bed_reader\tests\test_open_bed.py -k test_zero_files


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

# show docs
cargo doc --open