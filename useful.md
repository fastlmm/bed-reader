# Useful

## Release Checklist

Use this as the one deployment checklist for Python/PyPI, Rust/crates.io, docs, and GitHub artifacts.

### 1. Start Release Branch

```bash
git status
git switch main
git pull
git switch -c release-NEW_VERSION
```

Use the actual version number, for example `release-1.1.0`.

### 2. Bump Version

Update the release version in:

* `Cargo.toml`
* `pyproject.toml`
* `CHANGELOG.md`, if needed
* docs or README files, if they mention the release version

Use a normal final version, for example `1.1.0`, unless intentionally publishing a prerelease such as `1.1.0a1`.

`Cargo.toml` and `pyproject.toml` should have the same release version.

### 3. Preflight

```bash
git status # or use VSCode to see that everything is checked-in
rg -n "cmk" --glob '!doc/build/**' --glob '!docs/**' .
```

Update:

* `CHANGELOG.md`

Review:

* Consider which, if any Rust and Python versions to update to.

* Consider which, if any, Rust and Python dependencies should be updated for this release.

### 4. Rust Checks

Strict clippy is pinned so new monthly clippy lints do not unexpectedly fail CI.
After any Rust dependency or `Cargo.lock` update, rerun:

```bash
rustup run 1.83.0 cargo check --all-targets --all-features
rustup run 1.83.0 cargo clippy --all-targets --all-features -- -D warnings
```

Local stable Cargo may accept dependencies that CI's pinned Cargo cannot parse.

```bash
rustup toolchain install 1.83.0 --component clippy --component rustfmt --component rust-src
rustup run 1.83.0 cargo clippy --all-targets --all-features -- -D warnings
cargo fmt --check
cargo check
cargo check --no-default-features
cargo check --features extension-module
cargo test --no-default-features
cargo test --all-features
cargo test --doc supplemental_document_cloud_urls
```

Security and semver checks:

```bash
cargo install cargo-audit
cargo install cargo-semver-checks
cargo audit
# may not need the env set next time
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo semver-checks check-release
cargo publish --dry-run
```

Inspect the `.crate` path printed by `cargo publish --dry-run` and confirm it is not unexpectedly large. Around 510KB is fine.

Do not run `cargo update` as a routine release step. Use it only when a dependency update is deliberately part of the release.

### 5. Python Checks

Minimal dependency check, especially for Python 3.14:

```bash
uv python install 3.14
uv sync --python 3.14 --extra min_dev
source .venv/bin/activate
python --version
python -m pytest bed_reader/tests/test_opt_dep.py
```

Full Python test pass:

```bash
uv sync --python 3.14 --all-extras
python -m pytest .
python -m pytest --collect-only
```

Lint:

```bash
uvx ruff@0.15.16 check bed_reader
uvx ruff@latest check bed_reader
```

Aside, useful targeted checks (already covered above)

```bash
python -m pytest --doctest-modules bed_reader/_open_bed.py
python -m pytest --doctest-modules README.md
python -m pytest --doctest-modules doc/source/index.rst
python -m pytest bed_reader/tests/test_open_bed_cloud.py::test_http_two -s
python -m pytest bed_reader/tests/test_open_bed_cloud.py::test_http_cloud_urls_rst_1 -s
```

### 6. Docs

Rust docs:

```bash
cargo doc --no-deps --open
```

Python docs on Windows:

```cmd
cd doc
make.bat html
build\html\index.html
cd ..
if not exist docs mkdir docs
xcopy /c /e /s /h doc\build\html docs
```

Python docs on Linux/macOS:

```bash
cd doc
make html
python -m webbrowser build/html/index.html
cd ..
mkdir -p docs
cp -R doc/build/html/. docs/
```

Check the rendered README and docs pages before committing generated docs.
The PyPI version and Python-version badges are live badges. Before publishing, they may still show the previous release; that is expected.

### 7. Local Wheel Smoke Test

```bash
uv sync --all-extras
maturin build --release
python -m twine check target/wheels/*
```

Optional local editable build:

```bash
maturin develop
python -m pytest bed_reader/tests/test_opt_dep.py
```

### 8. Commit, Push, And Wait For CI

```bash
git status
git diff --check
git add Cargo.toml Cargo.lock pyproject.toml CHANGELOG.md README.md README-rust.md doc src bed_reader tests .github useful.md
test ! -d docs || git add docs
git commit -m "Support Python 3.14"
git push
```

Wait for GitHub Actions to finish. The `wheels-all` artifact should contain:

* macOS x86_64 wheels
* macOS aarch64 wheels
* Linux x86_64 wheels
* Linux aarch64 wheels
* Windows x64 wheels
* sdist

### 9. Inspect CI Artifacts

Download and extract the `wheels-all` artifact from GitHub Actions.

Windows example:

```cmd
cd /d "C:\Users\carlk\Downloads\wheels"
dir
```

Check the artifacts before upload:

```bash
python -m twine check bed_reader-*
```

Look for expected Python versions, platforms, and file names. For this release, expect Python 3.10 through 3.14 wheels.

### 10. TestPyPI

Upload to TestPyPI:

```bash
python -m twine upload --repository testpypi bed_reader-*
```

Install from TestPyPI in a clean environment:

```bash
uv venv --python 3.14 /tmp/bed-reader-testpypi
source /tmp/bed-reader-testpypi/bin/activate
python -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ bed-reader
python -c "from bed_reader import open_bed; print(open_bed)"
python -m pip uninstall -y bed-reader
deactivate
```

On Windows, use:

```cmd
uv venv --python 3.14 %TEMP%\bed-reader-testpypi
%TEMP%\bed-reader-testpypi\Scripts\activate
python -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ bed-reader
python -c "from bed_reader import open_bed; print(open_bed)"
python -m pip uninstall -y bed-reader
deactivate
```

### 11. Publish

After TestPyPI and CI look good, merge to `main` and wait for CI on `main`.

Publish Rust crate:

```bash
cargo publish
```

Rust dependency smoke test:

```bash
# In a small outside Rust project, such as InstallTestRust:
cargo add bed-reader@NEW_VERSION
cargo run
```

Publish Python package:

```bash
python -m twine upload bed_reader-*
```

Tag the release:

```bash
git tag -a v1.1.0 -m "v1.1.0"
git push origin v1.1.0
```

Use the actual version number.

### 12. Post-Release Verification

```bash
uv venv --python 3.14 /tmp/bed-reader-pypi
source /tmp/bed-reader-pypi/bin/activate
python -m pip install bed-reader
python -c "import bed_reader; from bed_reader import open_bed; print(bed_reader.__version__ if hasattr(bed_reader, '__version__') else open_bed)"
deactivate
```

Check:

* PyPI page: <https://pypi.org/project/bed-reader/>
* docs.rs page: <https://docs.rs/bed-reader/latest/bed_reader/>
* GitHub release/tag page
* project docs: <https://fastlmm.github.io/bed-reader>

## Development Notes

### PyO3 Debugging

To get some Rust PyO3 debugging to work, set this before starting VS Code:

```cmd
set PYO3_PYTHON=C:/Users/carlk/OneDrive/programs/bed-reader/.venv\Scripts\python.exe
```

### Benchmarking

```cmd
uv sync --extra samples --extra sparse --extra dev
uv pip install matplotlib
cd C:\Users\carlk\OneDrive\programs\bed-reader\bed_reader\tests\benchmark
python benchmark.py
```

Look in `M:\deldir\bench` for results.

### Windows Prompt

```cmd
set PROMPT=$P$G
```
