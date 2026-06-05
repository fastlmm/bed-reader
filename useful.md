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

Inspect the `ls -las target/package/tmp-crate/*.crate` path printed by `cargo publish --dry-run` and confirm it is not unexpectedly large. Around 510KB is fine.

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

* 25 wheels: Python 3.10 through 3.14 for macOS x86_64, macOS aarch64, Linux x86_64, Linux aarch64, and Windows x64.
* 1 sdist: `bed_reader-VERSION.tar.gz`.

If the GitHub Actions UI does not show the artifact list, download the merged artifact with `gh`:

```bash
gh run download RUN_ID --name wheels-all --dir wheels
```

### 9. Inspect CI Artifacts

Download and extract the `wheels-all` artifact from GitHub Actions. It may contain the files directly, or inside the directory passed to `--dir`.

Windows example after downloading:

```cmd
cd /d "C:\Users\carlk\Downloads\wheels"
dir
```

Linux/macOS example after downloading:

```bash
find wheels -type f | sort
```

Check the artifacts before upload:

```bash
uvx twine check wheels/bed_reader-1.1.0*
```

Use the actual version number. Look for expected Python versions, platforms, and file names:

* `cp310`, `cp311`, `cp312`, `cp313`, `cp314`
* `macosx_10_12_x86_64`
* `macosx_11_0_arm64`
* `manylinux_2_17_x86_64.manylinux2014_x86_64`
* `manylinux_2_17_aarch64.manylinux2014_aarch64`
* `win_amd64`
* `bed_reader-VERSION.tar.gz`

### 10. Optional: TestPyPI

Upload to TestPyPI:

```bash
uvx twine upload --repository testpypi wheels/bed_reader-VERSION*
```

Install from TestPyPI in a clean environment:

```bash
uv venv --python 3.14 /tmp/bed-reader-testpypi
source /tmp/bed-reader-testpypi/bin/activate
uv pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ bed-reader==VERSION
cd /tmp
python -c "import bed_reader; from bed_reader import open_bed; print(bed_reader.__file__); print(open_bed)"
deactivate
```

On Windows, use:

```cmd
uv venv --python 3.14 %TEMP%\bed-reader-testpypi
%TEMP%\bed-reader-testpypi\Scripts\activate
uv pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ bed-reader==VERSION
cd /d %TEMP%
python -c "import bed_reader; from bed_reader import open_bed; print(bed_reader.__file__); print(open_bed)"
deactivate
```

### 11. Publish

After TestPyPI and CI look good, merge to `main` and wait for CI on `main`.

```bash
git status
git switch main
git pull
git merge --no-ff RELEASE_BRANCH
git push origin main
```

Use the actual release branch name. Wait for GitHub Actions to pass on `main` before publishing.

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
uvx twine upload wheels/bed_reader-VERSION*
```

If publishing from WSL but credentials live in Windows, copy the PyPI config once:

```bash
cp /mnt/c/Users/carlk/.pypirc ~/.pypirc
chmod 600 ~/.pypirc
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
uv pip install bed-reader==VERSION
cd /tmp
python -c "import bed_reader; from bed_reader import open_bed; print(bed_reader.__file__); print(open_bed)"
uv pip install "bed-reader[samples]==VERSION"
python -c "from bed_reader import sample_file; print(sample_file('small.bed'))"
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
