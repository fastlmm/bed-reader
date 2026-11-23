# Useful

## Release Checklist

Use this as the deployment checklist for the `1.3.0b3` Python beta, the
corresponding `1.3.0-beta.3` Rust beta, documentation, and release artifacts.

### 1. Prepare And Review The Release Branch

The release changes are prepared and reviewed on `implement/cpu-only`, and the
beta is tagged there. `main` keeps 1.2.1 (and its published documentation)
until the final 1.3.0, when `implement/cpu-only` is squash-merged into `main`.
Confirm that the versions remain:

* `version = "1.3.0-beta.3"` in `Cargo.toml`
* `version = "1.3.0b3"` in `pyproject.toml`

Also confirm that `CHANGELOG.md`, README files, and documentation describe the
beta accurately. Do not use the Python spelling in `Cargo.toml` or the Rust
SemVer spelling in `pyproject.toml`.

### 2. Preflight

```bash
git status # or use VSCode to see that everything is checked-in
rg -n "cmk" --glob '!doc/build/**' --glob '!docs/**' .
```

Update:

* `CHANGELOG.md`: in the release-preparation commit, before tagging, rename
  `## [Unreleased]` to `## [1.3.0b3] - YYYY-MM-DD` with the actual release
  date. Do not commit a placeholder date.

Review:

* Consider which, if any Rust and Python versions to update to.

* Consider which, if any, Rust and Python dependencies should be updated for this release.

### 3. Rust Checks

Strict clippy is pinned so new monthly clippy lints do not unexpectedly fail CI.
After any Rust dependency or `Cargo.lock` update, rerun:

```bash
rustup run 1.85.0 cargo check --locked --all-targets --all-features
rustup run 1.85.0 cargo clippy --locked --all-targets --all-features -- -D warnings
```

Local stable Cargo may accept dependencies that CI's pinned Cargo cannot parse.

```bash
rustup toolchain install 1.85.0 --component clippy --component rustfmt --component rust-src
rustup run 1.85.0 cargo clippy --locked --all-targets --all-features -- -D warnings
cargo fmt --check
cargo check --locked
cargo check --locked --no-default-features
cargo check --locked --features extension-module
cargo test --locked --no-default-features
cargo test --locked --all-features
cargo test --doc supplemental_document_cloud_urls
```

Security and semver checks:

```bash
cargo install cargo-audit
cargo install cargo-semver-checks
cargo audit
# may not need the env set next time
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo semver-checks check-release
rm -rf target/package
cargo publish --dry-run
```

Inspect the `ls -las target/package/tmp-crate/*.crate` path printed by `cargo publish --dry-run` and confirm it is not unexpectedly large. Around 510KB is fine.
If verification fails with a missing file that is present in the `.crate`, remove stale Cargo package output and retry:

```bash
rm -rf target/package
cargo publish --dry-run
```

Do not run `cargo update` as a routine release step. Use it only when a dependency update is deliberately part of the release.

### 4. Python Checks

Minimal dependency check, including the newest supported Python:

```bash
uv python install 3.15
uv sync --frozen --python 3.15 --extra min_dev
source .venv/bin/activate
python --version
python -m pytest bed_reader/tests/test_opt_dep.py
```

Full Python test pass:

```bash
uv sync --frozen --python 3.15 --all-extras
python -m pytest .
python -m pytest --collect-only
```

Lint:

```bash
uvx ruff@0.17.5 check bed_reader
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

### 5. Docs

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
```

Python docs on Linux/macOS:

```bash
cd doc
make html
python -m webbrowser build/html/index.html
cd ..
```

Check the rendered README and docs pages. Generated HTML is not committed:
docs.yml builds the site on each push to main and deploys it to GitHub Pages.
The PyPI version and Python-version badges are live badges. Before publishing, they may still show the previous release; that is expected.

### 6. Local Distribution Smoke Test

```bash
uv sync --all-extras
maturin build --release --locked --out dist
maturin sdist --out dist
uvx twine check --strict dist/*
```

Optional local editable build:

```bash
maturin develop
python -m pytest bed_reader/tests/test_opt_dep.py
```

### 7. Commit, Push, And Wait For CI

```bash
git status
git diff --check
git add Cargo.toml Cargo.lock pyproject.toml CHANGELOG.md README.md README-rust.md doc src bed_reader tests .github useful.md
test ! -d docs || git add docs
git commit -m "Prepare bed-reader 1.3.0 beta 3"
git push origin implement/cpu-only
gh workflow run ci.yml --ref implement/cpu-only
```

Wait for the manually dispatched CI run to finish. The `wheels-all` artifact
should contain:

* 30 wheels: Python 3.11 through 3.15 for macOS x86_64, macOS aarch64,
  Linux x86_64, Linux aarch64, Windows x64, and Windows ARM64.
* 1 source distribution: `bed_reader-1.3.0b3.tar.gz`.

If the GitHub Actions UI does not show the artifact list, download the merged artifact with `gh`:

```bash
gh run download RUN_ID --name wheels-all --dir wheels
```

### 8. Inspect CI Artifacts

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
uvx twine check --strict wheels/bed_reader-1.3.0b3*
```

Use the actual version number. Look for expected Python versions, platforms, and file names:

* `cp311`, `cp312`, `cp313`, `cp314`, `cp315`
* `macosx_10_12_x86_64`
* `macosx_11_0_arm64`
* `manylinux_2_17_x86_64.manylinux2014_x86_64`
* `manylinux_2_17_aarch64.manylinux2014_aarch64`
* `win_amd64`
* `win_arm64`
* `bed_reader-1.3.0b3.tar.gz`

The workflows do not upload to TestPyPI. Do not run a separate `twine upload`:
the release workflow publishes the qualified files directly to PyPI.

### 9. Tag The Beta

After review and CI succeed, tag `implement/cpu-only`; do not merge the beta into
`main`, because `docs.yml` deploys the documentation site from `main`.
`release.yml` runs on any `v*` tag, whatever the branch. Confirm that the exact
commit has the beta versions and dated release notes, then create the
Python-version tag expected by `release.yml`:

```bash
git status
git switch implement/cpu-only
git pull
git show --stat
git tag -a v1.3.0b3 -m "v1.3.0b3"
git push origin v1.3.0b3
```

Pushing `v1.3.0b3` starts `.github/workflows/release.yml`. Its first job verifies
that the tag is `v1.3.0b3`, the Python version is `1.3.0b3`, and the Rust
version normalizes from `1.3.0-beta.3`. It then calls the complete CI workflow,
including wheel and source-distribution builds and installed-artifact tests.

The release workflow does not create a GitHub Release.

### 10. Approve The PyPI Deployment

After every qualification job passes, download and inspect the release run's
`wheels-all` artifact as described above. The publish job uses the protected
`pypi` GitHub environment and waits for its required reviewer. Approve that
deployment only after the artifacts are accepted. The job then publishes the
same qualified files directly to PyPI using trusted publishing.

There is no manual GPU qualification: this release and its wheels are CPU-only.

### 11. Publish The Rust Crate Manually

No GitHub Actions workflow publishes to crates.io. From the exact tagged commit,
inspect the dry run again and then publish `1.3.0-beta.3` manually:

```bash
rm -rf target/package
cargo publish --dry-run --locked
cargo publish --locked
```

Rust dependency smoke test:

```bash
# In a small outside Rust project, such as InstallTestRust:
cargo add bed-reader@1.3.0-beta.3
cargo run
```

### 12. Post-Release Verification

```bash
uv venv --python 3.15 /tmp/bed-reader-pypi
source /tmp/bed-reader-pypi/bin/activate
uv pip install bed-reader==1.3.0b3
cd /tmp
python -c "import bed_reader; from bed_reader import open_bed; print(bed_reader.__file__); print(open_bed)"
uv pip install "bed-reader[samples]==1.3.0b3"
python -c "from bed_reader import sample_file; print(sample_file('small.bed'))"
deactivate
```

Check:

* PyPI page: <https://pypi.org/project/bed-reader/>
* docs.rs page: <https://docs.rs/bed-reader/latest/bed_reader/>
* GitHub tag page
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
