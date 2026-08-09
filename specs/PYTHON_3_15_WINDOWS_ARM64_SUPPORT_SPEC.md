# Python 3.15, Windows ARM64, and Release Modernization Specification

## Status

Proposed.

This specification covers the work tracked by:

- [bed-reader issue #33](https://github.com/fastlmm/bed-reader/issues/33),
  Python 3.15 support and release modernization.
- [bed-reader issue #34](https://github.com/fastlmm/bed-reader/issues/34),
  Windows ARM64 wheels and tests.

Implementation must not begin by mechanically changing classifiers or expected
test output. Each supported interpreter, dependency boundary, artifact, and
native platform must be qualified as described below.

## Objective

Release `bed-reader` with:

- Standard, GIL-enabled CPython 3.12, 3.13, 3.14, and 3.15 support.
- Native wheels and installed-artifact tests on Windows ARM64 for every
  supported Python version.
- The existing Linux x86_64, Linux ARM64, Windows x86_64, Intel macOS, and
  Apple Silicon macOS support preserved.
- A current uv-based Python development, CI, artifact, documentation, and
  PyPI publication workflow consistent with FaST-LMM and PySnpTools.
- Explicit Rust minimum-version, feature, package, and crates.io release
  qualification.
- Independent Python and Rust release versions and tag streams.

The support claim covers installation, the Python and Rust test suites, all
published Python extras, native wheel use, source-distribution builds,
documentation, and publication controls. A classifier alone is not support.

## Decisions

The following decisions are settled for this work:

1. Python 3.10 and 3.11 support will be dropped. The supported Python range is
   3.12 through 3.15.
2. Python support will follow the Scientific Python three-year window: add the
   newest annual Python and drop versions outside the window.
3. While Python 3.15 is a prerelease, CI will use the newest available 3.15
   candidate. Its jobs are required, not allowed failures.
4. Prerelease dependency resolution is permitted only in Python 3.15 jobs and
   only when required. Python 3.12 through 3.14 continue to use stable
   dependencies. Remove the exception when stable dependencies are available.
5. Free-threaded CPython 3.15 (`3.15t`) is out of scope. It requires separate
   PyO3 concurrency and safety qualification.
6. Windows ARM64 support covers Python 3.12 through 3.15 and every published
   Python extra. Partial support will not be advertised.
7. Windows ARM64 builds and tests will run natively on GitHub's
   `windows-11-arm` hosted runner. Cross-compilation or emulation alone is not
   acceptance evidence.
8. The existing manylinux scope is preserved. Adding musllinux wheels is out
   of scope.
9. Pull requests will build and clean-install the complete release wheel
   matrix, including Windows ARM64.
10. Python and Rust versions remain independent. Python-only work does not
    force a Rust crate version bump or crates.io release.
11. Python release tags use `py-vX.Y.Z`; Rust release tags use
    `rust-vX.Y.Z`. Exact next version numbers are release-time decisions.
12. Python development and packaging orchestration use uv. Cargo and rustup
    remain the required tools for Rust compilation and crates.io publication.
    The PyPI publishing action remains the appropriate Trusted Publishing
    boundary.
13. Generated Sphinx HTML will no longer be tracked in the default branch.
    GitHub Pages will deploy the exact documentation artifact built in CI.

## Current State

- `pyproject.toml` advertises Python 3.10 through 3.14 but does not declare
  `requires-python`.
- `uv.lock` currently declares Python 3.14 or newer, which conflicts with the
  advertised range.
- The Python distribution version is 1.1.0; the Rust crate version is 1.1.1.
- Python packaging uses maturin and enables the PyO3 `extension-module`
  feature.
- CI builds wheels for Python 3.10 through 3.14 on Linux x86_64 and ARM64,
  Windows x86_64, Intel macOS, and Apple Silicon macOS.
- CI tests source checkouts but does not clean-install and exercise every
  wheel it produces.
- Windows ARM64 is absent from the build and test matrices.
- CI uses a mixture of uv, direct virtual-environment activation, direct
  pytest calls, unpinned action tags, and an old maturin-generated layout.
- There is no tag-only PyPI release workflow or crates.io release workflow.
- Generated Sphinx HTML is tracked under `docs/` while a separate workflow
  also builds and publishes `doc/build/html`.
- Python documentation currently obtains its displayed version from
  `Cargo.toml`, which is incorrect when Python and Rust versions are
  intentionally independent.

## Supported Matrix

### Interpreters

- CPython 3.12
- CPython 3.13
- CPython 3.14
- CPython 3.15, including the newest release candidate during qualification

Only standard GIL-enabled CPython is included.

### Wheel platforms

| Platform | Architecture | GitHub runner | Required wheel tag |
| --- | --- | --- | --- |
| Linux | x86_64 | `ubuntu-latest` | compatible manylinux x86_64 |
| Linux | ARM64 | `ubuntu-24.04-arm` | compatible manylinux aarch64 |
| Windows | x86_64 | `windows-2022` | `win_amd64` |
| Windows | ARM64 | `windows-11-arm` | `win_arm64` |
| macOS | x86_64 | `macos-15-intel` | compatible macOS x86_64 |
| macOS | ARM64 | `macos-15` | compatible macOS arm64 |

The release contains 24 CPython wheels: four Python versions on six native
platform/architecture combinations, plus one source distribution.

### Source-distribution coverage

The sdist must be clean-installed and tested with Python 3.15 on all six
platform/architecture combinations. Linux additionally tests sdist installs
at the Python 3.12 floor. An sdist test must build from the downloaded sdist,
not from the repository checkout or a sibling build directory.

## Python Metadata and Dependencies

1. Add `requires-python = ">=3.12"`.
2. Remove the Python 3.10 and 3.11 classifiers and add Python 3.15.
3. Preserve the `samples` and `sparse` extras as public installation paths.
4. Continue declaring every direct runtime import as a direct dependency.
5. Replace historical dependency floors with versions that genuinely support
   the applicable Python and architecture combination.
6. Qualify NumPy minimums separately where Python-version wheel availability
   requires markers. Do not retain a lower bound that cannot install on the
   Python version to which it applies.
7. Qualify SciPy, pooch, pandas, documentation dependencies, and test tools on
   Windows ARM64. All required acceptance environments must install binary
   dependency wheels; silently compiling third-party numerical dependencies
   from source is not evidence of normal user support.
8. Test declared lowest direct dependencies on Linux with Python 3.12 and
   3.15. Test the lowest versions with native wheels on Windows ARM64.
9. Test a monthly latest-stable dependency solve. While 3.15 requires
   prereleases, isolate that solve from stable-version jobs.
10. Replace the legacy `[tool.uv].dev-dependencies` declaration with the
    standardized dependency-group form used by FaST-LMM and PySnpTools.
    Preserve published extras where removing them would break documented
    downstream workflows.
11. Regenerate and commit `uv.lock` with a Python range consistent with
    `requires-python`. Ordinary CI and release jobs use `--frozen`.

No dependency bound may be raised solely because it is convenient. Each new
minimum must be supported by installation or test evidence, and upper bounds
must be used only for an identified unqualified incompatibility.

## uv Tooling Policy

Use the same initial pins as the current FaST-LMM and PySnpTools workflows:

- `astral-sh/setup-uv` at the same immutable commit.
- uv 0.12.2.
- Matching immutable commits for checkout and artifact actions.

Subsequent coordinated tool upgrades are separate maintenance work.

Use uv for:

- Installing and selecting Python interpreters.
- Creating project and isolated artifact environments.
- Locking and syncing Python dependencies.
- Running Python tests, Ruff, Sphinx, and maturin from locked environments.
- Building the Python sdist and native wheels through the maturin backend.
- Installing the exact wheel or sdist selected for artifact tests.
- Running metadata checks through pinned tools, for example `uvx twine`.

Do not manually activate virtual environments in workflows. Use `uv run` or
an explicit interpreter selected by uv. Do not use `latest` tool specifiers in
required jobs.

Cargo, rustup, platform linkers, GitHub Pages actions, and package-index
Trusted Publishing actions are deliberate exceptions because uv does not
replace those responsibilities.

## Rust Toolchain and Feature Qualification

1. Add an explicit `rust-version` to `Cargo.toml`.
2. First qualify the existing Rust 1.83 CI baseline. If the locked dependency
   graph no longer supports it, raise the minimum deliberately and record the
   dependency that requires the change.
3. CI must test both the declared minimum Rust version and current stable
   Rust. Strict lint results must be reproducible against a pinned toolchain;
   a moving stable toolchain may provide advisory coverage but may not silently
   replace the pinned required check.
4. Run formatting and clippy without suppressing failures in required jobs.
5. Test the minimal Rust crate with `--no-default-features`.
6. Test the default Rust feature set.
7. Qualify the `extension-module` feature through maturin builds and Python
   tests. Do not substitute a problematic direct `cargo test --all-features`
   invocation for actual extension loading.
8. Run `cargo package --locked` and `cargo publish --dry-run --locked` against
   the exact source intended for a Rust release.
9. Inspect the packaged crate so required sources, license, README, and Cargo
   metadata are present and repository-only fixtures or build output are not.

## Continuous Integration

Replace the generated monolithic workflow with explicit, bounded jobs.

### Required pull-request jobs

1. **Python lint**
   - Run pinned Ruff through uv once on Linux.
2. **Rust format and lint**
   - Check formatting.
   - Run strict clippy on the declared minimum/pinned Rust toolchain.
   - Run advisory or required stable-toolchain coverage as explicitly chosen
     during implementation; failures may not be hidden with shell fallbacks.
3. **Rust tests**
   - Run minimal and default feature tests.
   - Cover Linux, Windows, and macOS on their supported architectures where
     platform behavior differs.
4. **Python source tests**
   - Test Python 3.12 through 3.15 on Linux x86_64, Windows x86_64, Intel
     macOS, and Apple Silicon macOS.
   - Add native Windows ARM64 source tests for all four Python versions.
   - Use the locked full-extra environment.
5. **Minimum dependencies**
   - Test Python 3.12 and 3.15 lowest-direct environments on Linux.
   - Test the lowest wheel-installable full-extra environment on Windows
     ARM64.
6. **Documentation**
   - Build from clean Sphinx sources with warnings treated as errors.
   - Do not compare with or depend on checked-in generated HTML.
7. **Wheel builds and artifact tests**
   - Build all four CPython wheels in one job per native platform runner where
     practical.
   - Upload immutable per-platform artifacts with missing-file failures.
   - Clean-install and smoke-test every wheel with its matching interpreter.
8. **Source-distribution tests**
   - Build the sdist once from clean source.
   - Download, build, install, and test it on every target as specified above.
9. **Cloud integration**
   - Run network-dependent cloud-reader tests as a separately named required
     job.
   - Use bounded retries only for identified transient operations.
   - Preserve the original failure and diagnostics after retries are
     exhausted.

### Scheduled dependency freshness

Run a monthly workflow that refreshes the lock resolution without modifying
the repository, installs the result, and runs the complete applicable suite.
Stable Python jobs must not consume prerelease dependencies. Python 3.15 may
do so only while the candidate ecosystem requires it.

### Artifact smoke tests

Each wheel smoke test must:

- Verify `platform.machine()` and Python version match the matrix target.
- Verify the wheel filename/tag matches the intended architecture, including
  `win_arm64` on Windows ARM64.
- Import `bed_reader` from the isolated artifact environment, not the source
  checkout.
- Create, write, reopen, and read a small BED dataset through the native
  extension and verify its shape, values, ordering, and missing-value behavior.
- Exercise representative `samples` and `sparse` functionality.
- Confirm required package metadata and license files are present.

Windows ARM64 dependency installation must reject x86_64 wheels and avoid
emulation. Record the interpreter architecture in workflow logs.

## Windows ARM64 Implementation Requirements

1. Use the native `windows-11-arm` runner.
2. Install CPython 3.12 through 3.15 with uv and verify each interpreter reports
   ARM64.
3. Install or select the native MSVC build tools required by the runner image.
4. Install and verify the `aarch64-pc-windows-msvc` Rust target when it is not
   already the host target.
5. Build `win_arm64` wheels for all supported Python versions.
6. Run the full Python suite with all extras on native ARM64.
7. Clean-install and smoke-test every ARM64 wheel.
8. Build and test the sdist on native ARM64 with Python 3.15.
9. Run the Windows ARM64 lowest-dependency qualification described above.
10. Treat absence of a required third-party ARM64 wheel as a release blocker.
    Resolve it by selecting an honest supported minimum or coordinating the
    prerequisite release, not by silently omitting an extra.

## Documentation Modernization

1. Treat `doc/source` and its configuration as the documentation source of
   truth.
2. Change the Python documentation version lookup to read
   `pyproject.toml`, not `Cargo.toml`.
3. Build documentation with the locked uv environment and fail on warnings.
4. Configure GitHub Pages to deploy an uploaded Actions artifact.
5. Preserve existing documentation URLs and any custom-domain behavior.
6. Verify an artifact-based deployment before removing generated files under
   `docs/` from the default branch.
7. Remove stale-output comparisons after generated HTML is no longer tracked.
8. Keep the documentation build as a required pull-request and Python-release
   gate.
9. Pin the official Pages configuration, upload, and deployment actions to
   immutable revisions.

Documentation deployment occurs from the default branch. Release workflows
build documentation as a gate but do not deploy a second copy.

## Independent Versioning and Tags

### Python distribution

- Version source: `[project].version` in `pyproject.toml`.
- Tag format: `py-vX.Y.Z`, including normal PEP 440 prerelease versions such
  as `py-v1.2.0rc1`.
- The Python release workflow must reject a tag that does not exactly match
  the normalized project version.
- A Python release embeds the repository's Rust extension code but does not
  imply publication of the Rust crate to crates.io.

### Rust crate

- Version source: `[package].version` in `Cargo.toml`.
- Tag format: `rust-vX.Y.Z`.
- The Rust release workflow must reject a tag that does not exactly match the
  Cargo package version.
- A Rust release does not imply publication of a Python distribution.

Python documentation displays the Python version. Rust documentation and
crates.io metadata display the Cargo version. Cross-project text must not imply
that the independently versioned packages always have equal versions.

## Python Release Workflow

Create a tag-triggered workflow for `py-v*` tags with these gates:

1. Checkout the exact tagged commit using pinned actions.
2. Verify the tag matches `pyproject.toml` and confirm the Cargo metadata is
   internally valid without requiring version equality.
3. Verify `uv.lock` is current with `uv lock --check`.
4. Run required Rust, Python, documentation, and cloud integration checks on
   the tagged source.
5. Build the complete 24-wheel native matrix and one sdist using locked uv and
   maturin tooling.
6. Run every wheel and sdist artifact test required by this specification.
7. Run metadata and package-content checks.
8. Merge only the tested distributions into the release artifact. Never
   rebuild distributions in the publication job.
9. Publish through the protected GitHub `pypi` environment with OIDC Trusted
   Publishing and attestations enabled.
10. Grant `id-token: write` only to the publication job. All build and test
    jobs retain read-only contents permission.

The first Python 3.15 candidate release should be a PEP 440 prerelease. It may
be published only after the complete required matrix passes. Downstream
FaST-LMM and PySnpTools qualification must consume that published prerelease,
not an unpublished sibling checkout. A stable `bed-reader` release waits for
stable Python 3.15-compatible dependency artifacts.

## Rust Release Workflow

Create a separate tag-triggered workflow for `rust-v*` tags:

1. Checkout and verify the tag against `Cargo.toml`.
2. Use the declared minimum Rust version and current stable Rust as required
   qualification environments.
3. Run format, strict clippy, minimal-feature tests, and default-feature tests.
4. Build and test the Python extension on the newest supported standard Python
   with uv, because shared Rust changes can regress the Python distribution.
5. Run `cargo package --locked` and inspect the crate artifact.
6. Run `cargo publish --dry-run --locked` before the publication boundary.
7. Publish with `cargo publish --locked` only from a protected `crates-io`
   environment using crates.io Trusted Publishing or the narrowest supported
   repository credential mechanism.

The Rust publication workflow never modifies or publishes the Python package.
Agent-driven qualification may prepare this workflow and dry-run artifacts,
but a person must authorize an actual release.

## Manual Repository Configuration

The implementation is incomplete until a maintainer configures and verifies:

- GitHub Pages source set to GitHub Actions.
- A protected `pypi` environment with required approval and tag restrictions
  for `py-v*`.
- A matching PyPI Trusted Publisher for the Python release workflow.
- A protected `crates-io` environment with required approval and tag
  restrictions for `rust-v*`.
- A matching crates.io Trusted Publisher or narrowly scoped publication
  credential.
- Branch protection requiring all non-scheduled pull-request jobs, including
  Windows ARM64 artifacts and documentation.

Publication workflows must fail clearly when this external configuration is
missing. They must not fall back to long-lived broad credentials.

## Acceptance Criteria

### Metadata and resolution

- [ ] `pyproject.toml` requires Python 3.12 or newer and advertises 3.12 through
      3.15 only.
- [ ] `uv.lock` covers the declared Python range and passes `uv lock --check`.
- [ ] Direct dependency lower bounds are installable and tested on every
      applicable Python/architecture boundary.
- [ ] Stable Python jobs resolve only stable dependencies.
- [ ] Any temporary Python 3.15 prerelease dependency exception is explicit
      and isolated.

### Rust

- [ ] `Cargo.toml` declares a tested minimum Rust version.
- [ ] Formatting and strict clippy pass without ignored failures.
- [ ] Minimal and default Rust feature tests pass.
- [ ] The PyO3 extension builds and loads on every required target.
- [ ] `cargo package --locked` and `cargo publish --dry-run --locked` pass.

### Python and artifacts

- [ ] The complete Python suite passes for CPython 3.12 through 3.15.
- [ ] All public extras pass on Windows ARM64.
- [ ] Exactly the intended 24 native wheels and one sdist are produced.
- [ ] Every wheel is clean-installed and smoke-tested on its native target.
- [ ] The sdist builds and passes tests on all six native targets with Python
      3.15 and on Linux with Python 3.12.
- [ ] Windows ARM64 jobs verify native interpreters and `win_arm64` artifacts.
- [ ] Cloud integration tests pass as a distinct required gate.
- [ ] Numerical, dtype, ordering, missing-value, and file-format behavior is
      unchanged unless an independently justified behavior change is approved.

### Documentation

- [ ] Sphinx builds warning-free from clean sources in the locked environment.
- [ ] GitHub Pages deploys the uploaded build artifact successfully.
- [ ] Existing important documentation URLs continue to resolve.
- [ ] Generated Sphinx HTML is no longer tracked in the default branch.
- [ ] CI and release workflows do not compare against checked-in generated
      HTML.

### Publication

- [ ] `py-v*` and `rust-v*` tags validate against their independent version
      sources.
- [ ] The Python workflow publishes only artifacts that passed the complete
      native matrix.
- [ ] PyPI publication uses the protected environment, Trusted Publishing, and
      attestations.
- [ ] Rust publication is independently gated through the protected
      `crates-io` environment.
- [ ] A published Python 3.15 prerelease is clean-installed by downstream
      qualification without repository or sibling-checkout leakage.

## Implementation Order

1. Correct metadata, dependency declarations, dependency groups, version
   lookup, and lockfile coverage.
2. Qualify Python 3.15 and dependency minimums locally and on Linux CI.
3. Refactor source, Rust, documentation, and cloud CI into explicit jobs.
4. Add the six-platform native wheel matrix and installed-wheel tests.
5. Add Windows ARM64 full-extra, minimum-dependency, wheel, and sdist tests.
6. Migrate documentation deployment and verify public URLs before deleting
   tracked generated HTML.
7. Add and dry-run the independent Python and Rust release workflows.
8. Configure protected environments, Trusted Publishers, Pages, and required
   checks.
9. Publish a Python prerelease and qualify downstream projects against its
   published artifacts.
10. Publish a stable Python release only after Python 3.15 and its required
    dependency artifacts are stable.

## Out of Scope

- Free-threaded CPython wheels or safety claims.
- Python 3.10 or 3.11 compatibility after this release.
- Musllinux wheels.
- New CPU-specific optimizations or numerical algorithm changes.
- Combining Python and Rust versions or forcing coordinated releases.
- Publishing a real PyPI or crates.io release without explicit maintainer
  authorization.
