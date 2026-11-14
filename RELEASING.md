# Releasing bed-reader

Python and Rust releases share a version in `pyproject.toml` and `Cargo.toml`.
The maintainer chooses the version and performs or explicitly authorizes package
publication. Preparing this workflow does not authorize creating a release tag.

## Qualify the candidate

1. Update both versions and the changelog, refresh `uv.lock` and `Cargo.lock`,
   and review dependency API changes. Version 1.2 requires Python 3.11 and Rust
   1.85. Its cloud-file migration requires Rust callers to use cloud-file 0.3.
2. Open a PR against `main`. Require all CI jobs: Rust tests with and without
   default features, pinned-Rust strict Clippy, source Python tests on 3.11–3.15,
   native wheels for Linux/macOS/Windows on x86-64 and ARM64, installed wheel
   tests on 3.11 and 3.15, and source-distribution tests on both boundary versions.
   Do not waive a platform because its runners are queued. Branch protection
   requires the PR to be up to date with `main` and to pass the `CI result`
   check, so CI does not run again after merging; it reruns monthly and for
   the release tag.
3. Require metadata validation and complete wheel comparison between the Git
   source and sdist. `scripts/verify_sdist_wheel.py` uses the same source path,
   timestamp, and fresh build output for both builds, including native binaries.
4. Build strict Sphinx documentation from the frozen environment:

   ```console
   uv sync --frozen --all-extras --python 3.14
   uv run --frozen --no-sync sphinx-build -W --keep-going -E -a -b html doc/source doc/build/html
   ```

5. From the clean candidate commit, run the Rust package validation:

   ```console
   cargo publish --dry-run --locked
   ```

   Inspect the packaged file list and license, and verify that this builds the
   packaged crate rather than depending on unpublished sibling checkouts.

## Publish after maintainer approval

Before creating a tag, configure PyPI trusted publishing for owner `fastlmm`,
repository `bed-reader`, workflow `release.yml`, and environment `pypi`. Configure
that GitHub environment with the maintainer as a required reviewer. Review the
release workflow before authorizing it: pushing `v1.2.0` starts full qualification
and then the PyPI publishing job.

After qualification, update the changelog date and merge the reviewed candidate.
The maintainer may then create and push the matching version tag. The release
workflow checks Python/Rust/tag agreement, runs the canonical CI workflow, and
publishes the exact `wheels-all` artifacts only after every qualification job
passes and the `pypi` environment is approved. It does not rebuild after approval
or silently skip existing package files.

The Rust crate is published separately by the maintainer, from the same clean
tagged source, with `cargo publish --locked`. Never substitute a local cloud-file
checkout for the published dependency during final acceptance.

## Verify published prerequisites

Download and install the published PyPI wheels in clean environments on the
supported boundary interpreters and architectures, and verify the Rust crate
from crates.io. Only then raise PySnpTools/FaST-LMM's bed-reader minimum to the
new version and run their final acceptance. GitHub artifacts are useful candidate
evidence; they do not establish that the published dependency graph works.
