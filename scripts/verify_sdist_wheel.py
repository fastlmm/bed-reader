"""Compare wheels built from a clean Git tree and its source distribution.

Both builds use the same unpacked source path and SOURCE_DATE_EPOCH so that
native compiler path metadata and archive timestamps cannot conceal a source
distribution mismatch. This validates complete wheel bytes, including binaries.
"""
import argparse
import hashlib
import io
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path


def build(package, destination, interpreter, environment):
    subprocess.run(
        ["uv", "build", "--wheel", "--no-sources", "--python", str(interpreter),
         "--out-dir", str(destination)],
        cwd=package, env=environment, check=True,
    )
    wheels = list(destination.glob("*.whl"))
    if len(wheels) != 1:
        raise RuntimeError(f"Expected exactly one wheel in {destination}: {wheels}")
    return wheels[0]


def verify(source, sdist, interpreter):
    subprocess.run(["git", "diff", "--quiet"], cwd=source, check=True)
    subprocess.run(["git", "diff", "--cached", "--quiet"], cwd=source, check=True)
    epoch = subprocess.check_output(
        ["git", "show", "-s", "--format=%ct", "HEAD"], cwd=source, text=True,
    ).strip()
    archive = subprocess.check_output(["git", "archive", "HEAD"], cwd=source)
    environment = dict(os.environ, SOURCE_DATE_EPOCH=epoch)
    # A target directory outside the source would reuse the first build and
    # could hide a missing source file. Each build must start from fresh output.
    environment.pop("CARGO_TARGET_DIR", None)
    with tempfile.TemporaryDirectory(prefix="bed-reader-repro-") as directory:
        root = Path(directory)
        package = root / "package"
        package.mkdir()
        with tarfile.open(fileobj=io.BytesIO(archive)) as contents:
            contents.extractall(package, filter="data")
        direct = build(package, root / "direct", interpreter, environment)
        shutil.rmtree(package)
        unpacked = root / "unpacked"
        with tarfile.open(sdist) as contents:
            contents.extractall(unpacked, filter="data")
        roots = list(unpacked.iterdir())
        if len(roots) != 1 or not roots[0].is_dir():
            raise RuntimeError("Source distribution must contain one package directory")
        roots[0].rename(package)
        derived = build(package, root / "derived", interpreter, environment)
        if direct.name != derived.name or direct.read_bytes() != derived.read_bytes():
            with zipfile.ZipFile(direct) as first, zipfile.ZipFile(derived) as second:
                paths = sorted(set(first.namelist()) | set(second.namelist()))
                mismatches = [name for name in paths if name not in first.namelist()
                              or name not in second.namelist()
                              or first.read(name) != second.read(name)]
            raise AssertionError(
                f"Direct and sdist-derived wheels differ: {direct.name}, {derived.name}; "
                f"payload differences: {mismatches}. Archive-only differences also fail."
            )
        print(f"Identical wheels: {direct.name}; SHA-256={hashlib.sha256(direct.read_bytes()).hexdigest()}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sdist", type=Path)
    parser.add_argument("--source", type=Path, default=Path.cwd())
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    args = parser.parse_args()
    verify(args.source.resolve(), args.sdist.resolve(), args.python.resolve())


if __name__ == "__main__":
    main()
