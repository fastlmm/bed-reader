<!-- TODO0 consider deleting this spec once the work below is implemented and released. -->

# Python GPU wheel qualification

The initial Python target is NVIDIA on Linux x86-64 and Windows x86-64.
WSL2 uses the Linux wheel but requires separate runtime validation. This scope
does not restrict future Rust GPU runtime experiments.

The proposed installation extra is `bed-reader[cupy]`, using CuPy's CUDA 13
distribution, following CuPy's existing driver/toolkit prerequisites. Ordinary
CPU installs must continue to require only NumPy.
The same native wheel should serve CPU and GPU callers on the target platforms.
An extra installs dependencies; it does not change a wheel's Cargo features.

## Qualification gates

- Build CUDA-enabled Linux and Windows wheels with the existing locked Rust
  dependencies. Keep generated files and environments outside the repository.
- Install the resulting wheel with both pip and uv in isolated environments,
  outside the source checkout. Inspect the wheel metadata and dependencies.
- Verify CPU imports and reads without CuPy, and check that CPU operations do
  not load NVIDIA libraries. A machine without a driver is a separate test.
- Exercise GPU reads, ownership, allocation, and Unit/Beta standardization,
  including the fractional-input regression tests, using installed artifacts.
- Record runtime dependencies precisely. The current CubeCL CUDA build uses
  LLVM to compile kernels; its optional C++ compiler uses NVRTC and headers.
  CuPy also needs NVRTC for ordinary array operations. CUDA components may
  come from a system toolkit or CuPy's optional `[ctk]` installation path;
  bed-reader does not mandate a separate CUDA installation mechanism. Rust
  and Maturin are build tools, not prerequisites for a prebuilt wheel.
- Qualify the Python versions and Linux compatibility baseline before changing
  the release artifact matrix. A local CPython 3.13 wheel alone is insufficient
  evidence for the complete supported Python/platform matrix.

## Release blockers

<!-- TODO0 qualify both native platforms, the supported Python matrix, and the
Linux compatibility baseline before enabling CUDA in published wheels. -->

The existing release workflow still builds CPU-only wheels. Adding `[cupy]`
alone does not make those artifacts GPU-capable. Do not advertise the new
extra as a published GPU installation until the artifact builds are updated
and validated. No version changes, publication, or release tags are part of
this experiment.

## Initial local results (November 17, 2026)

- Built CPython 3.13 CUDA-enabled wheels for WSL2/Linux
  (`manylinux_2_35_x86_64`, about 42 MB) and native Windows (`win_amd64`).
- Both pip and uv installed the artifacts outside the source tree. The
  `[cupy]` extra selects `cupy-cuda13x>=14.2,<15`; CUDA setup follows CuPy.
- CPU reads passed without CuPy installed. Further checks confirmed CPU reads
  load no CUDA libraries on either platform, even with GPU dependencies present.
  The Linux extension has no direct CUDA shared-library dependency.
- The existing installed-artifact read/API/allocation/fractional Unit/Beta
  regression selection passed 87 tests on each platform using the RTX 3060.
- Native Windows initially lacked NVRTC despite having `nvcc`. Installing
  CuPy's optional `[ctk]` components with `cuda-toolkit==13.3.*` supplied the
  missing runtime libraries in the test environments. This remains an optional
  CuPy setup, not a dependency imposed by bed-reader.
- WSL2 also passed with pip-supplied CUDA components, system toolkit discovery
  disabled, and a fresh CuPy kernel cache.

These checks do not qualify native Linux, machines without NVIDIA drivers,
all supported Python versions, or older glibc baselines. CuPy 14.2 wheels cover
Python through 3.14; Python 3.15 GPU installation is not yet qualified.

Reusable local scripts and Windows artifacts are under
`/mnt/c/deldir/bed-reader-packaging` (`C:\deldir\bed-reader-packaging`). Linux
artifacts, environments, and logs are under
`~/.cache/bed-reader-gpu/packaging`. No existing user environment was modified.
