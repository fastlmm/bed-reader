.. currentmodule:: bed_reader

.. _gpu-reading:

Experimental GPU reads
======================

Python GPU reading currently supports NVIDIA CUDA and CuPy. It requires a
CUDA-enabled bed-reader build and a compatible NVIDIA driver.
:func:`has_gpu_support` checks the installed build; it does not
check whether a GPU is available. Stable releases do not yet include CUDA-enabled
wheels. The ``1.3.0a1`` alpha is being prepared for Linux x86-64 and native
Windows x86-64; WSL2 uses the Linux wheel. Python AMD
support remains deferred; Rust's separate HIP Cargo feature is unaffected.

Installing the alpha
--------------------

Publication is pending. After ``1.3.0a1`` is published, explicitly select it:

.. code-block:: bash

    pip install --only-binary=bed-reader "bed-reader[cupy]==1.3.0a1"
    # Or, in an existing uv environment:
    uv pip install --only-binary=bed-reader "bed-reader[cupy]==1.3.0a1"

The alpha targets CPython 3.11 through 3.15 on Windows and Linux x86-64.
CuPy 14.2 supports Python through 3.14, so Python 3.15 is currently a CPU-only
installation path. Linux alpha wheels require glibc 2.35 or later. Other
architectures and macOS are outside this initial alpha wheel set.
``--only-binary=bed-reader`` avoids building a CPU-only extension from the
source distribution when a matching GPU wheel is unavailable. For CPU-only
use of the alpha, omit ``[cupy]``. Normal unpinned installs continue selecting
the stable release.

Installing an experimental wheel
--------------------------------

If you have a CUDA-enabled wheel built from this checkout, its ``cupy`` extra
installs CuPy's CUDA 13 distribution. It follows CuPy's prerequisites: a
compatible NVIDIA driver and a CUDA 13 toolkit, including NVRTC and CUDA
headers. For example, replace the filename below
with the wheel for your Python version and platform:

.. code-block:: bash

    pip install "./bed_reader-1.3.0a1-cp313-cp313-manylinux_2_35_x86_64.whl[cupy]"
    # Alternatively, install into an existing uv environment:
    uv pip install "./bed_reader-1.3.0a1-cp313-cp313-manylinux_2_35_x86_64.whl[cupy]"

For native Windows, use a ``win_amd64.whl`` file instead. A prebuilt wheel does
not require Rust, Maturin, or a C++ compiler on the user's machine. The NVIDIA
driver and CUDA toolkit remain separate prerequisites. CuPy's installation guide explains its
`CUDA setup <https://docs.cupy.dev/en/stable/install.html>`_.

Initial artifact checks use WSL2 and native Windows with an NVIDIA RTX 3060.
Native Linux GPU testing is still needed; the Linux wheel is currently tested
under WSL2. Additional hardware reports are welcome.

The extra adds Python dependencies; it cannot enable GPU support in a CPU-only
wheel. Do not use ``pip install bed-reader[cupy]`` against the public release
until CUDA-enabled wheels are published. GPU-capable wheels still support CPU
reads without importing CuPy or initializing CUDA. GPU dependencies are not
included in an ordinary ``pip install bed-reader`` installation.

CuPy also offers pip-managed CUDA components for users who prefer them to a
system toolkit. For example, the following installs CuPy with the CUDA 13.3
components used in our initial wheel checks:

.. code-block:: bash

    pip install "cupy-cuda13x[ctk]>=14.2,<15" "cuda-toolkit==13.3.*"

This is an alternative CUDA setup, not an additional bed-reader requirement.
It still requires a compatible NVIDIA driver. See CuPy's installation guide
for supported configurations and avoid installing multiple CuPy distributions
in the same environment.

Building from source
--------------------

Build with the Cargo features ``cubecl-cuda`` (CUDA support) and
``extension-module`` (Python extension). A source build requires Rust 1.95 or
later and a CUDA toolkit providing ``nvcc``. From this checkout, use ``uv`` and
Maturin (omit the clone commands if you already have the repository):

.. code-block:: bash

    git clone https://github.com/fastlmm/bed-reader.git
    cd bed-reader
    git switch gpu-cubecl  # alpha development branch
    uv venv --python 3.13 ~/.cache/bed-reader-gpu/venv
    uv pip install --python ~/.cache/bed-reader-gpu/venv/bin/python numpy pooch
    uv tool run --from maturin maturin build --release --locked \
        --features cubecl-cuda,extension-module \
        --interpreter ~/.cache/bed-reader-gpu/venv/bin/python \
        --target-dir ~/.cache/bed-reader-gpu/target \
        --out ~/.cache/bed-reader-gpu/wheels
    uv pip install --python ~/.cache/bed-reader-gpu/venv/bin/python \
        --reinstall --no-index --no-deps \
        --find-links ~/.cache/bed-reader-gpu/wheels bed-reader

Install the matching CuPy distribution into the same virtual environment before
running this example with ``~/.cache/bed-reader-gpu/venv/bin/python``:

.. code-block:: python

    import cupy as cp
    import numpy as np
    from bed_reader import has_gpu_support, open_bed, sample_file

    if not has_gpu_support():
        raise RuntimeError("Rebuild bed-reader with the cubecl-cuda Cargo feature before requesting CuPy")

    with open_bed(sample_file("small.bed")) as bed:
        values = bed.read(index=np.s_[[2, 0], 1:4], dtype="float32", order="F", xp=cp)
        assert isinstance(values, cp.ndarray)
        assert values.shape == (2, 3)
        means = cp.nanmean(values, axis=0)  # computation stays on the GPU
        host_means = cp.asnumpy(means)  # explicitly copy only the result to the CPU
        cpu_values = bed.read(xp="numpy")  # overrides ARRAY_MODULE=cupy
        assert isinstance(cpu_values, np.ndarray)

``xp`` accepts the NumPy or CuPy module, or the strings ``"numpy"`` and ``"cupy"``.
With ``xp=None``, ``ARRAY_MODULE`` supplies the default; when unset, reads use
NumPy. An explicit ``xp`` overrides the environment. Requesting CuPy fails clearly
if CuPy or native GPU support is missing; there is no CPU decoding fallback.

GPU reads honor ``index``, ``dtype`` (int8, float32, float64), ``order`` (C/F),
and the reader's ``count_A1`` setting, just as dense CPU reads do. Rows are individuals and columns are
SNPs. Missing values are NaN for floats and -127 for int8. Packed input bytes
move to the GPU and are decoded there. The returned ``cupy.ndarray`` shares
Rust-owned memory without copying the decoded matrix; views keep it alive.
Releasing the last view waits for GPU consumers before native memory is recycled.
Use ordinary CuPy operations on the array; ``cp.asnumpy`` explicitly transfers
a result to the CPU.

GPU reads use device 0. They require a local file and reject
``force_python_only=True``. ``num_threads`` is ignored by GPU reads;
``max_concurrent_requests`` and ``max_chunk_bytes`` apply only to CPU cloud reads.
``read_sparse`` remains a CPU/SciPy operation and has no ``xp`` argument. A decoded GPU matrix is limited to ``u32::MAX`` elements
and available device memory. Rust supports separate CUDA and HIP Cargo features;
Python HIP reading is not available yet.

Allocating GPU arrays without a BED file
----------------------------------------

:func:`bed_reader.empty_gpu` allocates an uninitialized matrix through CubeCL
and returns an ordinary CuPy array. It has the same CUDA build requirements,
device 0 selection, dtypes, and element-count limit as GPU reads. It does not
change CuPy's global allocator. Fill the array before reading its values::

    import cupy as cp
    from bed_reader import empty_gpu, is_cubecl_backed

    values = empty_gpu((3, 2), dtype="float32", order="F")
    values.fill(1)
    assert isinstance(values, cp.ndarray)
    assert is_cubecl_backed(values)
    assert values.flags.f_contiguous
    assert float(values.sum()) == 6

:func:`bed_reader.is_cubecl_backed` identifies the underlying allocation without
copying data or synchronizing the GPU. Views retain that backing; results of
CuPy operations can use a different allocator. Check each result rather than
inferring ownership from its ``cupy.ndarray`` type. False means the allocation
is not tracked by bed-reader; it does not identify another allocator.

Views keep the native allocation alive. Owners may move between Python threads;
releasing the last reference waits for outstanding CUDA work before CubeCL can
reuse its memory. Allocation errors propagate rather than selecting another
allocator or moving computation to the CPU.

CubeCL and CuPy keep separate allocation pools. After dropping native-backed
arrays and all their views, call :func:`bed_reader.gpu_memory_cleanup` to return
unused CubeCL storage to CUDA for other libraries. Live arrays remain valid;
CuPy's own pool is unaffected. See that function's complete example.

Choosing float32 or float64
---------------------------

``dtype=None`` (or omitting dtype) uses ``ARRAY_DTYPE`` when set. Valid values
are ``float32`` and ``float64``; an explicit dtype always wins. Without the
environment variable, BED reads keep their float32 default and ``empty_gpu``
keeps its float64 default. These are Python defaults, separate from the Rust
builder's explicit type selection. NumPy/CuPy operations do not consult this
variable themselves. For example::

    from unittest.mock import patch
    from bed_reader import open_bed, sample_file

    with patch.dict("os.environ", {"ARRAY_MODULE": "cupy", "ARRAY_DTYPE": "float32"}):
        with open_bed(sample_file("small.bed")) as bed:
            values = bed.read()  # CuPy, float32
            precise = bed.read(dtype="float64")  # explicit override
        assert values.dtype.name == "float32"
        assert precise.dtype.name == "float64"
    # Both environment settings are restored, including on failure.


Float64 arithmetic can be much slower than float32 on consumer GPUs. Some
data-center GPUs have high float64 throughput; speed depends on the device and
operation. Allocation and BED decoding are not necessarily compute-bound.
Float64 arrays also use twice the storage. NVIDIA documents architecture-specific
`FP32/FP64 throughput ratios <https://docs.nvidia.com/cuda/cuda-programming-guide/05-appendices/compute-capabilities.html>`_.

GPU reads and allocations with float64 emit :class:`bed_reader.GpuFloat64Warning`.
The requested dtype and all defaults remain unchanged. Python's normal warning
filter displays the notice once per call site; the category allows applications
to handle it separately from precision-loss or correctness warnings.

Use ``bed.read(dtype="float32", xp="cupy")`` or
``empty_gpu(shape, dtype="float32")`` when single precision is sufficient.
Keep float64 when accuracy requires it, and compare representative results before
changing precision. CPU float64 arrays do not trigger the warning. An output
dtype is not a promise about all intermediates: GPU standardization may use
float64 reductions for numerical stability even with float32 storage.
