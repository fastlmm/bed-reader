.. currentmodule:: bed_reader

.. _gpu-reading:

Experimental GPU reads
======================

Python GPU reading currently supports NVIDIA CUDA and CuPy. It requires a
CUDA-enabled bed-reader build, a compatible NVIDIA driver and CUDA toolkit
(including ``nvcc``), Rust 1.95 or later, and a CuPy distribution matching your CUDA
installation. :func:`has_gpu_support` checks the installed build; it does not
check whether a GPU is available. The normal installation commands above do not
enable the optional Cargo GPU feature. GPU wheel distribution and Python AMD
support remain deferred.

Build with the Cargo features ``cubecl-cuda`` (CUDA support) and
``extension-module`` (Python extension). From this checkout, use ``uv`` and
Maturin (omit the clone commands if you already have the repository):

.. code-block:: bash

    git clone https://github.com/fastlmm/bed-reader.git
    cd bed-reader
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
