.. currentmodule:: bed_reader

.. _gpu-reading:

GPU reads (beta)
================

bed-reader can read genotype data directly into GPU memory as a
`CuPy <https://cupy.dev/>`_ array. The packed BED bytes are copied to the GPU
and decoded there. GPU reads are optional: without them, bed-reader reads into
NumPy arrays exactly as before.

GPU support is in beta. It currently requires an NVIDIA GPU.

Installing
----------

Install the beta with the ``cupy`` extra:

.. code-block:: bash

    pip install --only-binary=bed-reader "bed-reader[cupy]>=1.3.0b1"
    # Or, in a uv environment:
    uv pip install --only-binary=bed-reader "bed-reader[cupy]>=1.3.0b1"

Naming the beta version is needed while it is a pre-release; an ordinary
``pip install bed-reader`` still selects the latest stable release.
``--only-binary=bed-reader`` requires a prebuilt wheel instead of a source
build. Check the platform requirements below: older Linux systems can still
select a CPU-only wheel.

GPU wheels are available for:

* Linux x86-64 (glibc 2.34 or later), including WSL2.
* Windows x86-64.
* Python 3.11 through 3.15. GPU reads on Python 3.15 still require CuPy
  support for that Python version.

macOS and ARM64 get CPU-only wheels. Linux x86-64 also has CPU-only
manylinux2014 wheels for older glibc; pip selects the GPU wheel on glibc 2.34
or later, even on machines without a GPU. To check an installation, call
:func:`has_gpu_support`; it reports whether the installed build includes
GPU reading, not whether a GPU is present.

bed-reader's GPU reads themselves need only the NVIDIA driver, not the CUDA
toolkit. The ``cupy`` extra installs CuPy's CUDA 13 distribution, which needs a
compatible NVIDIA driver and the CUDA 13 toolkit; see CuPy's
`installation guide <https://docs.cupy.dev/en/stable/install.html>`_.
Alternatively, CuPy can install the CUDA components with pip, which still
requires an NVIDIA driver:

.. code-block:: bash

    pip install "cupy-cuda13x[ctk]>=14.2,<15" "cuda-toolkit==13.3.*"

Avoid installing more than one CuPy distribution in the same environment.

Reading into GPU memory
-----------------------

Pass ``xp="cupy"`` (or the ``cupy`` module) to :meth:`open_bed.read`:

.. code-block:: python

    import cupy as cp
    import numpy as np
    from bed_reader import open_bed, sample_file

    with open_bed(sample_file("small.bed")) as bed:
        values = bed.read(index=np.s_[[2, 0], 1:4], dtype="float32", xp=cp)
    # Confirm where the data is: a cupy.ndarray is always in GPU memory.
    print(type(values), values.device)  # <class 'cupy.ndarray'> <CUDA Device 0>
    assert cp.get_array_module(values) is cp  # numpy would mean CPU memory
    means = cp.nanmean(values, axis=0)  # computation stays on the GPU
    host_means = cp.asnumpy(means)  # copy only the result to the CPU
    assert isinstance(host_means, np.ndarray)

The result is an ordinary CuPy array; use CuPy operations on it, and
``cp.asnumpy`` to copy results to the CPU.
When the array and its views are dropped, bed-reader returns its GPU memory to
CUDA, where CuPy and other libraries can use it.

Requesting CuPy fails with an error if CuPy or GPU support is missing; there
is no silent fallback to the CPU. GPU reads have these limits:

* They use GPU device 0.
* They read local files only, not cloud URLs.
* They do not support ``force_python_only=True``. ``num_threads``,
  ``max_concurrent_requests``, and ``max_chunk_bytes`` do not apply.
* A result can have at most 2³²−1 elements, subject to GPU memory.
* :meth:`open_bed.read_sparse` remains CPU-only.

When to read into GPU memory
----------------------------

A GPU read takes about as long as a well-threaded CPU read. Its advantage is
skipping the copy to the GPU. If you will compute on the GPU, reading float32
data directly is about 2.3 times faster than reading on the CPU and then
copying (about 1.5 times for int8). If you will compute on the CPU, use a CPU
read.

These figures are from reading blocks of 1,000 to 4,000 SNPs for 250,000
individuals from a local file already in the operating system's cache, on an
NVIDIA RTX 3060. When reading from a slow drive, the drive dominates and all
methods take about the same time.

Setting defaults with ARRAY_MODULE and ARRAY_DTYPE
---------------------------------------------------

When ``xp`` or ``dtype`` is omitted (or ``None``), :meth:`open_bed.read` takes
its default from the environment variables ``ARRAY_MODULE`` (``numpy`` or
``cupy``) and ``ARRAY_DTYPE`` (``float32`` or ``float64``). An explicit
argument always wins. Without these variables, reads return NumPy float32
arrays, as before.

Because these are environment variables, a setting applies to everything a
program calls, including libraries such as PySnpTools and FaST-LMM, without
passing arguments through each layer. Scope a setting to a block of code with
``unittest.mock.patch.dict``, which restores the previous values afterward:

.. code-block:: python

    from unittest.mock import patch
    from bed_reader import open_bed, sample_file

    with patch.dict("os.environ", {"ARRAY_MODULE": "cupy", "ARRAY_DTYPE": "float32"}):
        with open_bed(sample_file("small.bed")) as bed:
            values = bed.read()  # CuPy, float32
            precise = bed.read(dtype="float64")  # explicit dtype wins
            on_cpu = bed.read(xp="numpy")  # explicit xp wins
    assert values.dtype.name == "float32"
    assert precise.dtype.name == "float64"
    assert type(on_cpu).__module__ == "numpy"

Choosing float32 or float64
---------------------------

On consumer GPUs, float64 arithmetic is often much slower than float32; some
data-center GPUs have high float64 throughput. NVIDIA documents the
`FP32/FP64 throughput ratios <https://docs.nvidia.com/cuda/cuda-programming-guide/05-appendices/compute-capabilities.html>`_
of each architecture. Float64 arrays also use twice the memory.

A float64 GPU read emits :class:`GpuFloat64Warning` as a reminder. The read
still returns float64; nothing is converted. Python shows the warning once per
call site, and its category lets an application filter it separately from
other warnings.

Use ``dtype="float32"`` when single precision is sufficient, and keep float64
when accuracy requires it. Compare representative results before switching.
