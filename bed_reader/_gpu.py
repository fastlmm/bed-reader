"""CubeCL allocation and ownership bookkeeping for ordinary CuPy arrays."""

from operator import index
from weakref import WeakValueDictionary
import sys
import warnings

import numpy as np

from bed_reader._array_module import array_dtype


# CuPy retains the native owner; this registry must not extend its lifetime.
_owners = WeakValueDictionary()


class GpuFloat64Warning(UserWarning):
    """Float64 GPU storage or computation may have a substantial performance cost.

    Float64 remains supported and is never silently converted. See :doc:`gpu`
    for precision choices and a float32 example.
    """


def _warn_gpu_float64(dtype):
    if np.dtype(dtype) != np.dtype("float64"):
        return
    # Attribute nested reader/workflow warnings to the user's call. Python's
    # normal warning filter then shows this message once per user call site.
    frame = sys._getframe(1)
    stacklevel = 2
    while frame is not None:
        modules = frame.f_globals.get("__name__", "").split(".")
        # pandas can consume our lazy map/reduce results; it is an intermediate
        # frame, not the user's precision choice or allocation call site.
        if modules[0] not in {"bed_reader", "pysnptools", "fastlmm", "pandas"} or any(
            part in {"test", "tests"} or part.startswith("test_") for part in modules[1:]
        ):
            break
        stacklevel += 1
        frame = frame.f_back
    del frame
    warnings.warn(
        "float64 arithmetic can be much slower than float32 on consumer GPUs. Float64 GPU arrays "
        "also use twice the storage. Consider float32 where its precision is "
        "sufficient; the requested dtype is unchanged.",
        GpuFloat64Warning,
        stacklevel=stacklevel,
    )


def _as_cupy(owner):
    import cupy as cp

    array = cp.asarray(owner)
    _owners[array.data.ptr] = owner
    return array


def is_cubecl_backed(array):
    """Return whether a CuPy array uses a bed-reader CubeCL allocation.

    Recognizes arrays returned by :func:`empty_gpu` and GPU reads, including
    their slices, transposes, and reshaped views. This describes the underlying
    allocation, not whether a view supports native standardization. A copy or
    calculation can allocate different storage; check the result separately.

    Returns False for NumPy arrays, other objects, and GPU allocations not
    tracked by bed-reader. Does not require CuPy to be installed, copy values,
    synchronize the GPU, or keep the allocation alive after returning.

    Parameters
    ----------
    array : object
        The array whose allocation to inspect.

    Returns
    -------
    bool
        Whether the array is backed by a live bed-reader CubeCL allocation.
        False does not identify which other allocator owns the memory.

    Examples
    --------
    With a CUDA-enabled build and a working GPU::

        import cupy as cp
        from bed_reader._gpu import empty_gpu, is_cubecl_backed

        values = empty_gpu((3, 2), dtype="float32")
        values.fill(1)
        assert is_cubecl_backed(values)
        assert is_cubecl_backed(values[1:, :])
        assert not is_cubecl_backed(cp.array(values, copy=True))
    """
    # An existing CuPy ndarray implies CuPy is already imported. CPU callers
    # need neither an optional import nor CUDA initialization for this query.
    cupy = sys.modules.get("cupy")
    if cupy is None or not isinstance(array, cupy.ndarray):
        return False
    if array.device.id != 0:
        return False
    # Views can have an offset data pointer. The public memory object's pointer
    # still identifies the original allocation registered by _as_cupy. The
    # registry is weak; CuPy itself retains the native allocation's owner.
    return _owners.get(array.data.mem.ptr) is not None


def empty_gpu(shape, dtype=None, order="F"):
    """Allocate an uninitialized 2-D CuPy array backed by CubeCL memory.

    Requires CuPy, CUDA device 0, and bed-reader built with the ``cubecl-cuda``
    Cargo feature. Supports int8, float32, float64, and C/F order. Each dimension
    and the element count must fit in uint32. Allocation errors propagate;
    there is no CPU fallback. Values are uninitialized: fill them before use.

    The array and its views retain their native owner, including across Python
    threads. Array destruction waits for outstanding CUDA consumers before
    recycling memory. CubeCL pools freed allocations separately from CuPy;
    :func:`gpu_memory_cleanup` releases its unused cached storage.

    Parameters
    ----------
    shape : pair of int
        Nonnegative row and column counts.
    dtype : data-type, optional
        int8, float32, or float64. None uses ARRAY_DTYPE, then float64;
        an explicit dtype overrides the environment.
        Float64 emits GpuFloat64Warning about potentially slow arithmetic on
        consumer GPUs. Choose float32 when its precision is sufficient.
    order : {'F', 'C'}, optional
        Output layout, default F (column major).

    Returns
    -------
    cupy.ndarray
        An ordinary CuPy array with CubeCL-owned storage on device 0.

    Examples
    --------
    With a CUDA-enabled build and a working GPU::

        import cupy as cp
        from bed_reader._gpu import empty_gpu

        values = empty_gpu((3, 2), dtype="float32", order="F")
        values.fill(1)
        assert isinstance(values, cp.ndarray)
        assert values.flags.f_contiguous
        assert float(values.sum()) == 6
    """
    try:
        dimensions = tuple(shape)
    except TypeError as error:
        raise ValueError("shape must contain two nonnegative integers") from error
    if len(dimensions) != 2 or any(isinstance(value, (bool, np.bool_)) for value in dimensions):
        raise ValueError("shape must contain two nonnegative integers")
    try:
        rows, columns = (index(value) for value in dimensions)
    except TypeError as error:
        raise ValueError("shape must contain two nonnegative integers") from error
    if min(rows, columns) < 0:
        raise ValueError("shape must contain two nonnegative integers")
    if max(rows, columns, rows * columns) > np.iinfo(np.uint32).max:
        raise ValueError("GPU dimensions and element count must fit in uint32")
    dtype = array_dtype(dtype)
    if dtype not in (np.dtype("int8"), np.dtype("float32"), np.dtype("float64")):
        raise ValueError("GPU dtype must be int8, float32, or float64")
    if order not in ("C", "F"):
        raise ValueError("GPU order must be 'F' or 'C'")
    from . import bed_reader as native

    if not hasattr(native, "allocate_gpu"):
        raise RuntimeError("GPU allocation requires bed-reader built with the cubecl-cuda Cargo feature")
    _require_device_zero()
    _warn_gpu_float64(dtype)
    return _as_cupy(native.allocate_gpu(rows, columns, dtype.name, order == "F"))


def _require_device_zero():
    """Raise if CuPy's current device is not CUDA device 0, which GPU reads and allocations use.

    Without this check, a caller inside ``with cupy.cuda.Device(1):`` would silently get an
    array on device 0 and only see a CuPy device-mismatch error later.
    """
    try:
        import cupy as cp
    except ImportError:
        return  # Without CuPy there is no current device to disagree with device 0.
    device = cp.cuda.Device().id
    if device != 0:
        raise ValueError(
            f"bed-reader GPU arrays use CUDA device 0, but CuPy's current device is {device}. "
            "To use another GPU, start Python with CUDA_VISIBLE_DEVICES set to it, "
            "which makes it device 0."
        )


def gpu_memory_cleanup():
    """Release all unused GPU memory cached by bed-reader's CubeCL runtime.

    Internal: arrays from GPU reads and :func:`empty_gpu` already return their
    memory to CUDA when their last view is dropped. This also frees CubeCL's
    small pooled temporaries. Live arrays remain valid. CuPy's own pool is
    unaffected. Requires a GPU-enabled build; otherwise raises RuntimeError.
    CUDA errors propagate.

    Examples
    --------
    With a GPU-enabled build and a working GPU::

        import cupy as cp
        from bed_reader import open_bed, sample_file
        from bed_reader._gpu import gpu_memory_cleanup

        with open_bed(sample_file("small.bed")) as bed:
            values = bed.read(xp=cp)
        del values
        gpu_memory_cleanup()
    """
    from . import bed_reader as native

    if not hasattr(native, "gpu_memory_cleanup"):
        raise RuntimeError("GPU memory cleanup requires bed-reader built with the cubecl-cuda Cargo feature")
    native.gpu_memory_cleanup()


def _native_owner(array):
    """Return the owner only for an unchanged complete native matrix view."""
    if array.device.id != 0:
        return None
    owner = _owners.get(array.data.ptr)
    if owner is None:
        return None
    interface = owner.__cuda_array_interface__
    if (array.data.ptr == interface["data"][0]
            and array.shape == interface["shape"]
            and array.strides == interface["strides"]
            and array.dtype == np.dtype(interface["typestr"])):
        return owner
    return None
