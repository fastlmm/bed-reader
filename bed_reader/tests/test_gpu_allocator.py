"""Allocation validation plus CUDA ownership, views, and stream interoperability."""

import gc
import sys
import weakref
from concurrent.futures import ThreadPoolExecutor
from types import SimpleNamespace

import numpy as np
import pytest

from bed_reader import empty_gpu, gpu_memory_cleanup, is_cubecl_backed, open_bed, sample_file
from bed_reader import bed_reader as native
from bed_reader._gpu import _native_owner


@pytest.fixture
def cupy():
    if not hasattr(native, "allocate_gpu"):
        pytest.skip("requires the cubecl-cuda Cargo feature")
    return pytest.importorskip("cupy")


@pytest.mark.parametrize("shape", [2, (2,), (2, 3, 4), (-1, 2), (1.5, 2), (True, 2), (2**32, 0), (65536, 65536)])
def test_invalid_shape(shape):
    with pytest.raises(ValueError):
        empty_gpu(shape)


@pytest.mark.parametrize("dtype,order", [("float16", "C"), ("int32", "F"), ("float32", "A")])
def test_invalid_type_or_order(dtype, order):
    with pytest.raises(ValueError):
        empty_gpu((2, 3), dtype=dtype, order=order)


def test_missing_native_allocator(monkeypatch):
    monkeypatch.delattr(native, "allocate_gpu", raising=False)
    with pytest.raises(RuntimeError, match="cubecl-cuda Cargo feature"):
        empty_gpu((2, 3))


def test_missing_native_cleanup(monkeypatch):
    monkeypatch.delattr(native, "gpu_memory_cleanup", raising=False)
    with pytest.raises(RuntimeError, match="cubecl-cuda Cargo feature"):
        gpu_memory_cleanup()


def test_ownership_query_without_cupy(monkeypatch):
    monkeypatch.setitem(sys.modules, "cupy", None)
    for value in (None, object(), [[1, 2]], np.ones((2, 3))):
        assert is_cubecl_backed(value) is False


@pytest.mark.parametrize("shape", [(4, 3), (0, 3), (4, 0)])
@pytest.mark.parametrize("order", ["C", "F"])
def test_ownership_query_views_and_copies(cupy, shape, order):
    values = empty_gpu(shape, dtype="float32", order=order)
    values.fill(7)
    owner = weakref.ref(_native_owner(values))
    views = (values[1:, 1:], values[::-1], values.T, values[:0], values.reshape(-1, order=order))
    assert is_cubecl_backed(values) is True
    for view in views:
        assert is_cubecl_backed(view) is True
    assert is_cubecl_backed(cupy.zeros((4, 3), dtype="float32")) is False
    assert is_cubecl_backed(np.ones((4, 3))) is False
    if values.size:
        assert is_cubecl_backed(values.copy()) is False
        assert is_cubecl_backed(values + 1) is False
    del values
    gc.collect()
    # The query recognizes offset views after the complete matrix is gone.
    assert is_cubecl_backed(views[0]) is True
    del views, view
    gc.collect()
    assert owner() is None  # querying must not extend the owner's lifetime


def test_gpu_read_ownership(cupy):
    with open_bed(sample_file("small.bed")) as bed:
        values = bed.read(dtype="float32", xp=cupy)
    assert is_cubecl_backed(values) is True
    assert is_cubecl_backed(values[1:, 1:]) is True
    assert is_cubecl_backed(cupy.asnumpy(values)) is False


@pytest.mark.parametrize("dtype", [np.int8, np.float32, np.float64])
@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize("shape", [(4, 3), (0, 3), (4, 0)])
def test_allocate_and_fill(cupy, dtype, order, shape):
    values = empty_gpu(shape, dtype=dtype, order=order)
    assert isinstance(values, cupy.ndarray) and values.dtype == dtype
    assert values.shape == shape
    assert values.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"]
    values.fill(7)
    np.testing.assert_array_equal(cupy.asnumpy(values), np.full(shape, 7, dtype=dtype))


def test_views_keep_owner_without_registry_leak(cupy):
    values = empty_gpu((8, 3), dtype="float32")
    values.fill(7)
    owner = weakref.ref(_native_owner(values))
    view = values[:, 1:]
    assert _native_owner(view) is None  # never standardize the larger owner through a slice
    del values
    gc.collect()
    assert owner() is not None
    for _ in range(4):
        other = empty_gpu((8, 3), dtype="float32")
        other.fill(9)
    gpu_memory_cleanup()  # a live view must keep its allocation protected
    np.testing.assert_array_equal(cupy.asnumpy(view), 7)
    del view
    gc.collect()
    assert owner() is None
    gpu_memory_cleanup()


@pytest.mark.parametrize("order", ["C", "F"])
def test_foreign_stream_writes_then_cubecl_then_foreign_read(cupy, order):
    values = empty_gpu((4, 2), dtype="float64", order=order)
    expected = np.array([[0., 1.], [1., 1.], [2., 1.], [1., 1.]])
    writer = cupy.cuda.Stream(non_blocking=True)
    with writer:
        values[...] = cupy.asarray(expected)
    native.standardize_gpu(_native_owner(values), False, 1., 25.)
    expected[:, 0] = (expected[:, 0] - 1) / np.sqrt(0.5)
    expected[:, 1] = 0
    reader = cupy.cuda.Stream(non_blocking=True)
    with reader:
        result = values.copy()
    reader.synchronize()
    np.testing.assert_allclose(cupy.asnumpy(result), expected, rtol=4 * np.finfo(float).eps, atol=0)


def test_async_release_before_pool_reuse(cupy):
    stream = cupy.cuda.Stream(non_blocking=True)
    values = empty_gpu((1024, 32), dtype="float32")
    with stream:
        values.fill(3)
        result = values.copy()
    del values  # must complete the foreign work before CubeCL can recycle storage
    gc.collect()
    replacement = empty_gpu((1024, 32), dtype="float32")
    replacement.fill(9)
    stream.synchronize()
    np.testing.assert_array_equal(cupy.asnumpy(result), 3)


def test_native_owner_can_be_used_and_released_on_another_thread(cupy):
    values = empty_gpu((3, 2), dtype="float64")
    values[...] = cupy.array([[0., 1.], [1., 1.], [2., 1.]])
    owner = weakref.ref(_native_owner(values))

    def consume(array):
        native.standardize_gpu(_native_owner(array), False, 1., 25.)
        return cupy.asnumpy(array)

    with ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(consume, values)
        del values
        result = future.result()
    # Joining the worker also ensures its last argument reference was released.
    gc.collect()
    assert owner() is None
    np.testing.assert_allclose(result[:, 0], [-np.sqrt(1.5), 0, np.sqrt(1.5)],
                               rtol=4 * np.finfo(float).eps, atol=0)
    np.testing.assert_array_equal(result[:, 1], 0)


def test_int8_standardization_still_rejected(cupy):
    values = empty_gpu((3, 2), dtype="int8")
    values.fill(1)
    with pytest.raises(ValueError, match="float32 or float64"):
        native.standardize_gpu(_native_owner(values), False, 1., 25.)


def test_allocator_does_not_change_active_device(cupy, monkeypatch):
    monkeypatch.setattr(cupy.cuda, "Device", lambda: SimpleNamespace(id=1))
    with pytest.raises(ValueError, match="device 0"):
        empty_gpu((2, 3))
