"""Environment precision defaults preserve explicit choices and legacy defaults."""
import numpy as np
import pytest

from bed_reader import has_gpu_support, open_bed, sample_file
from bed_reader._array_module import array_dtype
from bed_reader._gpu import empty_gpu


@pytest.mark.parametrize("configured", ["float32", "float64"])
def test_read_precision_precedence(monkeypatch, configured):
    monkeypatch.setenv("ARRAY_MODULE", "numpy")
    monkeypatch.setenv("ARRAY_DTYPE", configured)
    with open_bed(sample_file("small.bed")) as bed:
        assert bed.read().dtype == np.dtype(configured)
        assert bed.read(dtype=None).dtype == np.dtype(configured)
        for dtype in (np.float32, np.float64, np.int8):
            assert bed.read(dtype=dtype).dtype == dtype
        assert bed.read_sparse().dtype == np.dtype(configured)


def test_unconfigured_read_default_is_float32(monkeypatch):
    monkeypatch.delenv("ARRAY_DTYPE", raising=False)
    monkeypatch.setenv("ARRAY_MODULE", "numpy")
    assert array_dtype() == np.float64
    assert array_dtype(default="float32") == np.float32
    with open_bed(sample_file("small.bed")) as bed:
        assert bed.read().dtype == np.float32
        assert bed.read(dtype=None).dtype == np.float32
        assert bed.read_sparse().dtype == np.float32


@pytest.mark.parametrize("configured", ["", "float16", "int8", "FLOAT32", "unknown"])
def test_invalid_environment_is_rejected_but_explicit_dtype_wins(monkeypatch, configured):
    monkeypatch.setenv("ARRAY_DTYPE", configured)
    with pytest.raises(ValueError, match="ARRAY_DTYPE"):
        array_dtype()
    assert array_dtype(np.float32) == np.float32
    assert array_dtype(np.int8) == np.int8


@pytest.mark.parametrize("configured", ["float32", "float64"])
def test_gpu_read_and_allocation_defaults(monkeypatch, configured):
    if not has_gpu_support():
        pytest.skip("requires the cubecl-cuda Cargo feature")
    cupy = pytest.importorskip("cupy")
    monkeypatch.setenv("ARRAY_DTYPE", configured)
    with open_bed(sample_file("small.bed")) as bed:
        values = bed.read(xp=cupy)
    assert values.dtype == np.dtype(configured)
    assert empty_gpu((2, 3)).dtype == np.dtype(configured)
    assert empty_gpu((2, 3), dtype="float32").dtype == np.float32
