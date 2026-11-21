"""Environment precision defaults preserve explicit choices and legacy defaults."""
import numpy as np
import pytest

from bed_reader import open_bed, sample_file
from bed_reader._array_module import array_dtype


@pytest.mark.parametrize("configured", ["float32", "float64"])
@pytest.mark.parametrize("force_python_only", [False, True])
def test_read_precision_precedence(monkeypatch, configured, force_python_only):
    monkeypatch.setenv("ARRAY_DTYPE", configured)
    with open_bed(sample_file("small.bed")) as bed:
        assert bed.read(force_python_only=force_python_only).dtype == np.dtype(configured)
        assert bed.read(dtype=None, force_python_only=force_python_only).dtype == np.dtype(configured)
        for dtype in (np.float32, np.float64, np.int8):
            assert bed.read(dtype=dtype, force_python_only=force_python_only).dtype == dtype
        assert bed.read_sparse().dtype == np.dtype(configured)
        assert bed.read_sparse(dtype=None).dtype == np.dtype(configured)


def test_unconfigured_read_default_is_float32(monkeypatch):
    monkeypatch.delenv("ARRAY_DTYPE", raising=False)
    assert array_dtype() == np.float64
    assert array_dtype(default="float32") == np.float32
    with open_bed(sample_file("small.bed")) as bed:
        assert bed.read().dtype == np.float32
        assert bed.read(dtype=None).dtype == np.float32
        assert bed.read_sparse().dtype == np.float32
        assert bed.read_sparse(dtype=None).dtype == np.float32


@pytest.mark.parametrize("configured", ["", "float16", "int8", "FLOAT32", "unknown"])
def test_invalid_environment_is_rejected_but_explicit_dtype_wins(monkeypatch, configured):
    monkeypatch.setenv("ARRAY_DTYPE", configured)
    with pytest.raises(ValueError, match="ARRAY_DTYPE"):
        array_dtype()
    assert array_dtype(np.float32) == np.float32
    assert array_dtype(np.int8) == np.int8

    with open_bed(sample_file("small.bed")) as bed:
        for read in (bed.read, bed.read_sparse):
            with pytest.raises(ValueError, match="ARRAY_DTYPE"):
                read()
            assert read(dtype="float32").dtype == np.float32
