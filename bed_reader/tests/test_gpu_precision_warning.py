"""FP64 notices preserve requested precision and point to the calling code."""

from pathlib import Path
import warnings

import numpy as np
import pytest

from bed_reader import GpuFloat64Warning, empty_gpu, has_gpu_support, open_bed, sample_file


@pytest.fixture
def cupy():
    if not has_gpu_support():
        pytest.skip("requires the cubecl-cuda Cargo feature")
    return pytest.importorskip("cupy")


def test_float64_allocation_warns_once_and_keeps_dtype(cupy):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        values = empty_gpu((3, 2), dtype="float64")
    notices = [item for item in caught if issubclass(item.category, GpuFloat64Warning)]
    assert len(notices) == 1
    assert Path(notices[0].filename).resolve() == Path(__file__).resolve()
    assert str(notices[0].message).startswith("float64 arithmetic can be much slower")
    assert values.dtype == np.float64


def test_read_warning_and_cpu_float32_silence(cupy):
    with open_bed(sample_file("small.bed")) as reader:
        with pytest.warns(GpuFloat64Warning, match="consumer GPUs"):
            values = reader.read(dtype="float64", xp=cupy)
        assert values.dtype == np.float64
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            reader.read(dtype="float32", xp=cupy)
            reader.read(dtype="int8", xp=cupy)
            reader.read(dtype="float64", xp=np)
            empty_gpu((3, 2), dtype="float32")
        assert not any(issubclass(item.category, GpuFloat64Warning) for item in caught)


def test_warning_can_be_selected_without_changing_dtype(cupy):
    with warnings.catch_warnings():
        warnings.simplefilter("error", GpuFloat64Warning)
        with pytest.raises(GpuFloat64Warning):
            empty_gpu((3, 2))


def test_warning_location_skips_pandas_bridge():
    from bed_reader._gpu import _warn_gpu_float64

    namespace = {"__name__": "pandas.core.reshape.concat", "warn": _warn_gpu_float64}
    exec(compile("def consume():\n    warn('float64')\n", "pandas_bridge.py", "exec"), namespace)
    with pytest.warns(GpuFloat64Warning) as caught:
        namespace["consume"]()
    assert caught[0].filename == __file__
