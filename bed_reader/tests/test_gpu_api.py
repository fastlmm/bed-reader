"""Public array-module selection, including CPU-only build diagnostics."""
import gc
from unittest.mock import patch

import numpy as np
import pytest

from bed_reader import has_gpu_support, open_bed, sample_file


def test_explicit_numpy_overrides_environment():
    with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
        assert isinstance(open_bed(sample_file("small.bed")).read(xp=np), np.ndarray)


def test_missing_gpu_build_raises():
    with patch("bed_reader.has_gpu_support", return_value=False):
        with pytest.raises(RuntimeError, match="cubecl-cuda"):
            open_bed(sample_file("small.bed"))._read_gpu()


@pytest.mark.parametrize("dtype", ["int8", "float32", "float64"])
@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize("count_A1", [True, False])
def test_public_gpu_read_and_view_lifetime(dtype, order, count_A1):
    if not has_gpu_support():
        pytest.skip("requires bed-reader built with the cubecl-cuda Cargo feature")
    cp = pytest.importorskip("cupy")
    bed = open_bed(sample_file("small.bed"), count_A1=count_A1)
    index = np.s_[::-1, [3, 0, 3]]
    expected = bed.read(index=index, dtype=dtype, order=order, xp=np)
    actual = bed.read(index=index, dtype=dtype, order=order, xp=cp,
                      num_threads=1, max_concurrent_requests=1, max_chunk_bytes=1)
    assert actual.dtype == expected.dtype
    assert actual.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"]
    np.testing.assert_array_equal(cp.asnumpy(actual), expected)
    view = actual[::2, :]
    del actual
    gc.collect()
    cp.zeros((300, 100))
    np.testing.assert_array_equal(cp.asnumpy(view), expected[::2, :])


def test_gpu_default_and_unsupported_options():
    if not has_gpu_support():
        pytest.skip("requires bed-reader built with the cubecl-cuda Cargo feature")
    cp = pytest.importorskip("cupy")
    with open_bed(sample_file("small.bed")) as bed:
        with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
            assert isinstance(bed.read(), cp.ndarray)
        with pytest.raises(ValueError, match="force_python_only"):
            bed.read(xp=cp, force_python_only=True)
    # Skip opening remote data: GPU requests must reject the URL before any read.
    with open_bed("https://example.com/cohort.bed", skip_format_check=True) as bed:
        with pytest.raises(ValueError, match="cloud"):
            bed.read(xp=cp)


def test_release_waits_for_nondefault_stream():
    if not has_gpu_support():
        pytest.skip("requires bed-reader built with the cubecl-cuda Cargo feature")
    cp = pytest.importorskip("cupy")
    bed = open_bed(sample_file("toydata.5chrom.bed"))
    values = bed.read(xp=cp)
    stream = cp.cuda.Stream(non_blocking=True)
    event = cp.cuda.Event()
    with stream:
        result = values @ values.T
        event.record()
    del values
    # The allocation can now return to CubeCL's pool. Consumer work must be
    # complete even though it ran on a nondefault CuPy stream.
    assert event.done
    assert result.shape == (bed.iid_count, bed.iid_count)


@pytest.mark.parametrize("index", [np.s_[:, []], np.s_[[], :], np.s_[[], []]])
def test_empty_gpu_selection(index):
    if not has_gpu_support():
        pytest.skip("requires bed-reader built with the cubecl-cuda Cargo feature")
    cp = pytest.importorskip("cupy")
    bed = open_bed(sample_file("small.bed"))
    actual = bed.read(index=index, xp=cp)
    expected = bed.read(index=index, xp=np)
    np.testing.assert_array_equal(cp.asnumpy(actual), expected)
