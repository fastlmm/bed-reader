"""Floating-point CubeCL standardization, with independent centered references."""
import numpy as np
import pytest
from scipy.stats import beta as beta_distribution

from bed_reader import has_gpu_support, open_bed, sample_file, standardize_f32, to_bed
from bed_reader import bed_reader as native


@pytest.fixture
def cupy():
    if not has_gpu_support():
        pytest.skip("requires bed-reader built with the cubecl-cuda Cargo feature")
    return pytest.importorskip("cupy")


def reference(values, beta=None):
    """Population moments in extended precision, then round statistics to the API dtype."""
    stats = np.empty((values.shape[1], 2), dtype=values.dtype)
    result = np.zeros_like(values)
    for column in range(values.shape[1]):
        observed = values[:, column][~np.isnan(values[:, column])].astype(np.longdouble)
        assert len(observed) > 0
        mean = observed.mean()
        variance = ((observed - mean) ** 2).mean()
        stats[column] = [mean, np.sqrt(variance) if variance > 0 else np.inf]
        if variance == 0:
            continue
        mean, deviation = stats[column]
        if beta is None:
            factor = values.dtype.type(1) / deviation
        else:
            maf = float(mean) / 2
            factor = values.dtype.type(beta_distribution.pdf(min(maf, 1 - maf), *beta))
        valid = ~np.isnan(values[:, column])
        result[valid, column] = (values[valid, column] - mean) * factor
    return result, stats


def dosage_values(dtype):
    return np.array([
        [0.25, 0.75, 0, 2, np.nan, 0.5],
        [0.5, 0.75, 0, 2, np.nan, 1.25],
        [np.nan, 0.75, 0, 2, 1.125, 1.5],
        [1.75, 0.75, 0, 2, np.nan, 0.875],
        [1.125, 0.75, 0, 2, np.nan, 1.5],
        [1.5, 0.75, 0, 2, np.nan, 0.25],
    ], dtype=dtype)


def owner_with_values(path, values, order, cupy):
    to_bed(path, np.zeros(values.shape, dtype=np.int8))
    owner = open_bed(path)._read_gpu(dtype=values.dtype, order=order)
    view = cupy.asarray(owner)
    # Native standardization must wait for writes from any CuPy stream.
    with cupy.cuda.Stream(non_blocking=True):
        view[...] = cupy.asarray(values)
    return owner, view


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize("beta", [None, (1, 25), (2, 3), (0.5, 0.5)])
def test_fractional_unit_beta(tmp_path, cupy, dtype, order, beta):
    values = dosage_values(dtype)
    if beta == (2, 3):
        # Beta restricts the column mean, not each individual observation.
        values[:, 0] = [-1.25, 0.25, np.nan, 2.5, 1.125, 1.5]
    owner, actual = owner_with_values(tmp_path / "fractional.bed", values, order, cupy)
    expected, expected_stats = reference(values, beta)
    stats = native.standardize_gpu(owner, beta is not None, *(beta or (1, 25)))
    # These binary-exact inputs have exact shifted sums. Allow 16 ulps worth
    # of relative/absolute rounding for sqrt, Beta PDF, reciprocal, and apply.
    tolerance = 16 * np.finfo(dtype).eps
    np.testing.assert_allclose(stats, expected_stats, rtol=tolerance, atol=tolerance)
    np.testing.assert_allclose(cupy.asnumpy(actual), expected, rtol=tolerance, atol=tolerance)
    assert actual.dtype == values.dtype
    assert actual.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"]


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("order", ["C", "F"])
def test_unit_negative_out_of_range_and_cancellation(tmp_path, cupy, dtype, order):
    values = dosage_values(dtype)
    values[:, 0] -= 4
    values[:, 1] += 3
    values[:, 2] = -3.25
    # E[x*x]-E[x]**2 loses this variance in both dtypes. The increments and
    # mean are binary-exact, so the centered reference needs no large tolerance.
    values[:, 3] = (1e4 if dtype == np.float32 else 1e12) + np.arange(6) * 0.125
    owner, actual = owner_with_values(tmp_path / "outside.bed", values, order, cupy)
    expected, expected_stats = reference(values)
    stats = native.standardize_gpu(owner, False, 1, 25)
    tolerance = 16 * np.finfo(dtype).eps
    np.testing.assert_allclose(stats, expected_stats, rtol=tolerance, atol=tolerance)
    np.testing.assert_allclose(cupy.asnumpy(actual), expected, rtol=tolerance, atol=tolerance)
    assert np.isfinite(stats[3, 1]) and stats[3, 1] > 0


def test_general_errors_and_int8(tmp_path, cupy):
    values = dosage_values(np.float64)
    values[:, 0] = np.nan
    owner, actual = owner_with_values(tmp_path / "missing.bed", values, "F", cupy)
    with pytest.raises(ValueError, match="No individual observed"):
        native.standardize_gpu(owner, False, 1, 25)
    values[:, 0] = -1.25
    owner, actual = owner_with_values(tmp_path / "beta_mean.bed", values, "F", cupy)
    with pytest.raises(ValueError, match="Illegal SNP mean"):
        native.standardize_gpu(owner, True, 1, 25)
    owner = open_bed(tmp_path / "beta_mean.bed")._read_gpu(dtype="int8")
    with pytest.raises(ValueError, match="not int8"):
        native.standardize_gpu(owner, False, 1, 25)


def test_fused_read_after_standardization_errors(cupy):
    # The preceding error cases can leave a decode queued without exporting its
    # allocation. A subsequent packed upload must respect CubeCL's buffer reuse.
    with open_bed(sample_file("toydata.5chrom.bed")) as bed:
        rows = np.arange(bed.iid_count)[::-3][::-1]
        columns = np.arange(bed.sid_count)[::97][::-1]
        expected = bed.read(index=(rows, columns), dtype="float32", xp=np)
        expected_stats = np.zeros((len(columns), 2), dtype=np.float32)
        standardize_f32(expected, False, 1.0, 25.0, True, False, expected_stats, 1)
        owner = bed._read_gpu(index=(rows, columns), dtype="float32")
        stats = native.standardize_gpu(owner, False, 1.0, 25.0)
        np.testing.assert_array_equal(stats, expected_stats)
        np.testing.assert_array_equal(cupy.asnumpy(cupy.asarray(owner)), expected)


def test_invalid_variance_is_runtime_error(tmp_path, cupy):
    # Finite observations whose centered squares overflow must retain the native diagnostic.
    values = np.array([[1e308], [-1e308]], dtype=np.float64)
    owner, actual = owner_with_values(tmp_path / "variance_overflow.bed", values, "F", cupy)
    with pytest.raises(RuntimeError, match="Invalid GPU variance:"):
        native.standardize_gpu(owner, False, 1, 25)
