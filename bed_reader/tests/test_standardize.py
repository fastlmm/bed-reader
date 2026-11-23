"""The native standardizer zeros constant SNPs, including fractional ones."""
import numpy as np
import pytest

from bed_reader import standardize_f32, standardize_f64


@pytest.mark.parametrize(
    "dtype, constant, standardize",
    [(np.float64, 0.1, standardize_f64), (np.float32, 1.7, standardize_f32)],
)
@pytest.mark.parametrize("beta", [False, True])
@pytest.mark.parametrize("order", ["F", "C"])
def test_constant_fractional_snp(dtype, constant, standardize, beta, order):
    # A rounded mean of 1,000 equal values leaves E[x*x] - E[x]**2 a tiny positive value,
    # which used to give nonzero results and, when training, a finite standard deviation.
    # Column 0 is the constant SNP (with one missing value); column 1 varies.
    iid_index = np.arange(1_001)
    val = np.empty((len(iid_index), 2), dtype=dtype, order=order)
    val[:, 0] = constant
    val[500, 0] = np.nan
    val[:, 1] = (iid_index % 3) * dtype(constant)
    stats = np.zeros((2, 2), dtype=dtype, order=order)
    standardize(val, beta, 1.0, 25.0, True, False, stats, 1)
    assert np.all(val[:, 0] == 0)
    assert stats[0, 1] == np.inf
    assert np.any(val[:, 1] != 0)
    assert np.isfinite(stats[1, 1])

    # Applying the trained statistics to new values of the SNP gives zero, too.
    new_val = np.empty_like(val)
    new_val[:, :] = ((iid_index % 5) * dtype(constant))[:, np.newaxis]
    standardize(new_val, beta, 1.0, 25.0, True, True, stats, 1)
    assert np.all(new_val[:, 0] == 0)
    assert np.any(new_val[:, 1] != 0)


@pytest.mark.parametrize("order", ["F", "C"])
def test_float32_std_of_rare_variant_counted_by_common_allele(order):
    # Values mostly 2 made float32 sums cancel in E[x*x] - E[x]**2 (relative std error up to
    # about 3e-4); sums are now accumulated in float64.
    rng = np.random.default_rng(0)
    genotype = 2.0 - rng.binomial(2, 0.0002, size=(14_925, 50))
    genotype = genotype[:, genotype.std(axis=0) > 0]
    exact = genotype.std(axis=0)  # float64, two-pass
    stats = np.zeros((genotype.shape[1], 2), dtype=np.float32)
    standardize_f32(np.array(genotype, dtype=np.float32, order=order), False, 1.0, 1.0, False, False, stats, 1)
    relative_error = np.abs(stats[:, 1] - exact) / exact
    # Conversion of the final std to float32 alone allows half an ulp (3e-8).
    assert relative_error.max() < 4 * np.finfo(np.float32).eps
