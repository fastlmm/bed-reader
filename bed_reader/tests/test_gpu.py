"""GPU reading (bed-reader built with the `cubecl-cuda` Cargo feature): compare `read_gpu` and
`standardize_gpu` with the CPU functions exactly, and check the zero-copy handoff to CuPy.

Skipped, with the reason shown, when bed-reader was built without GPU support or CuPy is missing.
"""

import gc

import numpy as np
import pytest

from bed_reader import open_bed, sample_file, standardize_f32, standardize_f64, to_bed

bed_reader_rust = pytest.importorskip("bed_reader.bed_reader")
if not hasattr(bed_reader_rust, "read_gpu"):
    pytest.skip("bed-reader was built without the cubecl-cuda Cargo feature", allow_module_level=True)
cp = pytest.importorskip("cupy", reason="CuPy is needed to receive GPU arrays")

FILES = [
    "small.bed",
    "smallmode0.bed",
    "some_missing.bed",
    "plink_sim_10s_100v_10pmiss.bed",
    "sparse.bed",
    "toydata.5chrom.bed",
]


def index_cases(count):
    cases = [np.arange(count), np.arange(count)[::-1], np.arange(0, count, 3), np.array([count // 2])]
    cases.append(np.array([-1, 0, count - 1, -count, 0]))
    return [case.astype(np.intp) for case in cases]


def same(a, b):
    """Bitwise equal, with NaN matching NaN."""
    return (
        a.shape == b.shape
        and a.dtype == b.dtype
        and np.array_equal(a, b, equal_nan=True)
        and np.array_equal(np.signbit(a), np.signbit(b))
    )


def cpu_standardize(val, beta):
    stats = np.zeros((val.shape[1], 2), dtype=val.dtype)
    func = standardize_f32 if val.dtype == np.float32 else standardize_f64
    func(val, beta is not None, *(beta or (1.0, 25.0)), True, False, stats, 1)
    return stats


def gpu_read(path, bed, iid_index, sid_index, is_a1_counted, dtype, order):
    return bed_reader_rust.read_gpu(
        str(path), bed.iid_count, bed.sid_count, is_a1_counted, iid_index, sid_index, dtype, order == "F"
    )


def gpu_standardize(array, beta):
    return bed_reader_rust.standardize_gpu(array, beta is not None, *(beta or (1.0, 25.0)))


@pytest.mark.parametrize("name", FILES)
def test_gpu_matches_cpu(name):
    path = sample_file(name)
    with open_bed(path, count_A1=True) as bed_a1, open_bed(path, count_A1=False) as bed_a2:
        for iid_index in index_cases(bed_a1.iid_count):
            for sid_index in index_cases(bed_a1.sid_count):
                for is_a1_counted, bed in [(True, bed_a1), (False, bed_a2)]:
                    for dtype in ["float32", "float64"]:
                        for order in ["F", "C"]:
                            expected = bed.read(index=(iid_index, sid_index), dtype=dtype, order=order)
                            array = gpu_read(path, bed, iid_index, sid_index, is_a1_counted, dtype, order)
                            view = cp.asarray(array)
                            assert view.flags["F_CONTIGUOUS" if order == "F" else "C_CONTIGUOUS"]
                            assert same(cp.asnumpy(view, order="A"), expected)
                            for beta in [None, (1.0, 25.0), (2.0, 3.0)]:
                                cpu_val = expected.copy(order="A")
                                array = gpu_read(path, bed, iid_index, sid_index, is_a1_counted, dtype, order)
                                try:
                                    cpu_stats = cpu_standardize(cpu_val, beta)
                                except ValueError as cpu_error:
                                    with pytest.raises(ValueError) as gpu_error:
                                        gpu_standardize(array, beta)
                                    assert str(gpu_error.value) == str(cpu_error)
                                    continue
                                gpu_stats = gpu_standardize(array, beta)
                                assert same(gpu_stats, cpu_stats)
                                assert same(cp.asnumpy(cp.asarray(array), order="A"), cpu_val)


def test_gpu_single_value_and_all_missing(tmp_path):
    path = tmp_path / "special.bed"
    # Columns 0 and 1 are single-value SNPs; column 2 is all missing.
    val = np.array([[0, 1, -127, 2], [0, 1, -127, 2], [0, 1, -127, 0], [0, 1, -127, 1]], dtype=np.int8)
    to_bed(path, val)
    with open_bed(path) as bed:
        for dtype in ["float32", "float64"]:
            for order in ["F", "C"]:
                sid_index = np.array([0, 1, 3], dtype=np.intp)
                expected = bed.read(index=np.s_[:, sid_index], dtype=dtype, order=order)
                array = gpu_read(path, bed, np.arange(4, dtype=np.intp), sid_index, True, dtype, order)
                assert same(gpu_standardize(array, None), cpu_standardize(expected, None))
                assert same(cp.asnumpy(cp.asarray(array), order="A"), expected)

                all_missing = bed.read(dtype=dtype, order=order)
                with pytest.raises(ValueError) as cpu_error:
                    cpu_standardize(all_missing, None)
                array = gpu_read(path, bed, np.arange(4, dtype=np.intp), np.arange(4, dtype=np.intp), True, dtype, order)
                with pytest.raises(ValueError) as gpu_error:
                    gpu_standardize(array, None)
                assert str(gpu_error.value) == str(cpu_error.value)


def test_gpu_zero_copy_and_lifetime():
    path = sample_file("toydata.5chrom.bed")
    with open_bed(path) as bed:
        expected = bed.read(dtype="float32")
        array = gpu_read(path, bed, np.arange(bed.iid_count, dtype=np.intp), np.arange(bed.sid_count, dtype=np.intp), True, "float32", "F")
    view = cp.asarray(array)
    assert view.data.ptr == array.__cuda_array_interface__["data"][0], "CuPy copied instead of wrapping"
    del array
    gc.collect()
    cp.zeros(10_000_000, dtype=cp.float32)  # churn CuPy's allocator
    assert same(cp.asnumpy(view), expected), "the CuPy view did not keep the GPU buffer alive"
