import pytest

from bed_reader import open_bed, sample_file


def test_optional_dependencies(shared_datadir):
    try:
        import scipy.sparse as sparse
    except ImportError:
        sparse = None

    file = shared_datadir / "plink_sim_10s_100v_10pmiss.bed"
    with open_bed(file) as bed:
        if sparse is None:
            with pytest.raises(ImportError):
                _ = bed.read_sparse(dtype="int8")
        else:
            _ = bed.read_sparse(dtype="int8")

    try:
        import pooch
    except ImportError:
        pooch = None

    if pooch is None:
        with pytest.raises(ImportError):
            _ = sample_file("small.bed")
    else:
        _ = sample_file("small.bed")
