import pytest  # noqa: D100

from bed_reader import open_bed, sample_file


def test_optional_dependencies(shared_datadir) -> None:  # noqa: ANN001, D103
    try:
        from scipy import sparse
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


def test_sample_download_retries_transient_failure(tmp_path) -> None:  # noqa: ANN001
    """A transient network error on the first attempt is retried."""
    pooch = pytest.importorskip("pooch")
    requests = pytest.importorskip("requests")
    from bed_reader import _sample_data

    store = pooch.create(
        path=tmp_path,
        base_url=_sample_data.POOCH.base_url,
        registry=_sample_data.POOCH.registry,
        retry_if_failed=_sample_data.POOCH.retry_if_failed,
    )
    attempts = []

    def flaky(url, output_file, pooch_obj) -> None:  # noqa: ANN001
        attempts.append(url)
        if len(attempts) == 1:
            msg = "simulated name-resolution failure"
            raise requests.exceptions.ConnectionError(msg)
        pooch.HTTPDownloader()(url, output_file, pooch_obj)

    path = store.fetch("small.fam", downloader=flaky)
    assert len(attempts) == 2
    assert path.endswith("small.fam")
