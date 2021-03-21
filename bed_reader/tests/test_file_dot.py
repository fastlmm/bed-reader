import logging
import os
import platform
from pathlib import Path

import numpy as np
import pytest

from bed_reader import file_dot_piece
from bed_reader._open_bed import get_num_threads, open_bed  # noqa


def file_dot(filename, offset, iid_count, sid_count, sid_step):
    ata = np.zeros((sid_count, sid_count))
    for sid_start in range(0, sid_count, sid_step):
        sid_range_len = min(sid_step, sid_count - sid_start)
        ata_piece = np.zeros((sid_count - sid_start, sid_range_len))
        file_dot_piece(
            str(filename),
            offset,
            iid_count,
            sid_start,
            ata_piece,
            num_threads=get_num_threads(None),
            log_frequency=sid_range_len,
        )
        ata[sid_start:, sid_start : sid_start + sid_range_len] = ata_piece
    return ata


def test_file_dot_small(shared_datadir):

    filename = shared_datadir / "small_array.memmap"

    out_val = file_dot(filename, 0, 2, 3, 2)
    print(out_val)

    expected = np.array([[17.0, 0, 0], [22.0, 29.0, 0], [27.0, 36.0, 45.0]])
    print(expected)
    assert np.allclose(expected, out_val, equal_nan=True)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    shared_datadir = Path(r"D:\OneDrive\programs\bed-reader\bed_reader\tests\data")
    tmp_path = Path(r"m:/deldir/tests")

    if False:
        small_array = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], order="F")
        mm = np.memmap(
            tmp_path / "small_array.memmap",
            dtype="float64",
            mode="w+",
            shape=(2, 3),
            order="F",
        )
        mm[:] = small_array[:]
        print(mm.offset)
        mm.flush()

    if False:
        mm = np.memmap(
            tmp_path / "100x1000_o640_array.memmap",
            dtype="float64",
            mode="w+",
            offset=640,
            shape=(100, 1000),
            order="F",
        )
        total = mm.shape[0] * mm.shape[1]
        lin = np.linspace(0, 1, total).reshape(mm.shape)

        mm[:] = lin[:]
        print(mm.offset)
        mm.flush()

    if False:
        mm = np.memmap(
            tmp_path / "1000x10000_o640_array.memmap",
            dtype="float64",
            mode="w+",
            offset=640,
            shape=(1000, 10_000),
            order="F",
        )
        total = mm.shape[0] * mm.shape[1]
        lin = np.linspace(0, 1, total).reshape(mm.shape)

        mm[:] = lin[:]
        print(mm.offset)
        mm.flush()

    if False:
        mm = np.memmap(
            tmp_path / "10_000x100_000_o6400_array.memmap",
            dtype="float64",
            mode="w+",
            offset=640,
            shape=(10_000, 100_000),
            order="F",
        )
        total = mm.shape[0] * mm.shape[1]
        lin = np.linspace(0, 1, total).reshape(mm.shape)

        mm[:] = lin[:]
        print(mm.offset)
        mm.flush()

    if False:
        mm = np.memmap(
            tmp_path / "100_000x10_000_o6400_array.memmap",
            dtype="float64",
            mode="w+",
            offset=640,
            shape=(100_000, 10_000),
            order="F",
        )
        total = mm.shape[0] * mm.shape[1]
        lin = np.linspace(0, 1, total).reshape(mm.shape)

        mm[:] = lin[:]
        print(mm.offset)
        mm.flush()

    # test_zero_files(tmp_path)
    # test_index(shared_datadir)
    # test_c_reader_bed(shared_datadir)
    # test_read1(shared_datadir)
    pytest.main([__file__])
