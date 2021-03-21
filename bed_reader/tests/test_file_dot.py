import logging
import os
import platform
from pathlib import Path

import numpy as np
import pytest

from bed_reader import file_dot_piece
from bed_reader._open_bed import get_num_threads, open_bed  # noqa


def file_dot(filename, offset, iid_count, sid_count, sid_step):
    ata = np.full((sid_count, sid_count), np.nan)
    for sid_start in range(0, sid_count, sid_step):
        sid_range_len = min(sid_step, sid_count - sid_start)
        ata_piece = np.full((sid_count - sid_start, sid_range_len), np.nan)
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
    for sid_index in range(sid_count):
        ata[sid_index, sid_index + 1 :] = ata[sid_index + 1 :, sid_index]
    return ata


def write_read_test(iid_count, sid_count, sid_step, tmp_path):
    offset = 640
    file_path = tmp_path / f"{iid_count}x{sid_count}_o{offset}_array.memmap"
    mm = np.memmap(
        file_path,
        dtype="float64",
        mode="w+",
        offset=offset,
        shape=(iid_count, sid_count),
        order="F",
    )
    mm[:] = np.linspace(0, 1, mm.size).reshape(mm.shape)
    mm.flush()

    out_val = file_dot(file_path, offset, iid_count, sid_count, sid_step)
    expected = mm.T.dot(mm)
    assert np.allclose(expected, out_val, equal_nan=True)


def test_file_dot_medium(tmp_path):
    write_read_test(100, 1000, 33, tmp_path)


def test_file_dot_small(shared_datadir):

    filename = shared_datadir / "small_array.memmap"

    out_val = file_dot(filename, 0, 2, 3, 2)
    print(out_val)

    expected = np.array([[17.0, 22.0, 27.0], [22.0, 29.0, 36.0], [27.0, 36.0, 45.0]])
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
