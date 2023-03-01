import logging
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

from bed_reader import (
    file_aat_piece_f64_orderf,
    file_ata_piece_f64_orderf,
    file_b_less_aatbx,
    file_dot_piece,
)
from bed_reader._open_bed import get_num_threads, open_bed  # noqa


def file_ata(filename, offset, iid_count, sid_count, sid_step):
    ata = np.full((sid_count, sid_count), np.nan)
    for sid_index in range(0, sid_count, sid_step):
        sid_range_len = min(sid_step, sid_count - sid_index)
        ata_piece = np.full((sid_count - sid_index, sid_range_len), np.nan)
        if sid_index % 2 == 0:  # test new and old method
            file_ata_piece_f64_orderf(
                str(filename),
                offset,
                iid_count,
                sid_count,
                sid_index,
                ata_piece,
                num_threads=get_num_threads(None),
                log_frequency=sid_range_len,
            )
        else:
            file_dot_piece(
                str(filename),
                offset,
                iid_count,
                sid_index,
                ata_piece,
                num_threads=get_num_threads(None),
                log_frequency=sid_range_len,
            )
        ata[sid_index:, sid_index : sid_index + sid_range_len] = ata_piece
    for sid_index in range(sid_count):
        ata[sid_index, sid_index + 1 :] = ata[sid_index + 1 :, sid_index]
    return ata


def file_aat(filename, offset, iid_count, sid_count, iid_step):
    aat = np.full((iid_count, iid_count), np.nan)
    for iid_index in range(0, iid_count, iid_step):
        iid_range_len = min(iid_step, iid_count - iid_index)
        aat_piece = np.full((iid_count - iid_index, iid_range_len), np.nan)
        file_aat_piece_f64_orderf(
            str(filename),
            offset,
            iid_count,
            sid_count,
            iid_index,
            aat_piece,
            num_threads=get_num_threads(None),
            log_frequency=iid_range_len,
        )
        aat[iid_index:, iid_index : iid_index + iid_range_len] = aat_piece
    for iid_index in range(iid_count):
        aat[iid_index, iid_index + 1 :] = aat[iid_index + 1 :, iid_index]

    return aat


def write_read_test_file_ata(iid_count, sid_count, sid_step, tmp_path):
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

    out_val = file_ata(file_path, offset, iid_count, sid_count, sid_step)
    expected = mm.T.dot(mm)
    assert np.allclose(expected, out_val, equal_nan=True)


def write_read_test_file_aat(
    iid_count, sid_count, iid_step, tmp_path, skip_if_there=False
):
    offset = 640
    file_path = tmp_path / f"{iid_count}x{sid_count}_o{offset}_array.memmap"
    if skip_if_there and file_path.exists():
        logging.info(f"'{file_path}' exists")
        mm = np.memmap(
            file_path,
            dtype="float64",
            mode="r",
            offset=offset,
            shape=(iid_count, sid_count),
            order="F",
        )
    else:
        logging.info(f"Creating '{file_path}'")
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
        logging.info(f"Finished creating '{file_path}'")

    start_time = datetime.now()
    out_val = file_aat(file_path, offset, iid_count, sid_count, iid_step)
    delta = datetime.now() - start_time
    expected = mm.dot(mm.T)
    assert np.allclose(expected, out_val, equal_nan=True)
    return delta


def test_file_ata_medium(tmp_path):
    write_read_test_file_ata(100, 1000, 33, tmp_path)


def test_file_aat_medium(tmp_path):
    write_read_test_file_aat(100, 1000, 34, tmp_path)


# # Too slow
# def test_file_ata_giant(tmp_path):
#     write_read_test_file_ata(100_000, 10_000, 1000, tmp_path)

# # Too slow
# def test_file_aat_giant(tmp_path):
#     write_read_test_file_ata(1_000, 10_000, 100, tmp_path)


def test_file_ata_small(shared_datadir):
    filename = shared_datadir / "small_array.memmap"

    out_val = file_ata(filename, 0, 2, 3, 2)
    print(out_val)

    expected = np.array([[17.0, 22.0, 27.0], [22.0, 29.0, 36.0], [27.0, 36.0, 45.0]])
    print(expected)
    assert np.allclose(expected, out_val, equal_nan=True)


def test_file_aat_small(shared_datadir):
    filename = shared_datadir / "small_array.memmap"

    out_val = file_aat(filename, 0, iid_count=2, sid_count=3, iid_step=1)
    print(out_val)

    expected = np.array([[14.0, 32.0], [32.0, 77.0]])
    print(expected)
    assert np.allclose(expected, out_val, equal_nan=True)


def mmultfile_b_less_aatb(a_snp_mem_map, b, log_frequency=0, force_python_only=False):
    # Without memory efficiency
    #   a=a_snp_mem_map.val
    #   aTb = np.dot(a.T,b)
    #   aaTb = b-np.dot(a,aTb)
    #   return aTb, aaTb

    if force_python_only:
        aTb = np.zeros(
            (a_snp_mem_map.shape[1], b.shape[1])
        )  # b can be destroyed. Is everything is in best order, i.e. F vs C
        aaTb = b.copy()
        b_mem = np.array(b, order="F")
        with open(a_snp_mem_map.filename, "rb") as U_fp:
            U_fp.seek(a_snp_mem_map.offset)
            for i in range(a_snp_mem_map.shape[1]):
                a_mem = np.fromfile(
                    U_fp, dtype=np.float64, count=a_snp_mem_map.shape[0]
                )
                if log_frequency > 0 and i % log_frequency == 0:
                    logging.info("{0}/{1}".format(i, a_snp_mem_map.shape[1]))
                aTb[i, :] = np.dot(a_mem, b_mem)
                aaTb -= np.dot(a_mem.reshape(-1, 1), aTb[i : i + 1, :])
    else:
        b1 = np.array(b, order="F")
        aTb = np.zeros((a_snp_mem_map.shape[1], b.shape[1]))
        aaTb = np.array(b1, order="F")

        file_b_less_aatbx(
            str(a_snp_mem_map.filename),
            a_snp_mem_map.offset,
            a_snp_mem_map.shape[0],  # row count
            b1,  # B copy 1 in "F" order
            aaTb,  # B copy 2 in "F" order
            aTb,  # result
            num_threads=get_num_threads(None),
            log_frequency=log_frequency,
        )

    return aTb, aaTb


def write_read_test_file_b_less_aatbx(
    iid_count, a_sid_count, b_sid_count, log_frequency, tmp_path, do_both=True
):
    offset = 640
    file_path = tmp_path / f"{iid_count}x{a_sid_count}_o{offset}_array.memmap"
    mm = np.memmap(
        file_path,
        dtype="float64",
        mode="w+",
        offset=offset,
        shape=(iid_count, a_sid_count),
        order="F",
    )
    mm[:] = np.linspace(0, 1, mm.size).reshape(mm.shape)
    mm.flush()

    b = np.array(
        np.linspace(0, 1, iid_count * b_sid_count).reshape(
            (iid_count, b_sid_count), order="F"
        )
    )
    b_again = b.copy()

    logging.info("Calling Rust")
    aTb, aaTb = mmultfile_b_less_aatb(
        mm, b_again, log_frequency, force_python_only=False
    )

    if do_both:
        logging.info("Calling Python")
        aTb_python, aaTb_python = mmultfile_b_less_aatb(
            mm, b, log_frequency, force_python_only=True
        )

        if (
            not np.abs(aTb_python - aTb).max() < 1e-8
            or not np.abs(aaTb_python - aaTb).max() < 1e-8
        ):
            raise AssertionError(
                "Expect Python and Rust to get the same mmultfile_b_less_aatb answer"
            )


def test_file_b_less_aatbx_medium(tmp_path):
    write_read_test_file_b_less_aatbx(500, 400, 100, 10, tmp_path, do_both=True)


def test_file_b_less_aatbx_medium2(tmp_path):
    write_read_test_file_b_less_aatbx(5_000, 400, 100, 100, tmp_path, do_both=True)


# Slow and doesn't check answer
# def test_file_b_less_aatbx_2(tmp_path):
#     write_read_test_file_b_less_aatbx(50_000, 4000, 1000, 100, tmp_path, do_both=False)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    shared_datadir = Path(r"D:\OneDrive\programs\bed-reader\bed_reader\tests\data")
    tmp_path = Path(r"m:/deldir/tests")

    if False:
        small2_array = np.array(
            [[1.0, 2.0, 3.0, 4.0], [5.0, 6.0, 7.0, 8.0], [9.0, 10.0, 11.0, 12.0]],
            order="F",
        )
        mm = np.memmap(
            tmp_path / "small2_array.memmap",
            dtype="float64",
            mode="w+",
            shape=(3, 4),
            order="F",
        )
        mm[:] = small2_array[:]
        print(mm.offset)
        mm.flush()

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

    # test_file_b_less_aatbx_2(tmp_path)
    # test_zero_files(tmp_path)
    # test_index(shared_datadir)
    # test_c_reader_bed(shared_datadir)
    # test_read1(shared_datadir)
    pytest.main([__file__])
