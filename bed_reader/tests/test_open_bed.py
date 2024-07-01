import logging
import os
import platform
from pathlib import Path

import numpy as np
import pytest

from bed_reader import create_bed, open_bed, subset_f64_f64, to_bed


def test_read1(shared_datadir):
    import math

    file = shared_datadir / "plink_sim_10s_100v_10pmiss.bed"
    with open_bed(file) as bed:
        assert bed.iid_count == 10
        assert bed.fid[-1] == "0"
        assert bed.iid[-1] == "9"
        assert bed.shape == (10, 100)

        val = bed.read(dtype="int8")
        # really shouldn't do mean on data where -127 represents missing
        assert val.mean() == -13.142
        val_sparse = bed.read_sparse(dtype="int8")
        assert math.isclose(val_sparse.mean(), -13.142, rel_tol=1e-9)
        assert bed.chromosome[-1] == "1"
        assert bed.bp_position[-1] == 100


def test_write(tmp_path, shared_datadir):
    in_file = shared_datadir / "plink_sim_10s_100v_10pmiss.bed"
    out_file = tmp_path / "out.bed"
    with open_bed(in_file) as bed:
        val0 = bed.read()
        properties0 = {
            "fid": bed.fid,
            "iid": bed.iid,
            "sid": bed.sid,
            "chromosome": bed.chromosome,
            "cm_position": bed.cm_position,
            "bp_position": bed.bp_position,
        }
        to_bed(out_file, val0, properties=properties0)
        with open_bed(out_file) as bed1:
            assert np.allclose(val0, bed1.read(), equal_nan=True)
            val_sparse = bed1.read_sparse()
            assert np.allclose(val0, val_sparse.toarray(), equal_nan=True)
            assert np.array_equal(bed1.fid, properties0["fid"])
            assert np.array_equal(bed1.iid, properties0["iid"])
            assert np.array_equal(bed1.sid, properties0["sid"])
            assert np.issubdtype(bed1.sid.dtype, np.str_)
            assert np.array_equal(bed1.chromosome, properties0["chromosome"])
            assert np.allclose(bed1.cm_position, properties0["cm_position"])
            assert np.allclose(bed1.bp_position, properties0["bp_position"])

    val_float = val0.astype("float")
    val_float[0, 0] = 0.5

    for force_python_only in [False, True]:
        with pytest.raises(ValueError):
            to_bed(
                out_file,
                val_float,
                properties=properties0,
                force_python_only=force_python_only,
            )
    val0[np.isnan(val0)] = 0  # set any nan to 0
    val_int8 = val0.astype("int8")
    val_int8[0, 0] = -1
    for force_python_only in [False, True]:
        with pytest.raises(ValueError):
            to_bed(
                out_file,
                val_int8,
                properties=properties0,
                force_python_only=force_python_only,
            )


def test_write_stream(tmp_path, shared_datadir):
    in_file = shared_datadir / "plink_sim_10s_100v_10pmiss.bed"
    out_file = tmp_path / "out.bed"
    with open_bed(in_file) as bed:
        val0 = bed.read()
        properties0 = {
            "fid": bed.fid,
            "iid": bed.iid,
            "sid": bed.sid,
            "chromosome": bed.chromosome,
            "cm_position": bed.cm_position,
            "bp_position": bed.bp_position,
        }
        with create_bed(
            out_file,
            iid_count=bed.iid_count,
            sid_count=bed.sid_count,
            properties=properties0,
        ) as bed_writer:
            for sid_index in range(bed.sid_count):
                bed_writer.write(val0[:, sid_index])
        with open_bed(out_file) as bed1:
            assert np.allclose(val0, bed1.read(), equal_nan=True)
            val_sparse = bed1.read_sparse()
            assert np.allclose(val0, val_sparse.toarray(), equal_nan=True)
            assert np.array_equal(bed1.fid, properties0["fid"])
            assert np.array_equal(bed1.iid, properties0["iid"])
            assert np.array_equal(bed1.sid, properties0["sid"])
            assert np.issubdtype(bed1.sid.dtype, np.str_)
            assert np.array_equal(bed1.chromosome, properties0["chromosome"])
            assert np.allclose(bed1.cm_position, properties0["cm_position"])
            assert np.allclose(bed1.bp_position, properties0["bp_position"])

    val_float = val0.astype("float")
    val_float[0, 0] = 0.5

    for force_python_only in [False, True]:
        with pytest.raises(ValueError):
            to_bed(
                out_file,
                val_float,
                properties=properties0,
                force_python_only=force_python_only,
            )
    val0[np.isnan(val0)] = 0  # set any nan to 0
    val_int8 = val0.astype("int8")
    val_int8[0, 0] = -1
    for force_python_only in [False, True]:
        with pytest.raises(ValueError):
            to_bed(
                out_file,
                val_int8,
                properties=properties0,
                force_python_only=force_python_only,
            )


def test_overrides(shared_datadir):
    with open_bed(shared_datadir / "some_missing.bed") as bed:
        fid = bed.fid
        iid = bed.iid
        father = bed.father
        mother = bed.mother
        sex = bed.sex
        pheno = bed.pheno
        chromosome = bed.chromosome
        sid = bed.sid
        cm_position = bed.cm_position
        bp_position = bed.bp_position
        allele_1 = bed.allele_1
        allele_2 = bed.allele_2
    # lock in the expected results:
    # np.savez(
    #     shared_datadir / "some_missing.properties.npz",
    #     fid=fid,
    #     iid=iid,
    #     father=father,
    #     mother=mother,
    #     sex=sex,
    #     pheno=pheno,
    #     chromosome=chromosome,
    #     sid=sid,
    #     cm_position=cm_position,
    #     bp_position=bp_position,
    #     allele_1=allele_1,
    #     allele_2=allele_2,
    # )
    property_dict = np.load(shared_datadir / "some_missing.properties.npz")
    assert np.array_equal(property_dict["fid"], fid)
    assert np.array_equal(property_dict["iid"], iid)
    assert np.array_equal(property_dict["father"], father)
    assert np.array_equal(property_dict["mother"], mother)
    assert np.array_equal(property_dict["sex"], sex)
    assert np.array_equal(property_dict["pheno"], pheno)
    assert np.array_equal(property_dict["chromosome"], chromosome)
    assert np.array_equal(property_dict["sid"], sid)
    assert np.array_equal(property_dict["cm_position"], cm_position)
    assert np.array_equal(property_dict["bp_position"], bp_position)
    assert np.array_equal(property_dict["allele_1"], allele_1)
    assert np.array_equal(property_dict["allele_2"], allele_2)

    with pytest.raises(KeyError):
        open_bed(shared_datadir / "some_missing.bed", properties={"unknown": [3, 4, 4]})
    with open_bed(
        shared_datadir / "some_missing.bed", properties={"iid": None}
    ) as bed1:
        assert bed1.iid is None
    with open_bed(shared_datadir / "some_missing.bed", properties={"iid": []}) as bed1:
        assert np.issubdtype(bed1.iid.dtype, np.str_)
        assert len(bed1.iid) == 0
        with pytest.raises(ValueError):
            bed1.father

    with open_bed(
        shared_datadir / "some_missing.bed",
        properties={"sid": [i for i in range(len(sid))]},
    ) as bed1:
        assert np.issubdtype(bed1.sid.dtype, np.str_)
        assert bed1.sid[0] == "0"
    with pytest.raises(ValueError):
        open_bed(
            shared_datadir / "some_missing.bed",
            properties={"sex": ["F" for i in range(len(sex))]},
        )  # Sex must be coded as a number
    with open_bed(
        shared_datadir / "some_missing.bed",
        properties={"sid": np.array([i for i in range(len(sid))])},
    ) as bed1:
        assert np.issubdtype(bed1.sid.dtype, np.str_)
        assert bed1.sid[0] == "0"
    with pytest.raises(ValueError):
        open_bed(
            shared_datadir / "some_missing.bed",
            properties={"sid": np.array([(i, i) for i in range(len(sid))])},
        )
    with open_bed(
        shared_datadir / "some_missing.bed", properties={"sid": [1, 2, 3]}
    ) as bed1:
        with pytest.raises(ValueError):
            bed1.chromosome


def test_str(shared_datadir):
    with open_bed(shared_datadir / "some_missing.bed") as bed:
        assert "open_bed(" in str(bed)


def test_bad_bed(shared_datadir):
    with pytest.raises(ValueError):
        open_bed(shared_datadir / "badfile.bed")
    open_bed(shared_datadir / "badfile.bed", skip_format_check=True)


def test_bad_dtype_or_order(shared_datadir):
    with pytest.raises(ValueError):
        open_bed(shared_datadir / "some_missing.bed").read(dtype=np.int32)
    with pytest.raises(ValueError):
        open_bed(shared_datadir / "some_missing.bed").read(order="X")
    with pytest.raises(ValueError):
        open_bed(shared_datadir / "some_missing.bed").read_sparse(dtype=np.int32)


def setting_generator(seq_dict, seed=9392):
    import itertools

    from numpy.random import RandomState

    longest = max((len(value_list) for value_list in seq_dict.values()))

    for test_index in range(longest):
        setting = {}
        for offset, (key, value_list) in enumerate(seq_dict.items()):
            val = value_list[(test_index + offset) % len(value_list)]
            if not (isinstance(val, str) and "leave_out" == val):
                setting[key] = val
        yield setting

    all_combo = list(itertools.product(*seq_dict.values()))

    random_state = RandomState(seed)
    random_state.shuffle(all_combo)
    for combo in all_combo:
        setting = {
            key: value_list
            for key, value_list in itertools.zip_longest(seq_dict, combo)
            if not (isinstance(value_list, str) and "leave_out" == value_list)
        }
        yield setting


def test_properties(shared_datadir):
    file = shared_datadir / "plink_sim_10s_100v_10pmiss.bed"
    with open_bed(file) as bed:
        iid_list = bed.iid.tolist()
        sid_list = bed.sid.tolist()
        chromosome_list = bed.chromosome.tolist()

    test_count = 75

    seq_dict = {
        "iid": ["leave_out", None, iid_list, np.array(iid_list)],
        "iid_count": ["leave_out", len(iid_list)],
        "iid_before_read": [False, True],
        "iid_after_read": [False, True],
        "sid": ["leave_out", None, sid_list, np.array(sid_list)],
        "sid_count": [None, len(sid_list)],
        "sid_before_read": [False, True],
        "sid_after_read": [False, True],
        "chromosome": ["leave_out", None, chromosome_list, np.array(chromosome_list)],
        "chromosome_before_read": [False, True],
        "chromosome_after_read": [False, True],
    }

    def _not_set_to_none(settings, key):
        return key not in settings or settings[key] is not None

    for test_index, settings in enumerate(setting_generator(seq_dict)):
        if test_index >= test_count:
            break
        with open_bed(
            file,
            iid_count=settings.get("iid_count"),
            sid_count=settings.get("sid_count"),
            properties={
                k: v for k, v in settings.items() if k in {"iid", "sid", "chromosome"}
            },
        ) as bed:
            logging.info(f"Test {test_count}")
            if settings["iid_before_read"]:
                if _not_set_to_none(settings, "iid"):
                    assert np.array_equal(bed.iid, iid_list)
                else:
                    assert bed.iid is None
            if settings["sid_before_read"]:
                if _not_set_to_none(settings, "sid"):
                    assert np.array_equal(bed.sid, sid_list)
                else:
                    assert bed.sid is None
            if settings["chromosome_before_read"]:
                if _not_set_to_none(settings, "chromosome"):
                    assert np.array_equal(bed.chromosome, chromosome_list)
                else:
                    assert bed.chromosome is None
            val = bed.read()
            assert val.shape == (
                len(iid_list),
                len(sid_list),
            )
            val_sparse = bed.read_sparse()
            assert np.allclose(val, val_sparse.toarray(), equal_nan=True)
            if settings["iid_after_read"]:
                if _not_set_to_none(settings, "iid"):
                    assert np.array_equal(bed.iid, iid_list)
                else:
                    assert bed.iid is None
            if settings["sid_after_read"]:
                if _not_set_to_none(settings, "sid"):
                    assert np.array_equal(bed.sid, sid_list)
                else:
                    assert bed.sid is None
            if settings["chromosome_after_read"]:
                if _not_set_to_none(settings, "chromosome"):
                    assert np.array_equal(bed.chromosome, chromosome_list)
                else:
                    assert bed.chromosome is None
            # bed._assert_iid_sid_chromosome()


def test_c_reader_bed(shared_datadir):
    for force_python_only, format in [(False, "csc"), (True, "csr")]:
        bed = open_bed(shared_datadir / "some_missing.bed", count_A1=False)

        val = bed.read(order="F", force_python_only=force_python_only)
        assert val.dtype == np.float32
        ref_val = reference_val(shared_datadir)
        ref_val = ref_val * -1 + 2
        assert np.allclose(ref_val, val, rtol=1e-05, atol=1e-05, equal_nan=True)

        val_sparse = bed.read_sparse(format=format)
        assert val_sparse.dtype == np.float32
        assert np.allclose(
            ref_val, val_sparse.toarray(), rtol=1e-05, atol=1e-05, equal_nan=True
        )

        val = bed.read(order="F", dtype="int8", force_python_only=False)
        assert val.dtype == np.int8
        ref_val[ref_val != ref_val] = -127
        ref_val = ref_val.astype("int8")
        assert np.all(ref_val == val)

        del bed

        with open_bed(shared_datadir / "some_missing.bed") as bed:
            val = bed.read(
                order="F", dtype="float64", force_python_only=force_python_only
            )
            ref_val = reference_val(shared_datadir)
            assert np.allclose(ref_val, val, rtol=1e-05, atol=1e-05, equal_nan=True)
            val_sparse = bed.read_sparse(dtype="float64")
            assert np.allclose(
                ref_val, val_sparse.toarray(), rtol=1e-05, atol=1e-05, equal_nan=True
            )


def reference_val(shared_datadir):
    val = np.load(shared_datadir / "some_missing.val.npy")
    return val


def test_bed_int8(tmp_path, shared_datadir):
    with open_bed(shared_datadir / "some_missing.bed") as bed:
        for force_python_only in [True, False]:
            for order, format in [("F", "csc"), ("C", "csr")]:
                val = bed.read(
                    dtype="int8", force_python_only=force_python_only, order=order
                )
                assert val.dtype == np.int8
                assert (val.flags["C_CONTIGUOUS"] and order == "C") or (
                    val.flags["F_CONTIGUOUS"] and order == "F"
                )
                ref_val = reference_val(shared_datadir)
                ref_val[ref_val != ref_val] = -127
                ref_val = ref_val.astype("int8")
                assert np.array_equal(ref_val, val)
                output = str(tmp_path / "int8.bed")
                for count_A1 in [False, True]:
                    for major in ["individual", "SNP"]:
                        to_bed(
                            output,
                            ref_val,
                            count_A1=count_A1,
                            major=major,
                            force_python_only=force_python_only,
                        )
                        with open_bed(output, count_A1=count_A1) as bed2:
                            val2 = bed2.read(
                                dtype="int8", force_python_only=force_python_only
                            )
                            assert np.array_equal(val2, ref_val)
                            val_sparse = bed2.read_sparse(dtype="int8", format=format)
                            assert np.allclose(val_sparse.toarray(), ref_val)


def test_write1_bed_f64cpp(tmp_path, shared_datadir):
    with open_bed(shared_datadir / "some_missing.bed") as bed:
        for iid_index in [0, 1, 5]:
            for force_python_only, format in [(False, "csc"), (True, "csr")]:
                val_sparse = bed.read_sparse(
                    np.s_[0:iid_index, :], dtype=np.float64, format=format
                )
                assert val_sparse.shape == (iid_index, 100)
                val = bed.read(
                    np.s_[0:iid_index, :],
                    order="F",
                    dtype=np.float64,
                    force_python_only=force_python_only,
                )
                assert val.shape == (iid_index, 100)
                output = str(tmp_path / f"toydata.F64cpp.{iid_index}.bed")
                to_bed(output, val, count_A1=False)
                val2 = open_bed(output, count_A1=False).read(dtype="float64")
                assert np.allclose(val, val2, equal_nan=True)
            assert np.allclose(val_sparse.toarray(), val2, equal_nan=True)


def test_write1_bed_f64cpp_stream(tmp_path, shared_datadir):
    with open_bed(shared_datadir / "some_missing.bed") as bed:
        for iid_index in [0, 1, 5]:
            for force_python_only, format in [(False, "csc"), (True, "csr")]:
                val_sparse = bed.read_sparse(
                    np.s_[0:iid_index, :], dtype=np.float64, format=format
                )
                assert val_sparse.shape == (iid_index, 100)
                val = bed.read(
                    np.s_[0:iid_index, :],
                    order="F",
                    dtype=np.float64,
                    force_python_only=force_python_only,
                )
                assert val.shape == (iid_index, 100)
                output = str(tmp_path / f"toydata.F64cpp.{iid_index}.bed")
                with create_bed(
                    output, iid_count=iid_index, sid_count=100, count_A1=False
                ) as bed_writer:
                    for sid_index in range(100):
                        bed_writer.write(val[:, sid_index])
                val2 = open_bed(output, count_A1=False).read(dtype="float64")
                assert np.allclose(val, val2, equal_nan=True)
            assert np.allclose(val_sparse.toarray(), val2, equal_nan=True)


def test_write1_x_x_cpp(tmp_path, shared_datadir):
    for count_A1 in [False, True]:
        with open_bed(shared_datadir / "some_missing.bed", count_A1=count_A1) as bed:
            for order, format in [("F", "csc"), ("C", "csr")]:
                for dtype in [np.float32, np.float64]:
                    val = bed.read(order=order, dtype=dtype)
                    properties = bed.properties
                    val[-1, 0] = float("NAN")
                    output = str(
                        tmp_path
                        / "toydata.{0}{1}.cpp".format(
                            order, "32" if dtype == np.float32 else "64"
                        )
                    )
                    to_bed(output, val, properties=properties, count_A1=count_A1)
                    val2 = open_bed(output, count_A1=count_A1).read(dtype=dtype)
                    assert np.allclose(val, val2, equal_nan=True)
                    val_sparse = open_bed(output, count_A1=count_A1).read_sparse(
                        dtype=dtype, format=format
                    )
                    assert np.allclose(val, val_sparse.toarray(), equal_nan=True)


def test_write1_x_x_cpp_stream(tmp_path, shared_datadir):
    for count_A1 in [False, True]:
        with open_bed(shared_datadir / "some_missing.bed", count_A1=count_A1) as bed:
            for order, format in [("F", "csc"), ("C", "csr")]:
                for dtype in [np.float32, np.float64]:
                    val = bed.read(order=order, dtype=dtype)
                    properties = bed.properties
                    val[-1, 0] = float("NAN")
                    output = str(
                        tmp_path
                        / "toydata.{0}{1}.cpp".format(
                            order, "32" if dtype == np.float32 else "64"
                        )
                    )
                    with create_bed(
                        output,
                        iid_count=val.shape[0],
                        sid_count=val.shape[1],
                        properties=properties,
                        count_A1=count_A1,
                    ) as bed_writer:
                        for sid_index in range(val.shape[1]):
                            bed_writer.write(val[:, sid_index])
                    val2 = open_bed(output, count_A1=count_A1).read(dtype=dtype)
                    assert np.allclose(val, val2, equal_nan=True)
                    val_sparse = open_bed(output, count_A1=count_A1).read_sparse(
                        dtype=dtype, format=format
                    )
                    assert np.allclose(val, val_sparse.toarray(), equal_nan=True)


def test_respect_read_inputs(shared_datadir):
    import scipy.sparse as sparse

    ref_val_float = reference_val(shared_datadir)
    ref_val_float2 = ref_val_float.copy()
    ref_val_float2[ref_val_float != ref_val_float] = -127
    ref_val_int8 = ref_val_float2.astype("int8")

    with open_bed(shared_datadir / "some_missing.bed") as bed:
        for order, format in [("F", "csc"), ("C", "csr")]:
            for dtype in [np.int8, np.float32, np.float64]:
                for force_python_only in [True, False]:
                    val = bed.read(
                        order=order, dtype=dtype, force_python_only=force_python_only
                    )
                    has_right_order = (order == "C" and val.flags["C_CONTIGUOUS"]) or (
                        order == "F" and val.flags["F_CONTIGUOUS"]
                    )
                    assert val.dtype == dtype and has_right_order
                    ref_val = ref_val_int8 if dtype == np.int8 else ref_val_float
                    assert np.allclose(ref_val, val, equal_nan=True)

                val_sparse = bed.read_sparse(dtype=dtype, format=format)
                has_right_format = (
                    format == "csc" and isinstance(val_sparse, sparse.csc_matrix)
                ) or (format == "csr" and isinstance(val_sparse, sparse.csr_matrix))
                assert val_sparse.dtype == dtype and has_right_format
                assert np.allclose(ref_val, val_sparse.toarray(), equal_nan=True)


def test_threads(shared_datadir):
    ref_val_float = reference_val(shared_datadir)
    ref_val_float2 = ref_val_float.copy()
    ref_val_float2[ref_val_float != ref_val_float] = -127
    ref_val_int8 = ref_val_float2.astype("int8")

    for num_threads in [1, 4]:
        with open_bed(
            shared_datadir / "some_missing.bed", num_threads=num_threads
        ) as bed:
            val = bed.read(dtype="int8")
            assert np.allclose(ref_val_int8, val, equal_nan=True)
            val_sparse = bed.read_sparse(dtype="int8")
            assert np.allclose(ref_val_int8, val_sparse.toarray(), equal_nan=True)


def test_write12(tmp_path):
    # ===================================
    #    Starting main function
    # ===================================
    logging.info("starting 'test_writes'")
    np.random.seed(0)
    output_template = str(tmp_path / "writes.{0}.bed")
    i = 0
    for row_count in [0, 5, 2, 1]:
        for col_count in [0, 4, 2, 1]:
            val = np.random.randint(0, 4, size=(row_count, col_count)) * 1.0
            val[val == 3] = np.nan
            row0 = ["0", "1", "2", "3", "4"][:row_count]
            row1 = ["0", "1", "2", "3", "4"][:row_count]
            col = ["s0", "s1", "s2", "s3", "s4"][:col_count]
            for is_none in [True, False]:
                properties = {"fid": row0, "iid": row1, "sid": col}
                if is_none:
                    col_prop012 = [x for x in range(5)][:col_count]
                    properties["chromosome"] = col_prop012
                    properties["bp_position"] = col_prop012
                    properties["cm_position"] = col_prop012
                else:
                    col_prop012 = None

                filename = output_template.format(i)
                logging.info(filename)
                i += 1
                to_bed(filename, val, properties=properties)
                for subsetter in [None, np.s_[::2, ::3]]:
                    with open_bed(filename) as bed:
                        val2 = bed.read(index=subsetter, order="C", dtype="float32")
                        if subsetter is None:
                            expected = val
                        else:
                            expected = val[subsetter[0], :][:, subsetter[1]]
                        assert np.allclose(val2, expected, equal_nan=True)
                        assert np.array_equal(bed.fid, np.array(row0, dtype="str"))
                        assert np.array_equal(bed.iid, np.array(row1, dtype="str"))
                        assert np.array_equal(bed.sid, np.array(col, dtype="str"))
                        if col_prop012 is not None:
                            assert np.array_equal(
                                bed.chromosome, np.array(col_prop012, dtype="str")
                            )
                            assert np.array_equal(
                                bed.bp_position, np.array(col_prop012)
                            )
                            assert np.array_equal(
                                bed.cm_position, np.array(col_prop012)
                            )
                try:
                    os.remove(filename)
                except Exception:
                    pass
    logging.info("done with 'test_writes'")


def test_writes_small(tmp_path):
    output_file = tmp_path / "small.bed"

    val = [[1.0, 0, np.nan, 0], [2, 0, np.nan, 2], [0, 1, 2, 0]]

    properties = {
        "fid": ["fid1", "fid1", "fid2"],
        "iid": ["iid1", "iid2", "iid3"],
        "father": ["iid23", "iid23", "iid22"],
        "mother": ["iid34", "iid34", "iid33"],
        "sex": [1, 2, 0],
        "pheno": ["red", "red", "blue"],
        "chromosome": ["1", "1", "5", "Y"],
        "sid": ["sid1", "sid2", "sid3", "sid4"],
        "cm_position": [100.4, 2000.5, 4000.7, 7000.9],
        "bp_position": [1, 100, 1000, 1004],
        "allele_1": ["A", "T", "A", "T"],
        "allele_2": ["A", "C", "C", "G"],
    }

    to_bed(output_file, val, properties=properties)

    with open_bed(output_file) as bed:
        assert np.allclose(bed.read(), val, equal_nan=True)
        for key, value in bed.properties.items():
            assert np.array_equal(value, properties[key]) or np.allclose(
                value, properties[key]
            )


def test_writes_small_stream(tmp_path):
    output_file = tmp_path / "small.bed"

    val = np.array([[1.0, 0, np.nan, 0], [2, 0, np.nan, 2], [0, 1, 2, 0]])

    properties = {
        "fid": ["fid1", "fid1", "fid2"],
        "iid": ["iid1", "iid2", "iid3"],
        "father": ["iid23", "iid23", "iid22"],
        "mother": ["iid34", "iid34", "iid33"],
        "sex": [1, 2, 0],
        "pheno": ["red", "red", "blue"],
        "chromosome": ["1", "1", "5", "Y"],
        "sid": ["sid1", "sid2", "sid3", "sid4"],
        "cm_position": [100.4, 2000.5, 4000.7, 7000.9],
        "bp_position": [1, 100, 1000, 1004],
        "allele_1": ["A", "T", "A", "T"],
        "allele_2": ["A", "C", "C", "G"],
    }

    with create_bed(
        output_file,
        iid_count=val.shape[0],
        sid_count=val.shape[1],
        properties=properties,
    ) as bed_writer:
        for column_data in val.T:
            bed_writer.write(column_data)

    with open_bed(output_file) as bed:
        assert np.allclose(bed.read(), val, equal_nan=True)
        for key, value in bed.properties.items():
            assert np.array_equal(value, properties[key]) or np.allclose(
                value, properties[key]
            )


def test_index(shared_datadir):
    ref_val_float = reference_val(shared_datadir)

    with open_bed(shared_datadir / "some_missing.bed") as bed:
        val = bed.read()
        assert np.allclose(ref_val_float, val, equal_nan=True)
        val_sparse = bed.read_sparse()
        assert np.allclose(ref_val_float, val_sparse.toarray(), equal_nan=True)

        val = bed.read(2)
        assert np.allclose(ref_val_float[:, [2]], val, equal_nan=True)
        val_sparse = bed.read_sparse(2)
        assert np.allclose(ref_val_float[:, [2]], val_sparse.toarray(), equal_nan=True)

        val = bed.read((2))
        assert np.allclose(ref_val_float[:, [2]], val, equal_nan=True)
        val_sparse = bed.read_sparse((2))
        assert np.allclose(ref_val_float[:, [2]], val_sparse.toarray(), equal_nan=True)

        val = bed.read((None, 2))
        assert np.allclose(ref_val_float[:, [2]], val, equal_nan=True)
        val_sparse = bed.read_sparse((None, 2))
        assert np.allclose(ref_val_float[:, [2]], val_sparse.toarray(), equal_nan=True)

        val = bed.read((1, 2))
        assert np.allclose(ref_val_float[[1], [2]], val, equal_nan=True)
        val_sparse = bed.read_sparse((1, 2))
        assert np.allclose(
            ref_val_float[[1], [2]], val_sparse.toarray(), equal_nan=True
        )

        val = bed.read([2, -2])
        assert np.allclose(ref_val_float[:, [2, -2]], val, equal_nan=True)
        val_sparse = bed.read_sparse([2, -2])
        assert np.allclose(
            ref_val_float[:, [2, -2]], val_sparse.toarray(), equal_nan=True
        )

        val = bed.read(([1, -1], [2, -2]))
        assert np.allclose(ref_val_float[[1, -1], :][:, [2, -2]], val, equal_nan=True)
        val_sparse = bed.read_sparse(([1, -1], [2, -2]))
        assert np.allclose(
            ref_val_float[[1, -1], :][:, [2, -2]], val_sparse.toarray(), equal_nan=True
        )

        iid_bool = ([False, False, True] * bed.iid_count)[: bed.iid_count]
        sid_bool = ([True, False, True] * bed.sid_count)[: bed.sid_count]
        val = bed.read(sid_bool)
        assert np.allclose(ref_val_float[:, sid_bool], val, equal_nan=True)
        val_sparse = bed.read_sparse(sid_bool)
        assert np.allclose(
            ref_val_float[:, sid_bool], val_sparse.toarray(), equal_nan=True
        )

        val = bed.read((iid_bool, sid_bool))
        assert np.allclose(ref_val_float[iid_bool, :][:, sid_bool], val, equal_nan=True)
        val_sparse = bed.read_sparse((iid_bool, sid_bool))

        val = bed.read((1, sid_bool))
        assert np.allclose(ref_val_float[[1], :][:, sid_bool], val, equal_nan=True)
        val_sparse = bed.read_sparse((1, sid_bool))
        assert np.allclose(
            ref_val_float[[1], :][:, sid_bool], val_sparse.toarray(), equal_nan=True
        )

        slicer = np.s_[::2, ::3]
        val = bed.read(slicer[1])
        assert np.allclose(ref_val_float[:, slicer[1]], val, equal_nan=True)
        val_sparse = bed.read_sparse(slicer[1])
        assert np.allclose(
            ref_val_float[:, slicer[1]], val_sparse.toarray(), equal_nan=True
        )

        val = bed.read(slicer)
        assert np.allclose(ref_val_float[slicer], val, equal_nan=True)
        val_sparse = bed.read_sparse(slicer)
        assert np.allclose(ref_val_float[slicer], val_sparse.toarray(), equal_nan=True)

        val = bed.read((1, slicer[1]))
        assert np.allclose(ref_val_float[[1], slicer[1]], val, equal_nan=True)
        val_sparse = bed.read_sparse((1, slicer[1]))
        assert np.allclose(
            ref_val_float[[1], slicer[1]], val_sparse.toarray(), equal_nan=True
        )


def test_shape(shared_datadir):
    with open_bed(shared_datadir / "plink_sim_10s_100v_10pmiss.bed") as bed:
        assert bed.shape == (10, 100)


def test_zero_files(tmp_path):
    for force_python_only, format in [(False, "csc"), (True, "csr")]:
        for iid_count in [3, 0]:
            for sid_count in [0, 5]:
                for dtype in [np.int8, np.float32, np.float64]:
                    val = np.zeros((iid_count, sid_count), dtype=dtype)
                    if iid_count * sid_count > 0:
                        val[0, 0] = 2
                        val[0, 1] = -127 if np.dtype(dtype) == np.int8 else np.nan
                    filename = str(tmp_path / "zero_files.bed")

                    # Write
                    to_bed(filename, val, force_python_only=force_python_only)

                    # Read
                    with open_bed(filename) as bed2:
                        val2 = bed2.read(dtype=dtype)
                        assert np.allclose(val, val2, equal_nan=True)
                        val_sparse = bed2.read_sparse(dtype=dtype, format=format)
                        assert np.allclose(val, val_sparse.toarray(), equal_nan=True)
                        properties2 = bed2.properties
                        for prop in properties2.values():
                            assert len(prop) in {iid_count, sid_count}

                    # Change properties and write again
                    if iid_count > 0:
                        properties2["iid"][0] = "iidx"
                    if sid_count > 0:
                        properties2["sid"][0] = "sidx"
                    to_bed(
                        filename,
                        val2,
                        properties=properties2,
                        force_python_only=force_python_only,
                    )

                    # Read again
                    with open_bed(filename) as bed3:
                        val3 = bed3.read(dtype=dtype)
                        assert np.allclose(val, val3, equal_nan=True)
                        val_sparse = bed3.read_sparse(dtype=dtype, format=format)
                        assert np.allclose(val, val_sparse.toarray(), equal_nan=True)
                        properties3 = bed3.properties
                        for key2, value_list2 in properties2.items():
                            value_list3 = properties3[key2]
                            assert np.array_equal(value_list2, value_list3)


def test_iid_sid_count(shared_datadir):
    iid_count_ref, sid_count_ref = open_bed(
        shared_datadir / "plink_sim_10s_100v_10pmiss.bed"
    ).shape
    assert (iid_count_ref, sid_count_ref) == open_bed(
        shared_datadir / "plink_sim_10s_100v_10pmiss.bed", iid_count=iid_count_ref
    ).shape
    assert (iid_count_ref, sid_count_ref) == open_bed(
        shared_datadir / "plink_sim_10s_100v_10pmiss.bed", sid_count=sid_count_ref
    ).shape
    assert (iid_count_ref, sid_count_ref) == open_bed(
        shared_datadir / "plink_sim_10s_100v_10pmiss.bed",
        iid_count=iid_count_ref,
        sid_count=sid_count_ref,
    ).shape


def test_sample_file():
    from bed_reader import open_bed, sample_file

    file_name = sample_file("small.bed")
    with open_bed(file_name) as bed:
        print(bed.iid)
        print(bed.sid)
        print(bed.read())
        print(bed.read_sparse())


def test_coverage2(shared_datadir, tmp_path):
    with open_bed(
        shared_datadir / "plink_sim_10s_100v_10pmiss.bed", properties={"iid": None}
    ) as bed:
        assert bed.iid is None
    with pytest.raises(ValueError):
        open_bed(
            shared_datadir / "plink_sim_10s_100v_10pmiss.bed",
            properties={"iid": [1, 2, 3], "mother": [1, 2]},
        )
    val = np.zeros((3, 5))[::2]
    assert not val.flags["C_CONTIGUOUS"] and not val.flags["F_CONTIGUOUS"]
    with pytest.raises(ValueError):
        to_bed(tmp_path / "ignore", val)
    val = np.zeros((3, 5), dtype=np.str_)
    with pytest.raises(ValueError):
        to_bed(tmp_path / "ignore", val)


def test_coverage3(shared_datadir, tmp_path):
    with open_bed(
        shared_datadir / "small.bed", properties={"sex": [1.0, np.nan, 1.0, 2.0]}
    ) as bed:
        assert np.array_equal(bed.sex, np.array([1, 0, 1, 2]))

    with open_bed(
        shared_datadir / "small.bed",
        properties={"cm_position": [1000.0, np.nan, 2000.0, 3000.0]},
    ) as bed:
        assert np.array_equal(bed.cm_position, np.array([1000, 0, 2000, 3000]))

    list = [1.0, 0, np.nan, 0]
    output_file = tmp_path / "1d.bed"
    with pytest.raises(ValueError):
        to_bed(output_file, list)
    to_bed(output_file, np.array([list], dtype=np.float16))
    with open_bed(output_file) as bed:
        assert np.allclose(bed.read(), [list], equal_nan=True)
        assert np.allclose(bed.read_sparse().toarray(), [list], equal_nan=True)


def test_nones(shared_datadir, tmp_path):
    properties = {
        "father": None,
        "mother": None,
        "sex": None,
        "pheno": None,
        "allele_1": None,
        "allele_2": None,
    }

    with open_bed(shared_datadir / "small.bed", properties=properties) as bed:
        assert np.array_equal(bed.iid, ["iid1", "iid2", "iid3"])
        assert bed.father is None

    val = [[1.0, 0, np.nan, 0], [2, 0, np.nan, 2], [0, 1, 2, 0]]
    out_file = tmp_path / "testnones.bed"
    to_bed(out_file, val, properties=properties)


def test_fam_bim_filepath(shared_datadir, tmp_path):
    with open_bed(shared_datadir / "small.bed") as bed:
        val = bed.read()
        properties = bed.properties
    output_file = tmp_path / "small.deb"
    fam_file = tmp_path / "small.maf"
    bim_file = tmp_path / "small.mib"
    to_bed(
        output_file,
        val,
        properties=properties,
        fam_filepath=fam_file,
        bim_filepath=bim_file,
    )
    assert output_file.exists() and fam_file.exists() and bim_file.exists()

    with open_bed(output_file, fam_filepath=fam_file, bim_filepath=bim_file) as deb:
        val2 = deb.read()
        assert np.allclose(val, val2, equal_nan=True)
        val_sparse = deb.read_sparse()
        assert np.allclose(val, val_sparse.toarray(), equal_nan=True)
        properties2 = deb.properties
        for key in properties:
            np.array_equal(properties[key], properties2[key])


def test_write_nan_properties(shared_datadir, tmp_path):
    with open_bed(shared_datadir / "small.bed") as bed:
        val = bed.read()
        properties = bed.properties
        chrom = bed.chromosome.copy()
        chrom[bed.chromosome == "Y"] = 0
        chrom = np.array(chrom, dtype="float")
        chrom2 = chrom.copy()
        chrom2[chrom2 == 0] = np.nan
        cm_p = bed.cm_position.copy()
        cm_p[cm_p < 3000] = 0
        cm_p2 = cm_p.copy()
        cm_p2[cm_p == 0] = np.nan
        properties["chromosome"] = chrom2
        properties["cm_position"] = cm_p2

    output_file = tmp_path / "nan.bed"
    to_bed(output_file, val, properties=properties)

    with open_bed(output_file) as bed2:
        assert np.array_equal(bed2.chromosome, ["1.0", "1.0", "5.0", "0"])
        assert np.array_equal(bed2.cm_position, cm_p)

    with open_bed(
        shared_datadir / "small.bed",
        properties={"chromosome": chrom2, "cm_position": cm_p2},
    ) as bed3:
        assert np.array_equal(bed3.chromosome, ["1.0", "1.0", "5.0", "0"])
        assert np.array_equal(bed3.cm_position, cm_p)


def test_write_nan_properties_stream(shared_datadir, tmp_path):
    with open_bed(shared_datadir / "small.bed") as bed:
        val = bed.read()
        properties = bed.properties
        chrom = bed.chromosome.copy()
        chrom[bed.chromosome == "Y"] = 0
        chrom = np.array(chrom, dtype="float")
        chrom2 = chrom.copy()
        chrom2[chrom2 == 0] = np.nan
        cm_p = bed.cm_position.copy()
        cm_p[cm_p < 3000] = 0
        cm_p2 = cm_p.copy()
        cm_p2[cm_p == 0] = np.nan
        properties["chromosome"] = chrom2
        properties["cm_position"] = cm_p2

    output_file = tmp_path / "nan.bed"

    with create_bed(
        output_file,
        iid_count=val.shape[0],
        sid_count=val.shape[1],
        properties=properties,
    ) as bed_writer:
        for column_data in val.T:
            bed_writer.write(column_data)

    with open_bed(output_file) as bed2:
        assert np.array_equal(bed2.chromosome, ["1.0", "1.0", "5.0", "0"])
        assert np.array_equal(bed2.cm_position, cm_p)

    with open_bed(
        shared_datadir / "small.bed",
        properties={"chromosome": chrom2, "cm_position": cm_p2},
    ) as bed3:
        assert np.array_equal(bed3.chromosome, ["1.0", "1.0", "5.0", "0"])
        assert np.array_equal(bed3.cm_position, cm_p)


def test_env(shared_datadir):
    if platform.system() == "Darwin":
        return

    key = "MKL_NUM_THREADS"
    original_val = os.environ.get(key)
    try:
        os.environ[key] = "1"
        with open_bed(shared_datadir / "some_missing.bed") as bed:
            _ = bed.read(np.s_[:100, :100])
            _ = bed.read_sparse(np.s_[:100, :100])
        os.environ[key] = "10"
        with open_bed(shared_datadir / "some_missing.bed") as bed:
            _ = bed.read(np.s_[:100, :100])
            _ = bed.read_sparse(np.s_[:100, :100])
        os.environ[key] = "BADVALUE"
        with pytest.raises(ValueError):
            with open_bed(shared_datadir / "some_missing.bed") as bed:
                _ = bed.read(np.s_[:100, :100])
        with pytest.raises(ValueError):
            with open_bed(shared_datadir / "some_missing.bed") as bed:
                _ = bed.read_sparse(np.s_[:100, :100])
    finally:
        if original_val is None:
            if key in os.environ:
                del os.environ[key]
        else:
            os.environ[key] = original_val


def test_noncontig_indexes(shared_datadir):
    with open_bed(shared_datadir / "some_missing.bed") as bed:
        whole_iid_index = np.arange(bed.iid_count)
        assert whole_iid_index.flags["C_CONTIGUOUS"]
        every_other = whole_iid_index[::2]
        assert not every_other.flags["C_CONTIGUOUS"]
        val = bed.read((every_other, -2))

        whole_iid_index = np.arange(val.shape[0])
        assert whole_iid_index.flags["C_CONTIGUOUS"]
        every_other = whole_iid_index[::2]
        assert not every_other.flags["C_CONTIGUOUS"]
        val_out = np.zeros((len(every_other), 0))
        with pytest.raises(ValueError):
            subset_f64_f64(
                val.reshape(-1, bed.sid_count, 1), every_other, [], val_out, 1
            )


def test_bed_reading_example():
    import numpy as np

    from bed_reader import open_bed, sample_file

    file_name = sample_file("small.bed")
    with open_bed(file_name, count_A1=False) as bed:
        val = bed.read(index=np.s_[:, :3], dtype="int8", order="C", num_threads=1)
        print(val.shape)


def test_sparse():
    import numpy as np

    from bed_reader import open_bed, sample_file

    file_name = sample_file("small.bed")
    with open_bed(file_name, count_A1=False) as bed:
        val_sparse = bed.read_sparse(index=np.s_[:, :3], dtype="int8")
        print(val_sparse.shape)


def test_convert_to_dtype():
    from bed_reader._open_bed import _convert_to_dtype

    input = [
        [["a", "b", "c"], ["a", "b", "c"], None, None],
        [["1.0", "2.0", "3.0"], ["1.0", "2.0", "3.0"], [1, 2, 3], [1.0, 2.0, 3.0]],
        [["1.0", "2.0", "3.5"], ["1.0", "2.0", "3.5"], None, [1.0, 2.0, 3.5]],
        [["1", "2", "3"], ["1", "2", "3"], [1, 2, 3], [1.0, 2.0, 3.0]],
        [["1", "A", "3"], ["1", "A", "3"], None, None],
    ]
    # convert all to np.array
    input = [
        [np.array(inner) if inner is not None else None for inner in outer]
        for outer in input
    ]

    for ori, exp_str, exp_int, exp_float in input:
        for dtype, exp in (
            [np.str_, exp_str],
            [np.int32, exp_int],
            [
                np.float32,
                exp_float,
            ],
        ):
            try:
                actual = _convert_to_dtype(ori, dtype)
                assert np.array_equal(actual, exp)
            except ValueError as e:
                print(e)
                assert exp is None


def test_stream_i8(shared_datadir, tmp_path):
    # test major = "SNP" vs "individual"
    with open_bed(shared_datadir / "some_missing.bed") as bed2:
        val = bed2.read(dtype=np.int8)
        properties = bed2.properties
    for count_A1 in [True, False]:
        for fam_location, bim_location in [
            (None, None),
            (tmp_path / "some_missing.fam2", tmp_path / "some_missing.bim2"),
        ]:
            for force_python_only in [True, False]:
                with create_bed(
                    tmp_path / "some_missing.bed",
                    iid_count=val.shape[0],
                    sid_count=val.shape[1],
                    properties=properties,
                    fam_location=fam_location,
                    bim_location=bim_location,
                    count_A1=count_A1,
                    force_python_only=force_python_only,
                ) as bed_writer:
                    for column_data in val.T:
                        bed_writer.write(column_data)
                with open_bed(
                    tmp_path / "some_missing.bed",
                    count_A1=count_A1,
                    fam_filepath=fam_location,
                    bim_filepath=bim_location,
                ) as bed2:
                    val2 = bed2.read(dtype=np.int8)
                    assert np.allclose(val, val2, equal_nan=True)
                    properties2 = bed2.properties
                    for key, value in properties.items():
                        assert np.array_equal(value, properties2[key]) or np.allclose(
                            value, properties2[key]
                        )


def test_stream(shared_datadir, tmp_path):
    # test major = "SNP" vs "individual"
    with open_bed(shared_datadir / "some_missing.bed") as bed2:
        val = bed2.read()
        properties = bed2.properties
    for count_A1 in [False, True]:
        for fam_location, bim_location in [
            (None, None),
            (tmp_path / "some_missing.fam2", tmp_path / "some_missing.bim2"),
        ]:
            with create_bed(
                tmp_path / "some_missing.bed",
                iid_count=val.shape[0],
                sid_count=val.shape[1],
                properties=properties,
                fam_location=fam_location,
                bim_location=bim_location,
                count_A1=count_A1,
            ) as bed_writer:
                for column_data in val.T:
                    bed_writer.write(column_data)
            with open_bed(
                tmp_path / "some_missing.bed",
                count_A1=count_A1,
                fam_filepath=fam_location,
                bim_filepath=bim_location,
            ) as bed2:
                val2 = bed2.read()
                assert np.allclose(val, val2, equal_nan=True)
                properties2 = bed2.properties
                for key, value in properties.items():
                    assert np.array_equal(value, properties2[key]) or np.allclose(
                        value, properties2[key]
                    )


def test_stream_individual_major(shared_datadir, tmp_path):
    with open_bed(shared_datadir / "small.bed") as bed:
        val = bed.read()
        properties = bed.properties
    with create_bed(
        tmp_path / "small.bed",
        iid_count=val.shape[0],
        sid_count=val.shape[1],
        properties=properties,
        major="individual",
    ) as bed_writer:
        for row_data in val:
            bed_writer.write(row_data)
    with open_bed(
        tmp_path / "small.bed",
    ) as bed2:
        val2 = bed2.read()
        assert np.allclose(val, val2, equal_nan=True)
        properties2 = bed2.properties
        for key, value in properties.items():
            assert np.array_equal(value, properties2[key]) or np.allclose(
                value, properties2[key]
            )


def test_stream_without_with(shared_datadir, tmp_path):
    with open_bed(shared_datadir / "some_missing.bed") as bed2:
        val = bed2.read()
        properties = bed2.properties
    bed_writer = create_bed(
        tmp_path / "some_missing.bed",
        iid_count=val.shape[0],
        sid_count=val.shape[1],
        properties=properties,
    )
    for column_data in val.T:
        bed_writer.write(column_data)
    bed_writer.close()
    with open_bed(
        tmp_path / "some_missing.bed",
    ) as bed2:
        val2 = bed2.read()
        assert np.allclose(val, val2, equal_nan=True)
        properties2 = bed2.properties
        for key, value in properties.items():
            assert np.array_equal(value, properties2[key]) or np.allclose(
                value, properties2[key]
            )


def test_stream_expected_errors(shared_datadir, tmp_path):
    # test major = "SNP" vs "individual"
    with open_bed(shared_datadir / "some_missing.bed") as bed2:
        val = bed2.read()
        properties = bed2.properties
    # expect value error when major is set wrong
    with pytest.raises(ValueError):
        with create_bed(
            tmp_path / "some_missing.bed",
            major="illegal value",
            iid_count=val.shape[0],
            sid_count=val.shape[1],
            properties=properties,
        ) as bed_writer:
            for column_data in val.T:
                bed_writer.write(column_data)
    # properties disagree with given counts
    with pytest.raises(ValueError):
        with create_bed(
            tmp_path / "some_missing.bed",
            iid_count=val.shape[0] + 1,
            sid_count=val.shape[1],
            properties=properties,
        ) as bed_writer:
            for column_data in val.T:
                bed_writer.write(column_data)
    with pytest.raises(ValueError):
        with create_bed(
            tmp_path / "some_missing.bed",
            iid_count=val.shape[0],
            sid_count=val.shape[1] - 1,
            properties=properties,
        ) as bed_writer:
            for column_data in val.T:
                bed_writer.write(column_data)

    # write to few lines
    with pytest.raises(ValueError):
        with create_bed(
            tmp_path / "some_missing.bed",
            iid_count=val.shape[0],
            sid_count=val.shape[1],
            properties=properties,
        ) as bed_writer:
            for column_data in val.T:
                bed_writer.write(column_data)
                break

    # write to many lines
    with pytest.raises(ValueError):
        with create_bed(
            tmp_path / "some_missing.bed",
            iid_count=val.shape[0],
            sid_count=val.shape[1],
            properties=properties,
        ) as bed_writer:
            for column_data in val.T:
                bed_writer.write(column_data)
                bed_writer.write(column_data)

    # try to write after close
    with pytest.raises(RuntimeError):
        bed_writer = create_bed(
            tmp_path / "some_missing.bed",
            iid_count=val.shape[0],
            sid_count=val.shape[1],
            properties=properties,
        )
        for column_data in val.T:
            bed_writer.write(column_data)
        bed_writer.close()
        bed_writer.write(val[:, 0])

    # vector is the wrong length
    with pytest.raises(ValueError):
        with create_bed(
            tmp_path / "some_missing.bed",
            iid_count=val.shape[0],
            sid_count=val.shape[1],
            properties=properties,
        ) as bed_writer:
            for column_data in val.T:
                bed_writer.write(column_data[:-1])
    # vector is 2-d or a scalar
    with pytest.raises(ValueError):
        with create_bed(
            tmp_path / "some_missing.bed",
            iid_count=val.shape[0],
            sid_count=val.shape[1],
            properties=properties,
        ) as bed_writer:
            for column_data in val.T:
                bed_writer.write(val)
    with pytest.raises(ValueError):
        with create_bed(
            tmp_path / "some_missing.bed",
            iid_count=val.shape[0],
            sid_count=val.shape[1],
            properties=properties,
        ) as bed_writer:
            for column_data in val.T:
                bed_writer.write(3)


# Show that can read files of both modes by reading smallmode0.bed and small.bed
def test_mode(shared_datadir):
    for force_python_only in [False, True]:
        with open_bed(shared_datadir / "small.bed") as bed:
            val = bed.read(force_python_only=force_python_only)
        with open_bed(shared_datadir / "smallmode0.bed") as bed2:
            val2 = bed2.read(force_python_only=force_python_only)
        assert np.allclose(val, val2.T, equal_nan=True)


# Show that can write via streaming in both modes
def test_stream_modes(shared_datadir, tmp_path):
    for force_python_only in [False, True]:
        with open_bed(shared_datadir / "small.bed") as bed:
            val = bed.read(force_python_only=force_python_only)

        # write in mode 1
        with create_bed(
            tmp_path / "mode1.bed",
            major="SNP",
            iid_count=bed.iid_count,
            sid_count=bed.sid_count,
            properties=bed.properties,
            force_python_only=force_python_only,
        ) as bed_writer:
            for column_data in val.T:
                bed_writer.write(column_data)

        with open_bed(tmp_path / "mode1.bed") as bed2:
            val2 = bed2.read(force_python_only=force_python_only)
        assert np.allclose(val, val2, equal_nan=True)
        properties2 = bed2.properties
        for key, value in bed.properties.items():
            assert np.array_equal(value, properties2[key]) or np.allclose(
                value, properties2[key]
            )

        # write in mode 0
        with create_bed(
            tmp_path / "mode0.bed",
            major="individual",
            iid_count=bed.iid_count,
            sid_count=bed.sid_count,
            properties=bed.properties,
            force_python_only=force_python_only,
        ) as bed_writer:
            for snp_data in val:
                bed_writer.write(snp_data)
        with open_bed(tmp_path / "mode0.bed") as bed3:
            val3 = bed3.read(force_python_only=force_python_only)
        assert np.allclose(val, val3, equal_nan=True)
        properties3 = bed3.properties
        for key, value in bed.properties.items():
            assert np.array_equal(value, properties3[key]) or np.allclose(
                value, properties3[key]
            )


# show can can write different minor-lengths in both modes (including 0)
def test_stream_modes2(shared_datadir, tmp_path):
    for iid_count in range(8, -1, -1):
        for sid_count in range(8, -1, -1):
            for force_python_only in [False, True]:
                for order in ["f", "c"]:
                    rng = np.random.default_rng(0)
                    # fill values with 0,1, 2
                    val = rng.integers(3, size=(iid_count, sid_count), dtype=np.int8)
                    val = np.ascontiguousarray(np.array(val, order=order))
                    # make random value missing
                    if iid_count > 0 and sid_count > 0:
                        val[rng.integers(iid_count), rng.integers(sid_count)] = -127
                    # write in mode 1
                    with create_bed(
                        tmp_path / "mode1.bed",
                        major="SNP",
                        iid_count=iid_count,
                        sid_count=sid_count,
                        force_python_only=force_python_only,
                    ) as bed_writer:
                        for column_data in val.T:
                            bed_writer.write(column_data)
                    with open_bed(tmp_path / "mode1.bed") as bed2:
                        val2 = bed2.read(
                            force_python_only=force_python_only, dtype=np.int8
                        )
                    assert np.allclose(val, val2, equal_nan=True)

                    # write in mode 0
                    with create_bed(
                        tmp_path / "mode0.bed",
                        major="individual",
                        iid_count=iid_count,
                        sid_count=sid_count,
                        force_python_only=force_python_only,
                    ) as bed_writer:
                        for snp_data in val:
                            bed_writer.write(snp_data)
                    with open_bed(tmp_path / "mode0.bed") as bed3:
                        val3 = bed3.read(
                            force_python_only=force_python_only, dtype=np.int8
                        )
                    assert np.allclose(val, val3, equal_nan=True)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    shared_datadir = Path(r"D:\OneDrive\programs\bed-reader\bed_reader\tests\data")
    tmp_path = Path(r"m:/deldir/tests")
    test_bed_reading_example()
    # test_zero_files(tmp_path)
    # test_index(shared_datadir)
    # test_c_reader_bed(shared_datadir)
    # test_read1(shared_datadir)
    pytest.main([__file__])
