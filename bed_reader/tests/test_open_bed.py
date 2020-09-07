import logging
import os
import platform
from pathlib import Path

import numpy as np
import pytest

from bed_reader import open_bed, to_bed


def test_read1(shared_datadir):
    file = shared_datadir / "plink_sim_10s_100v_10pmiss.bed"
    with open_bed(file) as bed:
        assert bed.iid_count == 10
        assert bed.fid[-1] == "0"
        assert bed.iid[-1] == "9"
        assert bed.shape == (10, 100)
        val = bed.read(dtype="int8")
        assert (
            val.mean() == -13.142
        )  # really shouldn't do mean on data where -127 represents missing
        assert bed.chromosome[-1] == "1"
        assert bed.bp_position[-1] == 100


def test_write1(tmp_path, shared_datadir):
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
            assert np.array_equal(bed.fid, properties0["fid"])
            assert np.array_equal(bed.iid, properties0["iid"])
            assert np.array_equal(bed.sid, properties0["sid"])
            assert np.array_equal(bed.chromosome, properties0["chromosome"])
            assert np.allclose(bed.cm_position, properties0["cm_position"])
            assert np.allclose(bed.bp_position, properties0["bp_position"])

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
    # lock in the expected results: np.savez(shared_datadir / "some_missing.properties.npz",fid=fid,iid=iid,father=father,mother=mother,sex=sex,pheno=pheno,chromosome=chromosome,sid=sid,cm_position=cm_position,bp_position=bp_position,allele_1=allele_1,allele_2=allele_2)
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


def setting_generator(seq_dict, seed=9392):
    import itertools
    from numpy.random import RandomState

    longest = max((len(value_list) for value_list in seq_dict.values()))

    for test_index in range(longest):
        setting = {}
        for offset, (key, value_list) in enumerate(seq_dict.items()):
            val = value_list[(test_index + offset) % len(value_list)]
            if "leave_out" is not val:
                setting[key] = val
        yield setting

    all_combo = list(itertools.product(*seq_dict.values()))

    random_state = RandomState(seed)
    random_state.shuffle(all_combo)
    for combo in all_combo:
        setting = {
            key: value_list
            for key, value_list in itertools.zip_longest(seq_dict, combo)
            if value_list is not "leave_out"
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
        "chromosome": ["leave_out", None, chromosome_list, np.array(chromosome_list),],
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
            assert val.shape == (len(iid_list), len(sid_list),)
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
    for force_python_only in [False, True]:
        bed = open_bed(shared_datadir / "some_missing.bed", count_A1=False)

        val = bed.read(order="F", force_python_only=force_python_only)
        assert val.dtype == np.float32
        ref_val = reference_val(shared_datadir)
        ref_val = ref_val * -1 + 2
        assert np.allclose(ref_val, val, rtol=1e-05, atol=1e-05, equal_nan=True)

        val = bed.read(order="F", dtype="int8", force_python_only=False)
        assert val.dtype == np.int8
        ref_val[ref_val != ref_val] = -127
        ref_val = ref_val.astype("int8")
        ref_val = ref_val.astype("int8")
        assert np.all(ref_val == val)

        del bed

        with open_bed(shared_datadir / "some_missing.bed") as bed:
            val = bed.read(
                order="F", dtype="float64", force_python_only=force_python_only
            )
            ref_val = reference_val(shared_datadir)
            assert np.allclose(ref_val, val, rtol=1e-05, atol=1e-05, equal_nan=True)


def reference_val(shared_datadir):
    val = np.load(shared_datadir / "some_missing.val.npy")
    return val


def test_bed_int8(tmp_path, shared_datadir):
    with open_bed(shared_datadir / "some_missing.bed") as bed:
        for force_python_only in [False, True]:
            for order in ["F", "C"]:
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
                    to_bed(
                        output,
                        ref_val,
                        count_A1=count_A1,
                        force_python_only=force_python_only,
                    )
                    with open_bed(output, count_A1=count_A1) as bed2:
                        assert np.array_equal(
                            bed2.read(
                                dtype="int8", force_python_only=force_python_only
                            ),
                            ref_val,
                        )


def test_write1_bed_f64cpp(tmp_path, shared_datadir):
    with open_bed(shared_datadir / "some_missing.bed") as bed:
        for iid_index in [0, 1, 5]:
            for force_python_only in [False, True]:
                val = bed.read(
                    np.s_[0:iid_index, :],
                    order="F",
                    dtype=np.float64,
                    force_python_only=force_python_only,
                )
                assert val.shape == (iid_index, 100)
                output = str(tmp_path / f"toydata.F64cpp.{iid_index}")
                to_bed(output, val, count_A1=False)
                val2 = open_bed(output, count_A1=False).read(dtype="float64")
                assert np.allclose(val, val2, equal_nan=True)


def test_write1_x_x_cpp(tmp_path, shared_datadir):
    for count_A1 in [False, True]:
        with open_bed(shared_datadir / "some_missing.bed", count_A1=count_A1) as bed:
            for order in ["C", "F"]:
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


def test_respect_read_inputs(shared_datadir):
    ref_val_float = reference_val(shared_datadir)

    ref_val_int8 = ref_val_float.astype("int8")
    ref_val_int8[ref_val_float != ref_val_float] = -127

    with open_bed(shared_datadir / "some_missing.bed") as bed:
        for order in ["F", "C"]:
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


def test_threads(shared_datadir):
    ref_val_float = reference_val(shared_datadir)
    ref_val_int8 = ref_val_float.astype("int8")
    ref_val_int8[ref_val_float != ref_val_float] = -127

    for num_threads in [1, 4]:
        with open_bed(
            shared_datadir / "some_missing.bed", num_threads=num_threads
        ) as bed:
            val = bed.read(dtype="int8")
            assert np.allclose(ref_val_int8, val, equal_nan=True)


def test_write12(tmp_path):
    # ===================================
    #    Starting main function
    # ===================================
    logging.info("starting 'test_writes'")
    np.random.seed(0)
    output_template = str(tmp_path / "writes.{0}.bed")
    i = 0
    for row_count in [0, 5, 2, 1]:
        for col_count in [4, 2, 1, 0]:
            val = np.random.randint(0, 4, size=(row_count, col_count)) * 1.0
            val[val == 3] = np.NaN
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
                except:
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


def test_index(shared_datadir):
    ref_val_float = reference_val(shared_datadir)

    with open_bed(shared_datadir / "some_missing.bed") as bed:
        val = bed.read()
        assert np.allclose(ref_val_float, val, equal_nan=True)

        val = bed.read(2)
        assert np.allclose(ref_val_float[:, [2]], val, equal_nan=True)

        val = bed.read((2))
        assert np.allclose(ref_val_float[:, [2]], val, equal_nan=True)

        val = bed.read((None, 2))
        assert np.allclose(ref_val_float[:, [2]], val, equal_nan=True)

        val = bed.read((1, 2))
        assert np.allclose(ref_val_float[[1], [2]], val, equal_nan=True)

        val = bed.read([2, -2])
        assert np.allclose(ref_val_float[:, [2, -2]], val, equal_nan=True)

        val = bed.read(([1, -1], [2, -2]))
        assert np.allclose(ref_val_float[[1, -1], :][:, [2, -2]], val, equal_nan=True)

        iid_bool = ([False, False, True] * bed.iid_count)[: bed.iid_count]
        sid_bool = ([True, False, True] * bed.sid_count)[: bed.sid_count]
        val = bed.read(sid_bool)
        assert np.allclose(ref_val_float[:, sid_bool], val, equal_nan=True)

        val = bed.read((iid_bool, sid_bool))
        assert np.allclose(ref_val_float[iid_bool, :][:, sid_bool], val, equal_nan=True)

        val = bed.read((1, sid_bool))
        assert np.allclose(ref_val_float[[1], :][:, sid_bool], val, equal_nan=True)

        slicer = np.s_[::2, ::3]
        val = bed.read(slicer[1])
        assert np.allclose(ref_val_float[:, slicer[1]], val, equal_nan=True)

        val = bed.read(slicer)
        assert np.allclose(ref_val_float[slicer], val, equal_nan=True)

        val = bed.read((1, slicer[1]))
        assert np.allclose(ref_val_float[[1], slicer[1]], val, equal_nan=True)


def test_shape(shared_datadir):
    with open_bed(shared_datadir / "plink_sim_10s_100v_10pmiss.bed") as bed:
        assert bed.shape == (10, 100)


def test_zero_files(tmp_path):
    for force_python_only in [False, True]:
        for iid_count in [3, 0]:
            for sid_count in [5, 0]:
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
    val = np.zeros((3, 5), dtype=np.str)
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
        properties2 = deb.properties
        assert np.allclose(val, val2, equal_nan=True)
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


def test_env(shared_datadir):
    if platform.system() == "Darwin":
        return

    key = "MKL_NUM_THREADS"
    original_val = os.environ.get(key)
    try:
        os.environ[key] = "1"
        with open_bed(shared_datadir / "some_missing.bed") as bed:
            _ = bed.read(np.s_[:100, :100])
        os.environ[key] = "10"
        with open_bed(shared_datadir / "some_missing.bed") as bed:
            _ = bed.read(np.s_[:100, :100])
        os.environ[key] = "BADVALUE"
        with pytest.raises(ValueError):
            with open_bed(shared_datadir / "some_missing.bed") as bed:
                _ = bed.read(np.s_[:100, :100])
    finally:
        if original_val is None:
            if key in os.environ:
                del os.environ[key]
        else:
            os.environ[key] = original_val


#if __name__ == "__main__":
#    logging.basicConfig(level=logging.INFO)

#    shared_datadir = Path(r"D:\OneDrive\programs\bed-reader\bed_reader\tests\data")
#    tmp_path = Path(r"m:/deldir/tests")
#    test_write12(tmp_path)
#    pytest.main([__file__])
