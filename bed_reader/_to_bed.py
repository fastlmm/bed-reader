import logging
import os
from pathlib import Path
from typing import Any, List, Mapping, Optional, Union

import numpy as np

from bed_reader import open_bed


#!!!cmk say something about support for snp-minor vs major
def to_bed(
    filepath: Union[str, Path],
    val: np.ndarray,
    metadata: Mapping[str, List[Any]] = {},
    count_A1: bool = True,
    force_python_only: bool = False,
):
    """
    Write values to a file in PLINK BED format.

    Parameters
    ----------
    filepath:
        Bed file path #!!!cmk Bed vs BED vs cmkstar.bed
    val: array-like:
        A two-dimension array (or array-like object) of values. The values should
        be (or be convertable to) all floats or all integers. The values should be 0, 1, 2, or missing.
        If floats, missing is ``np.nan``. If integers, missing is -127.
    metadata: dict, optional
        A dictionary of metadata of interest. The default is an empty dictionary.
        The keys of the dictionary are the names of the metadata of interest.
        The possible keys are:

             "fid" (family id), "iid" (individual or sample id), "father" (father id),
             "mother" (mother id), "sex", "pheno" (phenotype), "chromosome", "sid"
             (SNP or variant id), "cm_position" (centimorgan position), "bp_position"
             (base-pair position), "allele_1", "allele_2".
            
        The values are lists or arrays. Any metadata not given will be filled in with
        a default value. CMK see example
    count_A1: bool, optional
        True (default) to count the number of A1 alleles (the PLINK standard). False to count the number of A2 alleles.
    force_python_only
        If False (default), uses the faster C++ code; otherwise it uses the slower pure Python code.

    Examples
    --------

    .. doctest::

    In this example, full metadata is given.

        >>> import numpy as np
        >>> from bed_reader import to_bed, tmp_path
        >>>
        >>> output_file = tmp_path() / "small.bed"
        >>> val = [[1.0, 0.0, np.nan, 0.0], [2.0, 0.0, np.nan, 2.0], [0.0, 1.0, 2.0, 0.0]]
        >>> metadata = {
        ...    "fid": ["fid1", "fid1", "fid2"],
        ...    "iid": ["iid1", "iid2", "iid3"],
        ...    "father": ["iid23", "iid23", "iid22"],
        ...    "mother": ["iid34", "iid34", "iid33"],
        ...    "sex": [1, 2, 0],
        ...    "pheno": ["red", "red", "blue"],
        ...    "chromosome": ["1", "1", "5", "Y"],
        ...    "sid": ["sid1", "sid2", "sid3", "sid4"],
        ...    "cm_position": [100.4, 2000.5, 4000.7, 7000.9],
        ...    "bp_position": [1, 100, 1000, 1004],
        ...    "allele_1": ["A", "T", "A", "T"],
        ...    "allele_2": ["A", "C", "C", "G"],
        ... }
        >>> to_bed(output_file, val, metadata=metadata)


    .. doctest::

    In this example, no metadata is given, so default values as assigned.
    If we read the file and list the chromosome, it will be a list of '0's
    the default chromosome value. cmkref

        >>> output_file2 = tmp_path() / "small2.bed"
        >>> val = [[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]]
        >>> to_bed(output_file2, val)
        >>>
        >>> from bed_reader import open_bed
        >>> with open_bed(output_file2) as bed2:
        ...     print(bed2.chromosome)
        ['0' '0' '0' '0']
        
    """
    filepath = Path(filepath)

    val = _fix_up_val(val)
    iid_count = val.shape[0]
    sid_count = val.shape[1]

    metadata, _ = open_bed._fix_up_metadata(
        metadata, iid_count=iid_count, sid_count=sid_count, use_fill_sequence=True
    )

    open_bed._write_fam_or_bim(filepath, metadata, "fam")
    open_bed._write_fam_or_bim(filepath, metadata, "bim")

    bedfile = str(
        open_bed._name_of_other_file(filepath, remove_suffix="bed", add_suffix="bed")
    ).encode("ascii")

    if not force_python_only:
        from bed_reader import wrap_plink_parser_onep

        if val.flags["C_CONTIGUOUS"]:
            order = "C"
        elif val.flags["F_CONTIGUOUS"]:
            order = "F"
        else:
            raise ValueError(f"val must be contiguous.")

        iid_count, sid_count = val.shape
        try:
            if val.dtype == np.float64:
                if order == "F":
                    wrap_plink_parser_onep.writePlinkBedFile2doubleFAAA(
                        bedfile, iid_count, sid_count, count_A1, val,
                    )
                else:
                    wrap_plink_parser_onep.writePlinkBedFile2doubleCAAA(
                        bedfile, iid_count, sid_count, count_A1, val,
                    )
            elif val.dtype == np.float32:
                if order == "F":
                    wrap_plink_parser_onep.writePlinkBedFile2floatFAAA(
                        bedfile, iid_count, sid_count, count_A1, val,
                    )
                else:
                    wrap_plink_parser_onep.writePlinkBedFile2floatCAAA(
                        bedfile, iid_count, sid_count, count_A1, val,
                    )
            elif val.dtype == np.int8:
                if order == "F":
                    wrap_plink_parser_onep.writePlinkBedFile2int8FAAA(
                        bedfile, iid_count, sid_count, count_A1, val,
                    )
                else:
                    wrap_plink_parser_onep.writePlinkBedFile2int8CAAA(
                        bedfile, iid_count, sid_count, count_A1, val,
                    )
            else:
                raise ValueError(
                    f"dtype '{val.dtype}' not known, only 'int8', 'float32', and 'float64' are allowed."
                )
        except SystemError as system_error:
            try:
                os.unlink(bedfile)
            except Exception:
                pass
            raise system_error.__cause__
    else:
        if not count_A1:
            zero_code = 0b00
            two_code = 0b11
        else:
            zero_code = 0b11
            two_code = 0b00

        with open(bedfile, "wb") as bed_filepointer:
            # see http://zzz.bwh.harvard.edu/plink/binary.shtml
            bed_filepointer.write(bytes(bytearray([0b01101100])))  # magic numbers
            bed_filepointer.write(bytes(bytearray([0b00011011])))  # magic numbers
            bed_filepointer.write(bytes(bytearray([0b00000001])))  # snp major

            for sid_index in range(sid_count):
                if sid_index % 1 == 0:
                    logging.info(
                        "Writing snp # {0} to file '{1}'".format(sid_index, filepath)
                    )

                col = val[:, sid_index]
                for iid_by_four in range(0, iid_count, 4):
                    vals_for_this_byte = col[iid_by_four : iid_by_four + 4]
                    byte = 0b00000000
                    for val_index in range(len(vals_for_this_byte)):
                        val_for_byte = vals_for_this_byte[val_index]
                        if val_for_byte == 0:
                            code = zero_code
                        elif val_for_byte == 1:
                            code = 0b10  # backwards on purpose
                        elif val_for_byte == 2:
                            code = two_code
                        elif (
                            val.dtype == np.int8 and val_for_byte == -127
                        ) or np.isnan(val_for_byte):
                            code = 0b01  # backwards on purpose
                        else:
                            raise ValueError(
                                "Attempt to write illegal value to BED file. Only 0,1,2,missing allowed."
                            )
                        byte |= code << (val_index * 2)
                    bed_filepointer.write(bytes(bytearray([byte])))
    logging.info(f"Done writing {filepath}")


def _fix_up_val(input):

    if not isinstance(input, np.ndarray):
        return _fix_up_val(np.array(input))

    if np.issubdtype(input.dtype, np.integer) and input.dtype != np.int8:
        return _fix_up_val(np.array(input, dtype=np.int8))
    elif np.issubdtype(input.dtype, np.float) and input.dtype not in (
        np.float32,
        np.float64,
    ):
        return _fix_up_val(np.array(input, dtype=np.float32))

    if len(input.shape) != 2:
        raise ValueError("val should be two dimensional")

    return input


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import pytest

    pytest.main(["--doctest-modules", __file__])
