from pathlib import Path
import numpy as np
import logging
import os
from bed_reader import open_bed

#!!!cmk say something about support for snp-minor vs major
def to_bed(
    filepath, val, metadata={}, count_A1=True, force_python_only=False,
):
    """
    !!!cmk need doc string
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
        return open_bed._fix_up_val(np.array(input))

    if len(input.shape) != 2:
        raise ValueError("val should be two dimensional")

    return input
