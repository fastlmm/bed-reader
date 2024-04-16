import logging
import os
from pathlib import Path
from typing import Any, List, Mapping, Union
import numpy as np
from bed_reader import get_num_threads, open_bed
from .bed_reader import write_f32, write_f64, write_i8, encode1_i8  # type: ignore


class create_bed:
    """
    cmk update docs
    Write values to a file in PLINK .bed format.

    Parameters
    ----------
    filepath:
        .bed file to write to.
    val: array-like:
        A two-dimension array (or array-like object) of values. The values should
        be (or be convertible to) all floats or all integers.
        The values should be 0, 1, 2, or missing.
        If floats, missing is ``np.nan``. If integers, missing is -127.
    properties: dict, optional
        A dictionary of property names and values to write to the .fam and .bim files.
        Any properties not mentioned will be filled in with default values.

        The possible property names are:

             "fid" (family id), "iid" (individual or sample id), "father" (father id),
             "mother" (mother id), "sex", "pheno" (phenotype), "chromosome", "sid"
             (SNP or variant id), "cm_position" (centimorgan position), "bp_position"
             (base-pair position), "allele_1", "allele_2".

         The values are lists or arrays. See example, below.
    count_A1: bool, optional
        True (default) to count the number of A1 alleles (the PLINK standard).
        False to count the number of A2 alleles.
    fam_filepath: pathlib.Path or str, optional
        Path to the file containing information about each individual (sample).
        Defaults to replacing the .bed file’s suffix with .fam.
    bim_filepath: pathlib.Path or str, optional
        Path to the file containing information about each SNP (variant).
        Defaults to replacing the .bed file’s suffix with .bim.
    force_python_only
        If False (default), uses the faster Rust code; otherwise it uses the slower
        pure Python code.

    num_threads: None or int, optional
        The number of threads with which to write data.
        Defaults to all available processors.
        Can also be set with these environment variables (listed in priority order):
        'PST_NUM_THREADS', 'NUM_THREADS', 'MKL_NUM_THREADS'.


    Examples
    --------

    In this example, all properties are given.

    .. doctest::

        >>> import numpy as np
        >>> from bed_reader import to_bed, tmp_path
        >>>
        >>> output_file = tmp_path() / "small.bed"
        >>> val = [[1.0, 0.0, np.nan, 0.0],
        ...        [2.0, 0.0, np.nan, 2.0],
        ...        [0.0, 1.0, 2.0, 0.0]]
        >>> properties = {
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
        >>> to_bed(output_file, val, properties=properties)

    Here, no properties are given, so default values are assigned.
    If we then read the new file and list the chromosome property,
    it is an array of '0's, the default chromosome value.

    .. doctest::

        >>> output_file2 = tmp_path() / "small2.bed"
        >>> val = [[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]]
        >>> to_bed(output_file2, val)
        >>>
        >>> from bed_reader import open_bed
        >>> with open_bed(output_file2) as bed2:
        ...     print(bed2.chromosome)
        ['0' '0' '0' '0']

    """

    def __init__(
        self,
        location: Union[str, Path],
        iid_count: int,
        sid_count: int,
        properties: Mapping[str, List[Any]] = {},
        major: str = "SNP",
        count_A1: bool = True,
        fam_location: Union[str, Path] = None,
        bim_location: Union[str, Path] = None,
        force_python_only: bool = False,
        num_threads=None,
    ):
        # cmk set this well
        buffer_size = 8192_000
        self.location = Path(location)
        self.iid_count = iid_count
        self.sid_count = sid_count
        self.count_A1 = count_A1
        self.force_python_only = force_python_only
        self.num_threads = num_threads or os.cpu_count()

        if major == "SNP":
            self.major_count = sid_count
            self.minor_count = iid_count
        elif major == "individual":
            self.major_count = iid_count
            self.minor_count = sid_count
        else:
            raise ValueError(f"major must be 'SNP' or 'individual', not '{major}'")
        self.major = major
        self.minor_count_div4 = (self.minor_count - 1) // 4 + 1
        self.buffer = np.zeros(self.minor_count_div4, dtype=np.uint8)

        if not count_A1:
            self.zero_code = 0b00
            self.two_code = 0b11
        else:
            self.zero_code = 0b11
            self.two_code = 0b00

        fam_location = (
            Path(fam_location)
            if fam_location is not None
            else self._replace_extension(self.location, "fam")
        )
        bim_location = (
            Path(bim_location)
            if bim_location is not None
            else self._replace_extension(self.location, "bim")
        )

        properties, _ = open_bed._fix_up_properties(
            properties,
            iid_count,
            sid_count,
            use_fill_sequence=True,
        )

        open_bed._write_fam_or_bim(self.location, properties, "fam", fam_location)
        open_bed._write_fam_or_bim(self.location, properties, "bim", bim_location)

        self.temp_filepath = self.location.with_suffix(".bed_temp")
        self.file_pointer = open(self.temp_filepath, "wb", buffering=buffer_size)
        # see http://zzz.bwh.harvard.edu/plink/binary.shtml
        self.file_pointer.write(bytes(bytearray([0b01101100])))  # magic numbers
        self.file_pointer.write(bytes(bytearray([0b00011011])))  # magic numbers
        if self.major == "SNP":
            self.file_pointer.write(bytes(bytearray([0b00000001])))  # snp major
        else:
            assert self.major == "individual"  # real assert
            self.file_pointer.write(bytes(bytearray([0b00000000])))
        self.write_count = 0

        if os.path.exists(self.location):
            os.unlink(self.location)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        if self.file_pointer is not None:
            try:
                self.file_pointer.close()
                self.file_pointer = None
                if self.write_count != self.major_count:
                    raise ValueError(
                        f"Attempt to write fewer vectors ({self.write_count}) than expected ({self.major_count})"
                    )
                os.rename(self.temp_filepath, self.location)
            except Exception as e:
                self._clean_up_error()
                raise e
        else:
            self._clean_up_error()

    def _clean_up_error(self):
        try:
            os.unlink(self.temp_filepath)
        except Exception:
            pass

    @staticmethod
    def _replace_extension(location, extension):
        if open_bed._is_url(location):
            return location.parent / (location.stem + "." + extension)

    def write(self, vector):
        if self.file_pointer is None:
            raise RuntimeError("Attempt to write after file was closed for writing.")

        try:
            self.write_count += 1
            if self.write_count > self.major_count:
                raise ValueError(
                    f"Attempt to write more vectors ({self.write_count}) than expected ({self.major_count})"
                )

            vector = _fix_up_vector(vector)
            if len(vector) != self.minor_count:
                raise ValueError(
                    f"Expected vector to {self.minor_count} values, got {len(vector)} instead"
                )

            self._internal_write(vector)

        except Exception as e:
            self.file_pointer.close()
            self.file_pointer = None  # Prevent further writes
            raise e  # Re-raise the exception to handle it externally

    def _internal_write(self, vector):
        if not self.force_python_only and vector.dtype == np.int8:  # cmk
            if vector.dtype == np.int8:
                encode1_i8(self.count_A1, vector, self.buffer, self.num_threads)
            else:
                raise ValueError(
                    f"dtype '{vector.dtype}' not known, only 'int8' is allowed so far cmk."
                )
            self.file_pointer.write(self.buffer)
        else:
            for minor_by_four in range(0, self.minor_count, 4):
                vals_for_this_byte = vector[minor_by_four : minor_by_four + 4]
                byte = 0b00000000
                for val_index in range(len(vals_for_this_byte)):
                    val_for_byte = vals_for_this_byte[val_index]
                    if val_for_byte == 0:
                        code = self.zero_code
                    elif val_for_byte == 1:
                        code = 0b10  # backwards on purpose
                    elif val_for_byte == 2:
                        code = self.two_code
                    elif (vector.dtype == np.int8 and val_for_byte == -127) or np.isnan(
                        val_for_byte
                    ):
                        code = 0b01  # backwards on purpose
                    else:
                        raise ValueError(
                            "Attempt to write illegal value to .bed file. "
                            + "Only 0,1,2,missing allowed."
                        )
                    byte |= code << (val_index * 2)
                # cmk yikes, this is writing one byte at a time
                self.file_pointer.write(bytes(bytearray([byte])))


def _fix_up_vector(input):
    if not isinstance(input, np.ndarray):
        return _fix_up_vector(np.array(input))

    if np.issubdtype(input.dtype, np.integer) and input.dtype != np.int8:
        return _fix_up_vector(np.array(input, dtype=np.int8))
    elif np.issubdtype(input.dtype, np.floating) and input.dtype not in (
        np.float32,
        np.float64,
    ):
        return _fix_up_vector(np.array(input, dtype=np.float32))

    if len(input.shape) != 1:
        raise ValueError("vector should be one dimensional")

    return input


def to_bed(
    filepath: Union[str, Path],
    val: np.ndarray,
    properties: Mapping[str, List[Any]] = {},
    count_A1: bool = True,
    fam_filepath: Union[str, Path] = None,
    bim_filepath: Union[str, Path] = None,
    force_python_only: bool = False,
    num_threads=None,
):
    """
    Write values to a file in PLINK .bed format.

    Parameters
    ----------
    filepath:
        .bed file to write to.
    val: array-like:
        A two-dimension array (or array-like object) of values. The values should
        be (or be convertible to) all floats or all integers.
        The values should be 0, 1, 2, or missing.
        If floats, missing is ``np.nan``. If integers, missing is -127.
    properties: dict, optional
        A dictionary of property names and values to write to the .fam and .bim files.
        Any properties not mentioned will be filled in with default values.

        The possible property names are:

             "fid" (family id), "iid" (individual or sample id), "father" (father id),
             "mother" (mother id), "sex", "pheno" (phenotype), "chromosome", "sid"
             (SNP or variant id), "cm_position" (centimorgan position), "bp_position"
             (base-pair position), "allele_1", "allele_2".

         The values are lists or arrays. See example, below.
    count_A1: bool, optional
        True (default) to count the number of A1 alleles (the PLINK standard).
        False to count the number of A2 alleles.
    fam_filepath: pathlib.Path or str, optional
        Path to the file containing information about each individual (sample).
        Defaults to replacing the .bed file’s suffix with .fam.
    bim_filepath: pathlib.Path or str, optional
        Path to the file containing information about each SNP (variant).
        Defaults to replacing the .bed file’s suffix with .bim.
    force_python_only
        If False (default), uses the faster Rust code; otherwise it uses the slower
        pure Python code.

    num_threads: None or int, optional
        The number of threads with which to write data.
        Defaults to all available processors.
        Can also be set with these environment variables (listed in priority order):
        'PST_NUM_THREADS', 'NUM_THREADS', 'MKL_NUM_THREADS'.


    Examples
    --------

    In this example, all properties are given.

    .. doctest::

        >>> import numpy as np
        >>> from bed_reader import to_bed, tmp_path
        >>>
        >>> output_file = tmp_path() / "small.bed"
        >>> val = [[1.0, 0.0, np.nan, 0.0],
        ...        [2.0, 0.0, np.nan, 2.0],
        ...        [0.0, 1.0, 2.0, 0.0]]
        >>> properties = {
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
        >>> to_bed(output_file, val, properties=properties)

    Here, no properties are given, so default values are assigned.
    If we then read the new file and list the chromosome property,
    it is an array of '0's, the default chromosome value.

    .. doctest::

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

    properties, _ = open_bed._fix_up_properties(
        properties, iid_count=iid_count, sid_count=sid_count, use_fill_sequence=True
    )

    open_bed._write_fam_or_bim(filepath, properties, "fam", fam_filepath)
    open_bed._write_fam_or_bim(filepath, properties, "bim", bim_filepath)

    if not force_python_only:
        if not val.flags["C_CONTIGUOUS"] and not val.flags["F_CONTIGUOUS"]:
            raise ValueError("val must be contiguous.")

        num_threads = get_num_threads(num_threads)

        iid_count, sid_count = val.shape
        try:
            if val.dtype == np.float64:
                write_f64(
                    str(filepath),
                    is_a1_counted=count_A1,
                    val=val,
                    num_threads=num_threads,
                )
            elif val.dtype == np.float32:
                write_f32(
                    str(filepath),
                    is_a1_counted=count_A1,
                    val=val,
                    num_threads=num_threads,
                )
            elif val.dtype == np.int8:
                write_i8(
                    str(filepath),
                    is_a1_counted=count_A1,
                    val=val,
                    num_threads=num_threads,
                )
            else:
                raise ValueError(
                    f"dtype '{val.dtype}' not known, only "
                    + "'int8', 'float32', and 'float64' are allowed."
                )
        except SystemError as system_error:
            try:
                os.unlink(filepath)
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

        # cmk change to buffer writing
        with open(filepath, "wb") as bed_filepointer:
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
                                "Attempt to write illegal value to .bed file. "
                                + "Only 0,1,2,missing allowed."
                            )
                        byte |= code << (val_index * 2)
                    bed_filepointer.write(bytes(bytearray([byte])))
    logging.info(f"Done writing {filepath}")


def _fix_up_val(input):
    if not isinstance(input, np.ndarray):
        return _fix_up_val(np.array(input))

    if np.issubdtype(input.dtype, np.integer) and input.dtype != np.int8:
        return _fix_up_val(np.array(input, dtype=np.int8))
    elif np.issubdtype(input.dtype, np.floating) and input.dtype not in (
        np.float32,
        np.float64,
    ):
        return _fix_up_val(np.array(input, dtype=np.float32))

    if len(input.shape) != 2:
        raise ValueError("val should be two dimensional")

    return input


# if __name__ == "__main__":
#    logging.basicConfig(level=logging.INFO)

#    import pytest

#    pytest.main(["--doctest-modules", __file__])
