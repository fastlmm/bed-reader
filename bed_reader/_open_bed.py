#!!!cmk todo: Offer to ignore some or all fam bim fields
#!!!cmk add typing info
#!!!cmk run flake8, isort, etc
import os
import numpy as np
import numbers
import pandas as pd
import logging
from pathlib import Path
import multiprocessing
from dataclasses import dataclass
import sys
import platform

import math
from typing import Any, List, Optional, Tuple, Union

from itertools import takewhile, repeat

# https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python
def _rawincount(filepath):
    with open(filepath, "rb") as f:
        bufgen = takewhile(lambda x: x, (f.raw.read(1024 * 1024) for _ in repeat(None)))
        return sum(buf.count(b"\n") for buf in bufgen)


@dataclass
class _MetaMeta:
    suffix : str
    column : int
    dtype: type
    missing_value: object
    fill_sequence: object

def _all_same(key, length, missing, dtype):
    if np.issubdtype(dtype, np.str_):
        dtype = f"<U{len(missing)}"
    return np.full(length, missing, dtype=dtype)

def _sequence(key, length, missing, dtype):
    if np.issubdtype(dtype, np.str_):
        longest = len(f"{key}{length}")
        dtype = f"<U{longest}"
    return np.fromiter(
        (f"{key}{i+1}" for i in range(length)), dtype=dtype, count=length
    )

_delimiters = {"fam": r"\s+", "bim": "\t"}
_count_name = {"fam": "iid_count", "bim": "sid_count"}


_meta_meta = {
    # https://stackoverflow.com/questions/41921255/staticmethod-object-is-not-callable
    "fid": _MetaMeta("fam", 0, np.str_, "0", _all_same),
    "iid": _MetaMeta("fam", 1, np.str_, None, _sequence),
    "father": _MetaMeta("fam", 2, np.str_, "0", _all_same),
    "mother": _MetaMeta("fam", 3, np.str_, "0", _all_same),
    "sex": _MetaMeta("fam", 4, np.int32, 0, _all_same),
    "pheno": _MetaMeta("fam", 5, np.str_, "0", _all_same),
    "chromosome": _MetaMeta("bim", 0, np.str_, "0", _all_same),
    "sid": _MetaMeta("bim", 1, np.str_, None, _sequence),
    "cm_position": _MetaMeta("bim", 2, np.float32, 0, _all_same),
    "bp_position": _MetaMeta("bim", 3, np.int32, 0, _all_same),
    "allele_1": _MetaMeta("bim", 4, np.str_, "A1", _all_same),
    "allele_2": _MetaMeta("bim", 5, np.str_, "A2", _all_same),
    }

class open_bed:  #!!!cmk need doc strings everywhere
    """
    A NumPy-inspired class for fast opening and reading of PLINK \*.bed files.

    Parameters
    ----------
    filepath: string or path
        BED file path
    iid_count: int or ``none``
        Number of individuals (samples) in the BED file.
        Defaults to quickly counting the number itself.
    sid_count: int or ``none``
        Number of SNPs (variants) in the BED file.
        Defaults to quickly counting the number itself.
    metadata: dictionary
        A dictionary of replacement metadata. By default, this
        is empty and any queries about metadata
        are answered by reading either the \*.fam or the *.\bim file. Any replacement
        metadata given here takes precedent over the information in the file.
        The keys of the dictionary are the names of the metadata, specifically,
        "fid" (family id), "iid" (individual or sample id), "father" (father id),
        "mother" (mother id), "sex", "pheno" (phenotype), "chromosome", "sid"
        (SNP or variant id), "cm_position" (centimorgan position), "bp_position"
        (base-pair position), "allele_1", "allele_2". The values in the dictionary
        are list or array.
    count_A1: bool
        Tells if the reader should count the number of A1 alleles (the PLINK standard and the default) or the number of A2 alleles.
    num_threads: int
        Tells how many threads to use to read data. Defaults to all available threads.
        Can also be set the 'MKL_NUM_THREADS' environment variable.
    skip_format_check: b
        If False (default), will immediately check that '.bed' file has expected starting bytes
        on open. If True, will not check until (and if) the file is read.
        
    Returns
    -------
    an open_bed object : :class:`open_bed`

    .. _open_examples:

    Examples
    --------
    #!!!cmk give examples of metadata
    #!!!cmk talk about missing data
    #!!!cmk talk about multithreading
    With the `with <https://docs.python.org/3/reference/compound_stmts.html#grammar-token-with-stmt>`__ statement, list individual (sample) :attr:`iid` and SNP (variant) :attr:`sid`, then :meth:`read` the whole file.

    .. doctest::

        >>> from bed_reader import open_bed
        >>>
        >>> with open_bed("distributed_bed_test1_X.bed") as bed:
        ...     print(bed.iid)
        ...     print(bed.sid)
        ...     print(bed.read())
        ['SNP1' 'SNP2' 'SNP3' 'SNP4']
        ['sample_0' 'sample_1' 'sample_2' 'sample_3']
        [[[1. 0. 1. 0.]
          [0. 1. 1. 0.]
          [1. 0. 0. 1.]
          [0. 1. 0. 1.]]
        <BLANKLINE>
         [[0. 1. 1. 0.]
          [1. 0. 0. 1.]
          [0. 1. 0. 1.]
          [1. 0. 1. 0.]]
        <BLANKLINE>
         [[1. 0. 0. 1.]
          [0. 1. 0. 1.]
          [1. 0. 1. 0.]
          [0. 1. 1. 0.]]
        <BLANKLINE>
         [[0. 1. 0. 1.]
          [1. 0. 1. 0.]
          [0. 1. 1. 0.]
          [1. 0. 0. 1.]]]

    Open the file (without `with`) and read probabilities for one variant.

    .. doctest::

        >>> bed = open_bed("distributed_bed_test1_X.bed")
        >>> print(bed.read(2))
        [[[1. 0. 0. 1.]]
        <BLANKLINE>
         [[0. 1. 0. 1.]]
        <BLANKLINE>
         [[1. 0. 1. 0.]]
        <BLANKLINE>
         [[0. 1. 1. 0.]]]
        >>> del bed                 # close and delete object

    Open the file and then first read for a :class:`slice` of samples and variants, and then for a single sample and variant.

    .. doctest::

        >>> bed = open_bed(file, verbose=False)
        >>> print(bed.read((slice(1,3),slice(2,4))))
        [[[0. 1. 0. 1.]
          [1. 0. 1. 0.]]
        <BLANKLINE>
         [[1. 0. 1. 0.]
          [0. 1. 1. 0.]]]
        >>> print(bed.read((0,1)))
        [[[0. 1. 1. 0.]]]
        >>> del bed                 # close and delete object


        #!!!cmk need example of accessing the metadata
        #!!!cmk need exmaple of overriding the metadata
        #!!!cmk in README say: Documentation not API Documatnion

    .. _sample format: https://www.well.ox.ac.uk/~gav/qctool/documentation/sample_file_formats.html #!!!cmk

    """
    def __init__(
        self,
        filepath,
        iid_count=None,
        sid_count=None,
        metadata={},
        count_A1=True,
        num_threads=None,
        skip_format_check=False,
    ):  #!!!document these new optionals. they are here
        self.filepath = Path(filepath)
        self.count_A1 = count_A1
        self._num_threads = num_threads
        self.skip_format_check = skip_format_check

        self.metadata_dict, self._counts = open_bed._fixup_metadata(
            metadata, iid_count, sid_count, use_fill_sequence=False
        )
        self._iid_range = None
        self._sid_range = None

        if not self.skip_format_check:
            bedfile = open_bed._name_of_other_file(self.filepath, "bed", "bed")
            with open(bedfile, "rb") as filepointer:
                self._check_file(filepointer)


    @staticmethod
    def _fixup_metadata(metadata, iid_count, sid_count, use_fill_sequence):

        metadata_dict = {key: None for key in _meta_meta}
        count_dict = {"fam": iid_count, "bim": sid_count}

        for key, input in metadata.items():
            if key not in _meta_meta:
                raise KeyError(f"metadata key '{key}' not known")

        for key, mm in _meta_meta.items():
            count = count_dict[mm.suffix]

            input = metadata.get(key)
            if input is None:
                if use_fill_sequence:
                    output = mm.fill_sequence(key, count, mm.missing_value, mm.dtype)
                else:
                    continue
            elif len(input) == 0:
                output = np.zeros([0], dtype=mm.dtype)
            else:
                if not isinstance(input, np.ndarray) or not np.issubdtype(
                    input.dtype, mm.dtype
                ):
                    input = np.array(input, dtype=mm.dtype)
                if len(input.shape) != 1:
                    raise ValueError(f"Override {key} should be one dimensional")
                output = input

            if count is None:
                count_dict[mm.suffix] = len(output)
            else:
                if count != len(output):
                    raise ValueError(
                        f"The length of override {key}, {len(output)}, should not be different from the current {_count_name[mm.suffix]}, {count}"
                    )
            metadata_dict[key] = output
        return metadata_dict, count_dict 

    def _read_fam_or_bim(self, suffix):
        metafile = open_bed._name_of_other_file(self.filepath, "bed", suffix)
        logging.info("Loading {0} file {1}".format(suffix, metafile))

        count = self._counts[suffix]

        delimiter = _delimiters[suffix]
        if delimiter in {r"\s+"}:
            delimiter = None
            delim_whitespace = True
        else:
            delim_whitespace = False

        if os.path.getsize(metafile) == 0:
            fields = []
        else:
            fields = pd.read_csv(
                metafile,
                delimiter=delimiter,
                delim_whitespace=delim_whitespace,
                header=None,
                index_col=False,
                comment=None,
            )

        if count is None:
            self._counts[suffix] = len(fields)
        else:
            if count != len(fields):
                raise ValueError(
                    f"The number of lines in the *.{suffix} file, {len(fields)}, should not be different from the current {_count_name[suffix]}, {count}"
                )
        for key, mm in _meta_meta.items():
            if mm.suffix is not suffix:
                continue
            val = self.metadata_dict[key]
            if val is None:
                if len(fields) == 0:
                    output = np.array([], dtype=mm.dtype)
                elif mm.missing_value is None:
                    output = np.array(fields[mm.column], dtype=mm.dtype)
                else:
                    output = np.array(fields[mm.column].fillna(mm.missing_value), dtype=mm.dtype)
                self.metadata_dict[key] = output

    @staticmethod
    def _name_of_other_file(filepath, remove_suffix, add_suffix):
        if filepath.suffix.lower() == "." + remove_suffix:
            filepath = filepath.parent / filepath.stem
        return filepath.parent / (filepath.name + "." + add_suffix)

    def __str__(self):
        return f"{self.__class__.__name__}('{self.filepath}',...)"

    #!!!cmk make sure these are in a good order for the documentation

    @property
    def fid(self):
        """
        The family id (a :class:`numpy.ndarray` of ``str``).

        #!!!cmk tell that if needed will open and reader the *.fam file

        Example
        -------
        .. doctest::

            >>> from bed_reader import open_bed
            >>>
            >>> file = example_filepath("haplotypes.bed")
            >>> with open_bed(file) as bed:
            ...     print(bed.fid)
            ['sample_0' 'sample_1' 'sample_2' 'sample_3']

        """

        return self.metadata_item("fid")

    @property
    def iid(self):
        '''
        !!!cmk need doc string
        '''
        return self.metadata_item("iid")

    @property
    def father(self):
        '''
        !!!cmk need doc string
        '''
        return self.metadata_item("father")

    @property
    def mother(self):
        '''
        !!!cmk need doc string
        '''
        return self.metadata_item("mother")

    @property
    def sex(self):
        '''
        !!!cmk need doc string
        '''
        return self.metadata_item("sex")

    @property
    def pheno(self):
        '''
        !!!cmk need doc string
        '''
        return self.metadata_item("pheno")

    @property
    def metadata(self):
        '''
        !!!cmk need doc string
        !!!cmk tell that if needed, will open and read *.fam and *.bim files
        '''
        for key in _meta_meta:
            self.metadata_item(key)
        return self.metadata_dict

    def metadata_item(self, key):
        '''
        !!!cmk need doc string
        '''
        val = self.metadata_dict.get(key)
        if val is None:
            mm = _meta_meta[key]
            self._read_fam_or_bim(suffix=mm.suffix)
            return self.metadata_dict[key]
        else:
            return val

    @property
    def chromosome(self):
        '''
        !!!cmk need doc string
        '''
        return self.metadata_item("chromosome")

    @property
    def sid(self):
        '''
        !!!cmk need doc string
        '''
        return self.metadata_item("sid")

    @property
    def cm_position(self):
        '''
        !!!cmk need doc string
        '''
        return self.metadata_item("cm_position")

    @property
    def bp_position(self):
        '''
        !!!cmk need doc string
        '''
        return self.metadata_item("bp_position")

    @property
    def allele_1(self):
        '''
        !!!cmk need doc string
        '''
        return self.metadata_item("allele_1")

    @property
    def allele_2(self):
        '''
        !!!cmk need doc string
        '''
        return self.metadata_item("allele_2")

    @property
    def iid_count(self):
        '''
        !!!cmk need doc string
        '''
        return self._count("fam")

    @property
    def sid_count(self):
        '''
        !!!cmk need doc string
        '''
        return self._count("bim")

    def _count(self, suffix):
        count = self._counts[suffix]
        if count is None:
            metafile = open_bed._name_of_other_file(self.filepath, "bed", suffix)
            count = _rawincount(metafile)
            self._counts[suffix] = count
        return count

    @staticmethod
    def _check_file(filepointer):
        mode = filepointer.read(2)
        if mode != b"l\x1b":
            raise ValueError(
                "No valid binary BED file"
            )
        mode = filepointer.read(1)  # \x01 = SNP major \x00 = individual major
        if mode != b"\x01":
            raise ValueError(
                "only SNP-major is implemented"
            )  #!!!cmk should mention this

    def __del__(self):
        self.__exit__()

    def close(self):
        """
        !!!cmk doc this
            """
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        pass

    @staticmethod
    def _find_openmp():
        if "bed_reader.wrap_plink_parser_openmp" in sys.modules:
            return
        if platform.system() == "Windows":
            #print("cmk in windows _find_openmp")
            from ctypes import cdll
            from ctypes.util import find_library
            dllname = "libiomp5md.dll"
            location_list = [Path(__file__).parent / dllname, Path(__file__).parent.parent / "external/intel/windows/compiler/lib/intel64" / dllname]
            for location in location_list:
                #print(f"cmk looking for '{dllname}' at '{location}'")
                if location.exists():
                    #print(f"cmk found it")
                    cdll.LoadLibrary(str(location))
                    #print(f"cmk loaded it")
                    return
            raise Exception(f"Can't find '{dllname}'")


    #!!!cmk say something about support for snp-minor vs major
    @staticmethod
    def write(
        filepath, val, metadata={}, count_A1=True, force_python_only=False,
    ):
        '''
        !!!cmk need doc string
        '''
        filepath = Path(filepath)
        iid_count = val.shape[0]
        sid_count = val.shape[1]

        metadata, _ = open_bed._fixup_metadata(
            metadata, iid_count=iid_count, sid_count=sid_count, use_fill_sequence=True
        )

        open_bed._write_fam_or_bim(filepath, metadata, "fam")
        open_bed._write_fam_or_bim(filepath, metadata, "bim")

        bedfile = str(open_bed._name_of_other_file(
            filepath, remove_suffix="bed", add_suffix="bed"
        )).encode("ascii")

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
                            bedfile,
                            iid_count,
                            sid_count,
                            count_A1,
                            val,
                        )
                    else:
                        wrap_plink_parser_onep.writePlinkBedFile2doubleCAAA(
                            bedfile,
                            iid_count,
                            sid_count,
                            count_A1,
                            val,
                        )
                elif val.dtype == np.float32:
                    if order == "F":
                        wrap_plink_parser_onep.writePlinkBedFile2floatFAAA(
                            bedfile,
                            iid_count,
                            sid_count,
                            count_A1,
                            val,
                        )
                    else:
                        wrap_plink_parser_onep.writePlinkBedFile2floatCAAA(
                            bedfile,
                            iid_count,
                            sid_count,
                            count_A1,
                            val,
                        )
                elif val.dtype == np.int8:
                    if order == "F":
                        wrap_plink_parser_onep.writePlinkBedFile2int8FAAA(
                            bedfile,
                            iid_count,
                            sid_count,
                            count_A1,
                            val,
                        )
                    else:
                        wrap_plink_parser_onep.writePlinkBedFile2int8CAAA(
                            bedfile,
                            iid_count,
                            sid_count,
                            count_A1,
                            val,
                        )
                else:
                    raise ValueError(
                        f"dtype '{val.dtype}' not known, only 'int8', 'float32', and 'float64' are allowed."
                    )
            except SystemError as system_error:
                try:
                    bedfile.unlink()
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
                            "Writing snp # {0} to file '{1}'".format(
                                sid_index, filepath
                            )
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
                            elif (val.dtype == np.int8 and val_for_byte == -127) or np.isnan(
                                val_for_byte
                            ): 
                                code = 0b01  # backwards on purpose
                            else:
                                raise ValueError(
                                    "Attempt to write illegal value to BED file. Only 0,1,2,missing allowed."
                                )
                            byte |= code << (val_index * 2)
                        bed_filepointer.write(bytes(bytearray([byte])))
        logging.info(f"Done writing {filepath}")

    def _get_num_threads(self):
        if platform.system() == "Darwin":
            return 1
        if self._num_threads is not None:
            return self._num_threads
        if "MKL_NUM_THREADS" in os.environ:
            return int(os.environ["MKL_NUM_THREADS"])
        return multiprocessing.cpu_count()

    @staticmethod
    def _array_properties_are_ok(val, order, dtype):
        dtype = np.dtype(dtype)

        if val.dtype != dtype:
            return False
        if order == "F":
            return val.flags["F_CONTIGUOUS"]
        elif order == "C":
            return val.flags["C_CONTIGUOUS"]

        return True

    @property
    def shape(self):
        return (len(self.iid), len(self.sid))

    def read(
        self,
        index: Optional[Any] = None,
        dtype: Optional[Union[type, str]] = np.float32,
        order: Optional[str] = "F",
        force_python_only: bool = False,
    ) -> np.ndarray:
        '''
        !!!cmk talk about default dtype and missing and the various index methods
        '''

        iid_index_or_slice_etc, sid_index_or_slice_etc = self._split_index(index)

        dtype = np.dtype(dtype)
        if order == "A":
            order = "F"
        if order not in {"F", "C"}:
            raise ValueError(f"order '{order}' not known, only 'F', 'C', and 'A")

        # Later happy with _iid_range and _sid_range or could it be done with allocation them?
        if self._iid_range is None:
            self._iid_range = np.arange(self.iid_count)
        if self._sid_range is None:
            self._sid_range = np.arange(self.sid_count)

        iid_index = self._iid_range[iid_index_or_slice_etc]
        sid_index = self._sid_range[sid_index_or_slice_etc]

        if not force_python_only:
            num_threads = self._get_num_threads()
            if num_threads>1:
                open_bed._find_openmp()
                from bed_reader import wrap_plink_parser_openmp as wrap_plink_parser
            else:
                from bed_reader import wrap_plink_parser_onep as wrap_plink_parser

            val = np.zeros((len(iid_index), len(sid_index)), order=order, dtype=dtype)
            bed_file_ascii = str(open_bed._name_of_other_file(self.filepath, "bed", "bed")).encode("ascii")


            if self.iid_count > 0 and self.sid_count > 0:
                if dtype == np.int8:
                    if order == "F":
                        wrap_plink_parser.readPlinkBedFile2int8FAAA(
                            bed_file_ascii,
                            self.iid_count,
                            self.sid_count,
                            self.count_A1,
                            iid_index,
                            sid_index,
                            val,
                            num_threads,
                        )
                    elif order == "C":
                        wrap_plink_parser.readPlinkBedFile2int8CAAA(
                            bed_file_ascii,
                            self.iid_count,
                            self.sid_count,
                            self.count_A1,
                            iid_index,
                            sid_index,
                            val,
                            num_threads,
                        )
                    else:
                        assert False, "real assert"
                elif dtype == np.float64:
                    if order == "F":
                        wrap_plink_parser.readPlinkBedFile2doubleFAAA(
                            bed_file_ascii,
                            self.iid_count,
                            self.sid_count,
                            self.count_A1,
                            iid_index,
                            sid_index,
                            val,
                            num_threads,
                        )
                    elif order == "C":
                        wrap_plink_parser.readPlinkBedFile2doubleCAAA(
                            bed_file_ascii,
                            self.iid_count,
                            self.sid_count,
                            self.count_A1,
                            iid_index,
                            sid_index,
                            val,
                            num_threads,
                        )
                    else:
                        assert False, "real assert"
                elif dtype == np.float32:
                    if order == "F":
                        wrap_plink_parser.readPlinkBedFile2floatFAAA(
                            bed_file_ascii,
                            self.iid_count,
                            self.sid_count,
                            self.count_A1,
                            iid_index,
                            sid_index,
                            val,
                            num_threads,
                        )
                    elif order == "C":
                        wrap_plink_parser.readPlinkBedFile2floatCAAA(
                            bed_file_ascii,
                            self.iid_count,
                            self.sid_count,
                            self.count_A1,
                            iid_index,
                            sid_index,
                            val,
                            num_threads,
                        )
                    else:
                        assert False, "real assert"

                else:
                    raise ValueError(
                        f"dtype '{val.dtype}' not known, only 'int8', 'float32', and 'float64' are allowed."
                    )

        else:
            if not self.count_A1:
                byteZero = 0
                byteThree = 2
            else:
                byteZero = 2
                byteThree = 0
            if dtype == np.int8:
                missing = -127
            else:
                missing = np.nan
            # An earlier version of this code had a way to read consecutive SNPs of code in one read. May want
            # to add that ability back to the code.
            # Also, note that reading with python will often result in non-contiguous memory
            # logging.warn("using pure python plink parser (might be much slower!!)")
            val = np.zeros(
                ((int(np.ceil(0.25 * self.iid_count)) * 4), len(sid_index)),
                order=order,
                dtype=dtype,
            )  # allocate it a little big

            bedfile = self._name_of_other_file(self.filepath, "bed", "bed")
            with open(bedfile, "rb") as filepointer:
                for SNPsIndex, bimIndex in enumerate(sid_index):

                    startbit = int(np.ceil(0.25 * self.iid_count) * bimIndex + 3)
                    filepointer.seek(startbit)
                    nbyte = int(np.ceil(0.25 * self.iid_count))
                    bytes = np.array(bytearray(filepointer.read(nbyte))).reshape(
                        (int(np.ceil(0.25 * self.iid_count)), 1), order="F"
                    )

                    val[3::4, SNPsIndex : SNPsIndex + 1] = byteZero
                    val[3::4, SNPsIndex : SNPsIndex + 1][bytes >= 64] = missing
                    val[3::4, SNPsIndex : SNPsIndex + 1][bytes >= 128] = 1
                    val[3::4, SNPsIndex : SNPsIndex + 1][bytes >= 192] = byteThree
                    bytes = np.mod(bytes, 64)
                    val[2::4, SNPsIndex : SNPsIndex + 1] = byteZero
                    val[2::4, SNPsIndex : SNPsIndex + 1][bytes >= 16] = missing
                    val[2::4, SNPsIndex : SNPsIndex + 1][bytes >= 32] = 1
                    val[2::4, SNPsIndex : SNPsIndex + 1][bytes >= 48] = byteThree
                    bytes = np.mod(bytes, 16)
                    val[1::4, SNPsIndex : SNPsIndex + 1] = byteZero
                    val[1::4, SNPsIndex : SNPsIndex + 1][bytes >= 4] = missing
                    val[1::4, SNPsIndex : SNPsIndex + 1][bytes >= 8] = 1
                    val[1::4, SNPsIndex : SNPsIndex + 1][bytes >= 12] = byteThree
                    bytes = np.mod(bytes, 4)
                    val[0::4, SNPsIndex : SNPsIndex + 1] = byteZero
                    val[0::4, SNPsIndex : SNPsIndex + 1][bytes >= 1] = missing
                    val[0::4, SNPsIndex : SNPsIndex + 1][bytes >= 2] = 1
                    val[0::4, SNPsIndex : SNPsIndex + 1][bytes >= 3] = byteThree
                val = val[iid_index, :]  # reorder or trim any extra allocation
                if not open_bed._array_properties_are_ok(val, order, dtype):
                    val = val.copy(order=order)

        return val

    @staticmethod
    def _split_index(index):
        if not isinstance(index, tuple):
            index = (None, index)
        iid_index = open_bed._fix_up_index(index[0])
        sid_index = open_bed._fix_up_index(index[1])
        return iid_index, sid_index

    @staticmethod
    def _fix_up_index(index):
        if index is None:  # make a shortcut for None
            return slice(None)
        try:  # If index is an int, return it in an array
            index = index.__index__()  # (see
            # https://stackoverflow.com/questions/3501382/checking-whether-a-variable-is-an-integer-or-not)
            return [index]
        except Exception:
            pass
        return index

    @staticmethod
    def _write_fam_or_bim(basefilepath, metadata, suffix_of_interest):
        assert suffix_of_interest in {"fam", "bim"}, "real assert"

        filepath = open_bed._name_of_other_file(basefilepath, "bed", suffix_of_interest)

        fam_bim_list = []
        for key, mm in _meta_meta.items():
            if mm.suffix == suffix_of_interest:
                assert len(fam_bim_list) == mm.column, "real assert"
                fam_bim_list.append(metadata[key])

        sep = " " if suffix_of_interest == "fam" else "\t"

        with open(filepath, "w") as filepointer:
            for index in range(len(fam_bim_list[0])):
                filepointer.write(
                    sep.join(str(seq[index]) for seq in fam_bim_list) + "\n"
                )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    import os

    if True: #!!!cmk
        import numpy as np
        from bed_reader._open_bed import open_bed

        # Can get file from https://www.dropbox.com/sh/xluk9opjiaobteg/AABgEggLk0ZoO0KQq0I4CaTJa?dl=0
        bigfile = r"M:\deldir\genbgen\2\merged_487400x220000.1.bed"
        # bigfile = '/mnt/m/deldir/genbgen/2/merged_487400x220000.1.bed'
        with open_bed(bigfile, num_threads=20) as bed:
            sid_batch = 22 * 1000
            for sid_start in range(0, 10 * sid_batch, sid_batch):
                slicer = np.s_[:10000, sid_start : sid_start + sid_batch]
                print(slicer)
                val = bed.read(slicer)
            print(val.shape)

    if False:
        file = r"D:\OneDrive\programs\sgkit-plink\bed_reader\tests\data/plink_sim_10s_100v_10pmiss.bed"
        with open_bed(file) as bed:
            print(bed.iid)
            print(bed.shape)
            val = bed.read()
            print(val)

    if False:

        # bed_file = example_file('doc/ipynb/all.*','*.bed')
        bed_file = r"F:\backup\carlk4d\data\carlk\cachebio\genetics\onemil\id1000000.sid_1000000.seed0.byiid\iid990000to1000000.bed"
        bed = Bed(bed_file, count_A1=False)
        snpdata1 = bed[:, :1000].read()
        snpdata2 = bed[:, :1000].read(dtype="int8", _require_float32_64=False)
        print(snpdata2)
        snpdata3 = bed[:, :1000].read(
            dtype="int8", order="C", _require_float32_64=False
        )
        print(snpdata3)
        snpdata3.val = snpdata3.val.astype("float32")
        snpdata3.val.dtype

    if False:
        from bed_reader import Bed, SnpGen

        iid_count = 487409
        sid_count = 5000
        sid_count_max = 5765294
        sid_batch_size = 50

        sid_batch_count = -(sid_count // -sid_batch_size)
        sid_batch_count_max = -(sid_count_max // -sid_batch_size)
        snpgen = SnpGen(seed=234, iid_count=iid_count, sid_count=sid_count_max)

        for batch_index in range(sid_batch_count):
            sid_index_start = batch_index * sid_batch_size
            sid_index_end = (batch_index + 1) * sid_batch_size  # what about rounding
            filename = r"d:\deldir\rand\fakeukC{0}x{1}-{2}.bed".format(
                iid_count, sid_index_start, sid_index_end
            )
            if not os.path.exists(filename):
                Bed.write(
                    filename + ".temp", snpgen[:, sid_index_start:sid_index_end].read()
                )
                os.rename(filename + ".temp", filename)

    if False:
        from bed_reader import Pheno, Bed

        filename = r"m:\deldir\New folder (4)\all_chr.maf0.001.N300.bed"
        iid_count = 300
        iid = [["0", "iid_{0}".format(iid_index)] for iid_index in range(iid_count)]
        bed = Bed(filename, iid=iid, count_A1=False)
        print(bed.iid_count)

    if False:
        from pysnptools.util import example_file

        pheno_fn = example_file("pysnptools/examples/toydata.phe")

    if False:
        from bed_reader import Pheno, Bed

        print(os.getcwd())
        snpdata = Pheno("../examples/toydata.phe").read()  # Read data from Pheno format
        # pstutil.create_directory_if_necessary("tempdir/toydata.5chrom.bed")
        Bed.write(
            "tempdir/toydata.5chrom.bed", snpdata, count_A1=False
        )  # Write data in Bed format

    import doctest

    #!!!cmk put this back
    # doctest.testmod(
    #    optionflags=doctest.ELLIPSIS
    # )  #!!!cmk how do you doctest with PyTest?
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
