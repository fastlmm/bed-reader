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
import logging

import math
from typing import Any, List, Optional, Union, Mapping

from itertools import takewhile, repeat

# https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python
def _rawincount(filepath):
    with open(filepath, "rb") as f:
        bufgen = takewhile(lambda x: x, (f.raw.read(1024 * 1024) for _ in repeat(None)))
        return sum(buf.count(b"\n") for buf in bufgen)


@dataclass
class _MetaMeta:
    suffix: str
    column: int
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
    A NumPy-inspired class for fast opening and reading of PLINK cmkstar.bed files.

    Parameters
    ----------
    filepath:
        Bed file path #!!!cmk Bed vs BED vs cmkstar.bed
    iid_count: None or int, optional
        Number of individuals (samples) in the BED file.
        The default (``iid_count=None``) finds the number
        automatically by quickly scanning the FAM file.
    sid_count: None or int, optional
        Number of SNPs (variants) in the BED file.
        The default (``sid_count=None``) finds the number
        automatically by quickly scanning the BIM file.
    metadata: dict, optional
        A dictionary of any replacement metadata. The default is an empty dictionary.
        The keys of the dictionary are the names of the metadata to replace.
        The possible keys are:

             "fid" (family id), "iid" (individual or sample id), "father" (father id),
             "mother" (mother id), "sex", "pheno" (phenotype), "chromosome", "sid"
             (SNP or variant id), "cm_position" (centimorgan position), "bp_position"
             (base-pair position), "allele_1", "allele_2".
            
        The values are the replacement lists or arrays. CMK see example
    count_A1: bool, optional
        True (default) to count the number of A1 alleles (the PLINK standard). False to count the number of A2 alleles.
    num_threads: None or int, optional
        The number of threads with which to read data. Defaults to all available threads.
        Can also be set with the 'MKL_NUM_THREADS' environment variable.
        
        On MacOS, this
        parameter is ignored and all reads are singled threaded. On Windows, if reads create
        library-loading problems, setting this to 1 will fix the problems (but be slower).
    skip_format_check: bool, optional
        False (default) to immediately check for expected starting bytes in the BED file.
        True to delay the check until (and if) data is read.
        
    Returns
    -------
    an open_bed object : :class:`open_bed`

    .._open_examples:

    Examples
    --------

    With the `with <https://docs.python.org/3/reference/compound_stmts.html#grammar-token-with-stmt>`__ statement,
    list individual (sample) :attr:`iid` and SNP (variant) :attr:`sid`, then :meth:`read` the whole file.

    .. doctest::

        >>> from bed_reader import open_bed, sample_file
        >>>
        >>> file_name = sample_file("small.bed")
        >>> with open_bed(file_name) as bed:
        ...     print(bed.iid)
        ...     print(bed.sid)
        ...     print(bed.read())
        ['iid1' 'iid2' 'iid3']
        ['sid1' 'sid2' 'sid3' 'sid4']
        [[ 1.  0. nan  0.]
         [ 2.  0. nan  2.]
         [ 0.  1.  2.  0.]]

    Open the file (without `with`) and read probabilities for one SNP (variant)
    at index position 2.

    .. doctest::

        >>> bed = open_bed(file_name)
        >>> print(bed.read(index=2))
        [[nan]
         [nan]
         [ 2.]]
        >>> del bed  # optional: close and delete object

    Replace the sample :attr:`iid`.

        >>> with open_bed(file_name, metadata={"iid":["sample1","sample2","sample3"]}) as bed:
        ...     print(bed.iid)
        ...     print(bed.sid)
        ['sample1' 'sample2' 'sample3']
        ['sid1' 'sid2' 'sid3' 'sid4']

    Tell it the number of individuals (samples) and SNPs (variants). This lets it read data without
    the needing to ever open the cmkstar.fam and cmkstar.bim files.

        >>> with open_bed(file_name, iid_count=3, sid_count=4) as bed:
        ...     print(bed.read())
        [[ 1.  0. nan  0.]
         [ 2.  0. nan  2.]
         [ 0.  1.  2.  0.]]


    See the :meth:`read` for details of reading batches via slicing and fancy indexing.

        #!!!cmk in README say: Documentation not API Documatnion

    .. _sample format: https://www.well.ox.ac.uk/~gav/qctool/documentation/sample_file_formats.html #!!!cmk

    """

    def __init__(
        self,
        filepath: Union[str, Path],
        iid_count: Optional[int] = None,
        sid_count: Optional[int] = None,
        metadata: Mapping[str, List[Any]] = {},
        count_A1: bool = True,
        num_threads: Optional[int] = None,
        skip_format_check: bool = False,
    ):  #!!!document these new optionals. they are here
        self.filepath = Path(filepath)
        self.count_A1 = count_A1
        self._num_threads = num_threads
        self.skip_format_check = skip_format_check

        self.metadata_dict, self._counts = open_bed._fix_up_metadata(
            metadata, iid_count, sid_count, use_fill_sequence=False
        )
        self._iid_range = None
        self._sid_range = None

        if not self.skip_format_check:
            bedfile = open_bed._name_of_other_file(self.filepath, "bed", "bed")
            with open(bedfile, "rb") as filepointer:
                self._check_file(filepointer)

    def read(
        self,
        index: Optional[Any] = None,
        dtype: Optional[Union[type, str]] = np.float32,
        order: Optional[str] = "F",
        force_python_only: Optional[bool] = False,
    ) -> np.ndarray:
        """
        Read genotype information from an :class:`open_bed` object.

        Parameters
        ----------
        index:
            An optional expression specifying the individuals (samples) and SNPs (variants)
            to read. (See :ref:`read_examples`, below).
            Defaults to ``None``, meaning read all.
        dtype: {'float32' (default), 'float64', 'int8'}, optional
            The desired data-type for the returned array.
        order : {'F','C'}, optional
            The desired memory layout for the returned array.
            Defaults to ``F`` (Fortran order, which is SNP-major).
        force_python_only: bool, optional
            If False (default), uses the faster C++ code; otherwise it uses the slower pure Python code.

        Returns
        -------
        :class:`numpy.ndarray` of 0, 1, 2, or missing values with ``dtype`` and shape `(iid_count,sid_count)`

        .. _read_examples:

        Examples
        --------
        * Index Examples

        To read all data in a BED file, set ``index`` to ``None``. This is the default.

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.read())
            [[ 1.  0. nan  0.]
             [ 2.  0. nan  2.]
             [ 0.  1.  2.  0.]]

        To read selected SNPs (variants), set ``index`` to an ``int``, a list of ``int``, a :class:`slice`, or a list of ``bool``.
        Negative integers count from the end of the data.


        .. doctest::

            >>> bed = open_bed(file_name)
            >>> print(bed.read(2))  # read the SNPs indexed by 2.
            [[nan]
             [nan]
             [ 2.]]
            >>> print(bed.read([2,3,0]))  # read the SNPs indexed by 2, 3, and 0
            [[nan  0.  1.]
             [nan  2.  2.]
             [ 2.  0.  0.]]
            >>> print(bed.read(slice(2))) #read the first 2 SNPs
            [[1. 0.]
             [2. 0.]
             [0. 1.]]
            >>> print(bed.read(slice(1,4))) #read SNPs from 1 (inclusive) to 4 (exclusive)
            [[ 0. nan  0.]
             [ 0. nan  2.]
             [ 1.  2.  0.]]
            >>> print(bed.read(slice(2,None))) # read SNPs starting at index 2.
            [[nan  0.]
             [nan  2.]
             [ 2.  0.]]
            >>> print(bed.read(slice(None,None,2))) #read every 2nd SNPs
            [[ 1. nan]
             [ 2. nan]
             [ 0.  2.]]
            >>> print(np.unique(bed.chromosome)) # print unique chrom values
            ['1' '5' 'Y']
            >>> print(bed.read(bed.chromosome=='5')) # read all SNPs in chrom 1
            [[nan]
             [nan]
             [ 2.]]
            >>> print(bed.read(-1)) # read the last SNPs
            [[0.]
             [2.]
             [0.]]


        To read selected individuals (samples), set ``index`` to a tuple of the form ``(individual_index,None)``, where ``individual_index`` follows the form
        of ``SNP index``, above.

        .. doctest::

            >>> print(bed.read((0,None))) # Read 1st individual (across all SNPs)
            [[ 1.  0. nan  0.]]
            >>> print(bed.read((slice(None,None,2),None))) # Read every 2nd individual
            [[ 1.  0. nan  0.]
             [ 0.  1.  2.  0.]]


        To read selected individuals and selected SNPs, set ``index`` to a tuple of the form ``(individual_index,SNP_index)``,
        where ``individual_index`` and ``SNP_index`` follow the forms above.

        .. doctest::

            >>> # Read individuals 1 (inclusive) to 3 (exclusive) and the first 2 SNPs.
            >>> print(bed.read((slice(1,3),slice(2))))
            [[2. 0.]
             [0. 1.]]
            >>> #read last and 2nd-to-last individuals and the last SNPs
            >>> print(bed.read(([-1,-2],-1)))
            [[0.]
             [2.]]


        * dtype example

        You can give a dtype. For float32 and float64, NaN indicates missing values.
        For int8, -127 indicates missing values.

        .. doctest::

            >>> print(bed.read(dtype='int8'))
            [[   1    0 -127    0]
             [   2    0 -127    2]
             [   0    1    2    0]]
            >>> del bed  # optional: close and delete object

        """

        iid_index_or_slice_etc, sid_index_or_slice_etc = self._split_index(index)

        dtype = np.dtype(dtype)
        if order not in {"F", "C"}:
            raise ValueError(f"order '{order}' not known, only 'F', 'C'")

        # Later happy with _iid_range and _sid_range or could it be done with allocation them?
        if self._iid_range is None:
            self._iid_range = np.arange(self.iid_count)
        if self._sid_range is None:
            self._sid_range = np.arange(self.sid_count)

        iid_index = self._iid_range[iid_index_or_slice_etc]
        sid_index = self._sid_range[sid_index_or_slice_etc]

        if not force_python_only:
            num_threads = self._get_num_threads()
            if num_threads > 1:
                open_bed._find_openmp()
                from bed_reader import wrap_plink_parser_openmp as wrap_plink_parser
            else:
                from bed_reader import wrap_plink_parser_onep as wrap_plink_parser

            val = np.zeros((len(iid_index), len(sid_index)), order=order, dtype=dtype)
            bed_file_ascii = str(
                open_bed._name_of_other_file(self.filepath, "bed", "bed")
            ).encode("ascii")

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

    def __str__(self):
        return f"{self.__class__.__name__}('{self.filepath}',...)"

    #!!!cmk make sure these are in a good order for the documentation

    @property
    def fid(self):
        """
        The family id (a :class:`numpy.ndarray` of ``str``).

        If needed, will cause a one-time read of the cmkstar.fam file.

        Example
        -------
        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.fid)
            ['fid1' 'fid1' 'fid2']

        """

        return self.metadata_item("fid")

    @property
    def iid(self):
        """
        The individual id (a :class:`numpy.ndarray` of ``str``).

        If needed, will cause a one-time read of the cmkstar.fam file.

        Example
        -------
        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.iid)
            ['iid1' 'iid2' 'iid3']

        """
        return self.metadata_item("iid")

    @property
    def father(self):
        """
        The father id
       
        :rtype:  :class:`numpy.ndarray` of ``str``
        
        If needed, will cause a one-time read of the cmkstar.fam file.

        Example
        -------
        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.father)
            ['iid23' 'iid23' 'iid22']

        """
        return self.metadata_item("father")

    @property
    def mother(self):
        """
        !!!cmk need doc string
        """
        return self.metadata_item("mother")

    @property
    def sex(self):
        """
        !!!cmk need doc string
        """
        return self.metadata_item("sex")

    @property
    def pheno(self):
        """
        !!!cmk need doc string
        """
        return self.metadata_item("pheno")

    @property
    def metadata(self):
        """
        !!!cmk need doc string
        !!!cmk tell that if needed, will open and read cmkstar.fam and cmkstar.bim files
        """
        for key in _meta_meta:
            self.metadata_item(key)
        return self.metadata_dict

    def metadata_item(self, key):
        """
        !!!cmk need doc string
        """
        val = self.metadata_dict.get(key)
        if val is None:
            mm = _meta_meta[key]
            self._read_fam_or_bim(suffix=mm.suffix)
            return self.metadata_dict[key]
        else:
            return val

    @property
    def chromosome(self):
        """
        !!!cmk need doc string
        """
        return self.metadata_item("chromosome")

    @property
    def sid(self):
        """
        !!!cmk need doc string
        """
        return self.metadata_item("sid")

    @property
    def cm_position(self):
        """
        !!!cmk need doc string
        """
        return self.metadata_item("cm_position")

    @property
    def bp_position(self):
        """
        !!!cmk need doc string
        """
        return self.metadata_item("bp_position")

    @property
    def allele_1(self):
        """
        !!!cmk need doc string
        """
        return self.metadata_item("allele_1")

    @property
    def allele_2(self):
        """
        !!!cmk need doc string
        """
        return self.metadata_item("allele_2")

    @property
    def iid_count(self):
        """
        !!!cmk need doc string
        """
        return self._count("fam")

    @property
    def sid_count(self):
        """
        !!!cmk need doc string
        """
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
            raise ValueError("No valid binary BED file")
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
    def _get_version_number(filename):
        # http://timgolden.me.uk/python/win32_how_do_i/get_dll_version.html
        from win32api import GetFileVersionInfo, LOWORD, HIWORD

        info = GetFileVersionInfo(filename, "\\")
        ms = info["FileVersionMS"]
        ls = info["FileVersionLS"]
        return HIWORD(ms), LOWORD(ms), HIWORD(ls), LOWORD(ls)

    @staticmethod
    def _find_openmp():
        if "bed_reader.wrap_plink_parser_openmp" in sys.modules:
            return
        if platform.system() == "Windows":
            logging.info("in windows _find_openmp")
            from ctypes import cdll
            from ctypes.util import find_library

            dllname = "libomp.dll"
            find_location = find_library(dllname)
            if find_location is not None:  #!!!cmk
                logging.info(f"found '{dllname}' at '{find_library(dllname)}'")
                found_ver = open_bed._get_version_number(find_location)
                goal_ver = (5, 0, 2014, 926)
                logging.info(f"found ver is '{found_ver}'. Goal ver is '{goal_ver}'")
                if found_ver >= goal_ver:
                    logging.info("found version looks good, so load that")
                    cdll.LoadLibrary(str(find_location))
                    return
            location_list = [
                Path(__file__).parent / dllname,
                Path(__file__).parent.parent / "external/llvm/windows/bin" / dllname,
            ]
            for location in location_list:
                if location.exists():
                    logging.info(f"loading my own version from '{location}'")
                    cdll.LoadLibrary(str(location))
                    return
            raise Exception(f"Can't find '{dllname}'")

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

    @staticmethod
    def _fix_up_metadata_array(input, dtype, missing_value, key):
        if len(input) == 0:
            return np.zeros([0], dtype=dtype)

        if not isinstance(input, np.ndarray):
            return open_bed._fix_up_metadata_array(
                np.array(input), dtype, missing_value, key
            )

        if len(input.shape) != 1:
            raise ValueError(f"{key} should be one dimensional")

        if not np.issubdtype(input.dtype, dtype):
            output = np.array(input, dtype=dtype)
            # If converting float to non-float: Change NaN to new missing value
            if np.isrealobj(input) and not np.isrealobj(output):
                output[input != input] = missing_value
            return output

        return input

    @staticmethod
    def _fix_up_metadata(metadata, iid_count, sid_count, use_fill_sequence):

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
            else:
                output = open_bed._fix_up_metadata_array(
                    input, mm.dtype, mm.missing_value, key
                )

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
                    output = np.array(
                        fields[mm.column].fillna(mm.missing_value), dtype=mm.dtype
                    )
                self.metadata_dict[key] = output

    @staticmethod
    def _name_of_other_file(filepath, remove_suffix, add_suffix):
        if filepath.suffix.lower() == "." + remove_suffix:
            filepath = filepath.parent / filepath.stem
        return filepath.parent / (filepath.name + "." + add_suffix)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    import os

    if True:  #!!!cmk
        from bed_reader import open_bed, sample_file

        file_name = sample_file("small.bed")
        with open_bed(file_name) as bed:
            print(bed.iid)
            print(bed.sid)
            print(bed.read())

    # if False:  #!!!cmk
    #     import numpy as np
    #     from bed_reader._open_bed import open_bed

    #     # Can get file from https://www.dropbox.com/sh/xluk9opjiaobteg/AABgEggLk0ZoO0KQq0I4CaTJa?dl=0
    #     bigfile = r"M:\deldir\genbgen\2\merged_487400x220000.1.bed"
    #     # bigfile = '/mnt/m/deldir/genbgen/2/merged_487400x220000.1.bed'
    #     with open_bed(bigfile, num_threads=20) as bed:
    #         sid_batch = 22 * 1000
    #         for sid_start in range(0, 10 * sid_batch, sid_batch):
    #             slicer = np.s_[:10000, sid_start : sid_start + sid_batch]
    #             print(slicer)
    #             val = bed.read(slicer)
    #         print(val.shape)

    # if False:
    #     file = r"D:\OneDrive\programs\sgkit-plink\bed_reader\tests\data/plink_sim_10s_100v_10pmiss.bed"
    #     with open_bed(file) as bed:
    #         print(bed.iid)
    #         print(bed.shape)
    #         val = bed.read()
    #         print(val)

    # if False:

    #     # bed_file = example_file('doc/ipynb/all.*','*.bed')
    #     bed_file = r"F:\backup\carlk4d\data\carlk\cachebio\genetics\onemil\id1000000.sid_1000000.seed0.byiid\iid990000to1000000.bed"
    #     bed = Bed(bed_file, count_A1=False)
    #     snpdata1 = bed[:, :1000].read()
    #     snpdata2 = bed[:, :1000].read(dtype="int8", _require_float32_64=False)
    #     print(snpdata2)
    #     snpdata3 = bed[:, :1000].read(
    #         dtype="int8", order="C", _require_float32_64=False
    #     )
    #     print(snpdata3)
    #     snpdata3.val = snpdata3.val.astype("float32")
    #     snpdata3.val.dtype

    # if False:
    #     from bed_reader import Bed, SnpGen

    #     iid_count = 487409
    #     sid_count = 5000
    #     sid_count_max = 5765294
    #     sid_batch_size = 50

    #     sid_batch_count = -(sid_count // -sid_batch_size)
    #     sid_batch_count_max = -(sid_count_max // -sid_batch_size)
    #     snpgen = SnpGen(seed=234, iid_count=iid_count, sid_count=sid_count_max)

    #     for batch_index in range(sid_batch_count):
    #         sid_index_start = batch_index * sid_batch_size
    #         sid_index_end = (batch_index + 1) * sid_batch_size  # what about rounding
    #         filename = r"d:\deldir\rand\fakeukC{0}x{1}-{2}.bed".format(
    #             iid_count, sid_index_start, sid_index_end
    #         )
    #         if not os.path.exists(filename):
    #             Bed.write(
    #                 filename + ".temp", snpgen[:, sid_index_start:sid_index_end].read()
    #             )
    #             os.rename(filename + ".temp", filename)

    # if False:
    #     from bed_reader import Pheno, Bed

    #     filename = r"m:\deldir\New folder (4)\all_chr.maf0.001.N300.bed"
    #     iid_count = 300
    #     iid = [["0", "iid_{0}".format(iid_index)] for iid_index in range(iid_count)]
    #     bed = Bed(filename, iid=iid, count_A1=False)
    #     print(bed.iid_count)

    # if False:
    #     from pysnptools.util import example_file

    #     pheno_fn = example_file("pysnptools/examples/toydata.phe")

    # if False:
    #     from bed_reader import Pheno, Bed

    #     print(os.getcwd())
    #     snpdata = Pheno("../examples/toydata.phe").read()  # Read data from Pheno format
    #     # pstutil.create_directory_if_necessary("tempdir/toydata.5chrom.bed")
    #     Bed.write(
    #         "tempdir/toydata.5chrom.bed", snpdata, count_A1=False
    #     )  # Write data in Bed format

    import pytest
    pytest.main(["--doctest-modules", __file__])