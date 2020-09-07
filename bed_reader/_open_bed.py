# LATER: Offer to ignore some or all fam bim fields
import logging
import math
import multiprocessing
import numbers
import os
import platform
import sys
from dataclasses import dataclass
from itertools import repeat, takewhile
from pathlib import Path
from typing import Any, List, Mapping, Optional, Union

import numpy as np
import pandas as pd


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


class open_bed:
    """
    Open a PLINK .bed file for reading.

    Parameters
    ----------
    filepath: pathlib.Path or str
        File path to the .bed file.
    iid_count: None or int, optional
        Number of individuals (samples) in the .bed file.
        The default (``iid_count=None``) finds the number
        automatically by quickly scanning the .fam file.
    sid_count: None or int, optional
        Number of SNPs (variants) in the .bed file.
        The default (``sid_count=None``) finds the number
        automatically by quickly scanning the .bim file.
    properties: dict, optional
        A dictionary of any replacement properties. The default is an empty dictionary.
        The keys of the dictionary are the names of the properties to replace.
        The possible keys are:

             "fid" (family id), "iid" (individual or sample id), "father" (father id),
             "mother" (mother id), "sex", "pheno" (phenotype), "chromosome", "sid"
             (SNP or variant id), "cm_position" (centimorgan position), "bp_position"
             (base-pair position), "allele_1", "allele_2".            

          The values are replacement lists or arrays. A value can also be `None`,
          meaning do not read or offer this property. See examples, below.

          The list or array will be converted to a :class:`numpy.ndarray`
          of the appropriate dtype, if necessary. Any :class:`numpy.nan` values will converted to
          the appropriate missing value. The PLINK `.fam specification <https://www.cog-genomics.org/plink2/formats#fam>`_
          and `.bim specification <https://www.cog-genomics.org/plink2/formats#bim>`_ lists
          the dtypes and missing values for each property.

    count_A1: bool, optional
        True (default) to count the number of A1 alleles (the PLINK standard). False to count the number of A2 alleles.
    num_threads: None or int, optional
        The number of threads with which to read data. Defaults to all available processors.
        Can also be set with the 'MKL_NUM_THREADS' environment variable.
        
        On MacOS, this parameter is ignored and all reads are singled threaded.
    skip_format_check: bool, optional
        False (default) to immediately check for expected starting bytes in the .bed file.
        True to delay the check until (and if) data is read.
    fam_filepath: pathlib.Path or str, optional
        Path to the file containing information about each individual (sample).
        Defaults to replacing the .bed file’s suffix with .fam.
    bim_filepath: pathlib.Path or str, optional
        Path to the file containing information about each SNP (variant).
        Defaults to replacing the .bed file’s suffix with .bim.
        
    Returns
    -------
    open_bed
        an open_bed object

    Examples
    --------

    List individual (sample) :attr:`iid` and SNP (variant) :attr:`sid`, then :meth:`read` the whole file.

    .. doctest::

        >>> from bed_reader import open_bed, sample_file
        >>>
        >>> file_name = sample_file("small.bed")
        >>> bed = open_bed(file_name)
        >>> print(bed.iid)
        ['iid1' 'iid2' 'iid3']
        >>> print(bed.sid)
        ['sid1' 'sid2' 'sid3' 'sid4']
        >>> print(bed.read())
        [[ 1.  0. nan  0.]
         [ 2.  0. nan  2.]
         [ 0.  1.  2.  0.]]
        >>> del bed  # optional: delete bed object

    Open the file and read data for one SNP (variant)
    at index position 2.

    .. doctest::

        >>> import numpy as np
        >>> with open_bed(file_name) as bed:
        ...     print(bed.read(np.s_[:,2]))
        [[nan]
         [nan]
         [ 2.]]

    Replace :attr:`iid`.

        >>> bed = open_bed(file_name, properties={"iid":["sample1","sample2","sample3"]})
        >>> print(bed.iid) # replaced
        ['sample1' 'sample2' 'sample3']
        >>> print(bed.sid) # same as before
        ['sid1' 'sid2' 'sid3' 'sid4']

    Give the number of individuals (samples) and SNPs (variants) so that the .fam and .bim files need never be opened.

        >>> with open_bed(file_name, iid_count=3, sid_count=4) as bed:
        ...     print(bed.read())
        [[ 1.  0. nan  0.]
         [ 2.  0. nan  2.]
         [ 0.  1.  2.  0.]]

    Mark some properties as "don’t read or offer".

        >>> bed = open_bed(file_name, properties={
        ...    "father" : None, "mother" : None, "sex" : None, "pheno" : None,
        ...    "allele_1" : None, "allele_2":None })
        >>> print(bed.iid)        # read from file
        ['iid1' 'iid2' 'iid3']
        >>> print(bed.allele_2)   # not read and not offered
        None

    See the :meth:`read` for details of reading batches via slicing and fancy indexing.

    """

    def __init__(
        self,
        filepath: Union[str, Path],
        iid_count: Optional[int] = None,
        sid_count: Optional[int] = None,
        properties: Mapping[str, List[Any]] = {},
        count_A1: bool = True,
        num_threads: Optional[int] = None,
        skip_format_check: bool = False,
        fam_filepath: Union[str, Path] = None,
        bim_filepath: Union[str, Path] = None,
    ):  #!!!document these new optionals. they are here
        self.filepath = Path(filepath)
        self.count_A1 = count_A1
        self._num_threads = num_threads
        self.skip_format_check = skip_format_check
        self._fam_filepath = (
            Path(fam_filepath)
            if fam_filepath is not None
            else self.filepath.parent / (self.filepath.stem + ".fam")
        )
        self._bim_filepath = (
            Path(bim_filepath)
            if bim_filepath is not None
            else self.filepath.parent / (self.filepath.stem + ".bim")
        )

        self.properties_dict, self._counts = open_bed._fix_up_properties(
            properties, iid_count, sid_count, use_fill_sequence=False
        )
        self._iid_range = None
        self._sid_range = None

        if not self.skip_format_check:
            with open(self.filepath, "rb") as filepointer:
                self._check_file(filepointer)

    def read(
        self,
        index: Optional[Any] = None,
        dtype: Optional[Union[type, str]] = 'float32',
        order: Optional[str] = "F",
        force_python_only: Optional[bool] = False,
    ) -> np.ndarray:
        """
        Read genotype information.

        Parameters
        ----------
        index:
            An optional expression specifying the individuals (samples) and SNPs (variants)
            to read. (See examples, below).
            Defaults to ``None``, meaning read all.

            (If index is a tuple, the first component indexes the individuals and the second indexes
            the SNPs. If it is not a tuple and not None, it indexes SNPs.)

        dtype: {'float32' (default), 'float64', 'int8'}, optional
            The desired data-type for the returned array.
        order : {'F','C'}, optional
            The desired memory layout for the returned array.
            Defaults to ``F`` (Fortran order, which is SNP-major).
        force_python_only: bool, optional
            If False (default), uses the faster C++ code; otherwise it uses the slower pure Python code.

        Returns
        -------
        numpy.ndarray
            2-D array containing values of 0, 1, 2, or missing


        Rows represent individuals (samples). Columns represent SNPs (variants). 
        
        For ``dtype`` 'float32' and 'float64', NaN indicates missing values.
        For 'int8', -127 indicates missing values.

        Examples
        --------

        To read all data in a .bed file, set ``index`` to ``None``. This is the default.

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.read())
            [[ 1.  0. nan  0.]
             [ 2.  0. nan  2.]
             [ 0.  1.  2.  0.]]

        To read selected individuals (samples) and/or SNPs (variants), set each part of a :class:`numpy.s_` to an `int`, a list of `int`, a slice expression, or a list of `bool`.
        Negative integers count from the end of the list.


        .. doctest::

            >>> import numpy as np
            >>> bed = open_bed(file_name)
            >>> print(bed.read(np.s_[:,2]))  # read the SNPs indexed by 2.
            [[nan]
             [nan]
             [ 2.]]
            >>> print(bed.read(np.s_[:,[2,3,0]]))  # read the SNPs indexed by 2, 3, and 0
            [[nan  0.  1.]
             [nan  2.  2.]
             [ 2.  0.  0.]]
            >>> print(bed.read(np.s_[:,1:4])) #read SNPs from 1 (inclusive) to 4 (exclusive)
            [[ 0. nan  0.]
             [ 0. nan  2.]
             [ 1.  2.  0.]]
            >>> print(np.unique(bed.chromosome)) # print unique chrom values
            ['1' '5' 'Y']
            >>> print(bed.read(np.s_[:,bed.chromosome=='5'])) # read all SNPs in chrom 5
            [[nan]
             [nan]
             [ 2.]]
            >>> print(bed.read(np.s_[0,:])) # Read 1st individual (across all SNPs)
            [[ 1.  0. nan  0.]]
            >>> print(bed.read(np.s_[::2,:])) # Read every 2nd individual
            [[ 1.  0. nan  0.]
             [ 0.  1.  2.  0.]]
            >>> #read last and 2nd-to-last individuals and the last SNPs
            >>> print(bed.read(np.s_[[-1,-2],-1]))
            [[0.]
             [2.]]


        You can give a dtype for the output.

        .. doctest::

            >>> print(bed.read(dtype='int8'))
            [[   1    0 -127    0]
             [   2    0 -127    2]
             [   0    1    2    0]]
            >>> del bed  # optional: delete bed object

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
                from bed_reader import wrap_plink_parser_openmp as wrap_plink_parser
            else:
                from bed_reader import wrap_plink_parser_onep as wrap_plink_parser

            val = np.zeros((len(iid_index), len(sid_index)), order=order, dtype=dtype)
            bed_file_ascii = str(self.filepath).encode("ascii")

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
                    else:
                        assert order == "C" # real assert
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
                    else: 
                        assert order == "C" #real assert
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
                    else:
                        assert order == "C" # real assert
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

            with open(self.filepath, "rb") as filepointer:
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
                assert val.dtype == np.dtype(dtype)  # real assert
                if not open_bed._array_properties_are_ok(val, order):
                    val = val.copy(order=order)

        return val

    def __str__(self) -> str:
        return f"{self.__class__.__name__}('{self.filepath}',...)"

    @property
    def fid(self) -> np.ndarray:
        """
        Family id of each individual (sample).
       
        Returns
        -------
        numpy.ndarray
            array of str


        '0' represents a missing value.
        
        If needed, will cause a one-time read of the .fam file.

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

        return self.property_item("fid")

    @property
    def iid(self) -> np.ndarray:
        """
        Individual id of each individual (sample).
       
        Returns
        -------
        numpy.ndarray
            array of str


        If needed, will cause a one-time read of the .fam file.

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
        return self.property_item("iid")

    @property
    def father(self) -> np.ndarray:
        """
        Father id of each individual (sample).
       
        Returns
        -------
        numpy.ndarray
            array of str


        '0' represents a missing value.
        
        If needed, will cause a one-time read of the .fam file.

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
        return self.property_item("father")

    @property
    def mother(self) -> np.ndarray:
        """
        Mother id of each individual (sample).
       
        Returns
        -------
        numpy.ndarray
            array of str


        '0' represents a missing value.
        
        If needed, will cause a one-time read of the .fam file.

        Example
        -------

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.mother)
            ['iid34' 'iid34' 'iid33']

        """
        return self.property_item("mother")

    @property
    def sex(self) -> np.ndarray:
        """
        Sex of each individual (sample).
       
        Returns
        -------
        numpy.ndarray
            array of 0, 1, or 2


        0 is unknown, 1 is male, 2 is female

        If needed, will cause a one-time read of the .fam file.

        Example
        -------
        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.sex)
            [1 2 0]

        """
        return self.property_item("sex")

    @property
    def pheno(self) -> np.ndarray:
        """
        A phenotype for each individual (sample)
        (seldom used).
       
        Returns
        -------
        numpy.ndarray
            array of str


        '0' may represent a missing value.
        
        If needed, will cause a one-time read of the .fam file.

        Example
        -------

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.pheno)
            ['red' 'red' 'blue']

        """
        return self.property_item("pheno")

    @property
    def properties(self) -> Mapping[str, np.array]:
        """
        All the properties returned as a dictionary.
       
        Returns
        -------
        dict
            all the properties


        The keys of the dictionary are the names of the properties, namely:

             "fid" (family id), "iid" (individual or sample id), "father" (father id),
             "mother" (mother id), "sex", "pheno" (phenotype), "chromosome", "sid"
             (SNP or variant id), "cm_position" (centimorgan position), "bp_position"
             (base-pair position), "allele_1", "allele_2".
            
        The values are :class:`numpy.ndarray`.

        If needed, will cause a one-time read of the .fam and .bim file.

        Example
        -------

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(len(bed.properties)) #length of dict
            12

        """
        for key in _meta_meta:
            self.property_item(key)
        return self.properties_dict

    def property_item(self, name: str) -> np.ndarray:
        """
        Retrieve one property by name.
       
        Returns
        -------
        numpy.ndarray
            a property value


        The name is one of these:

             "fid" (family id), "iid" (individual or sample id), "father" (father id),
             "mother" (mother id), "sex", "pheno" (phenotype), "chromosome", "sid"
             (SNP or variant id), "cm_position" (centimorgan position), "bp_position"
             (base-pair position), "allele_1", "allele_2".
            
        If needed, will cause a one-time read of the .fam or .bim file.

        Example
        -------

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.property_item('chromosome'))
            ['1' '1' '5' 'Y']

        """
        if name not in self.properties_dict:
            mm = _meta_meta[name]
            self._read_fam_or_bim(suffix=mm.suffix)
        return self.properties_dict[name]

    @property
    def chromosome(self) -> np.ndarray:
        """
        Chromosome of each SNP (variant)
       
        Returns
        -------
        numpy.ndarray
            array of str


        '0' represents a missing value.
        
        If needed, will cause a one-time read of the .bim file.

        Example
        -------

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.chromosome)
            ['1' '1' '5' 'Y']

        """
        return self.property_item("chromosome")

    @property
    def sid(self) -> np.ndarray:
        """
        SNP id of each SNP (variant).
       
        Returns
        -------
        numpy.ndarray
            array of str


        If needed, will cause a one-time read of the .bim file.

        Example
        -------

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.sid)
            ['sid1' 'sid2' 'sid3' 'sid4']

        """
        return self.property_item("sid")

    @property
    def cm_position(self) -> np.ndarray:
        """
        Centimorgan position of each SNP (variant).
       
        Returns
        -------
        numpy.ndarray
            array of float


        0.0 represents a missing value.
        
        If needed, will cause a one-time read of the .bim file.

        Example
        -------

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.cm_position)
            [ 100.4 2000.5 4000.7 7000.9]

        """
        return self.property_item("cm_position")

    @property
    def bp_position(self) -> np.ndarray:
        """
        Base-pair position of each SNP (variant).
       
        Returns
        -------
        numpy.ndarray
            array of int


        0 represents a missing value.
        
        If needed, will cause a one-time read of the .bim file.

        Example
        -------

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.bp_position)
            [   1  100 1000 1004]

        """
        return self.property_item("bp_position")

    @property
    def allele_1(self) -> np.ndarray:
        """
        First allele of each SNP (variant).
       
        Returns
        -------
        numpy.ndarray
            array of str


        If needed, will cause a one-time read of the .bim file.

        Example
        -------

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.allele_1)
            ['A' 'T' 'A' 'T']

        """
        return self.property_item("allele_1")

    @property
    def allele_2(self) -> np.ndarray:
        """
        Second allele of each SNP (variant),
       
        Returns
        -------
        numpy.ndarray
            array of str


        If needed, will cause a one-time read of the .bim file.

        Example
        -------

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.allele_2)
            ['A' 'C' 'C' 'G']

        """
        return self.property_item("allele_2")

    @property
    def iid_count(self) -> np.ndarray:
        """
        Number of individuals (samples).
       
        Returns
        -------
        int
            number of individuals


        If needed, will cause a fast line-count of the .fam file.

        Example
        -------

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.iid_count)
            3

        """
        return self._count("fam")

    @property
    def sid_count(self) -> np.ndarray:
        """
        Number of SNPs (variants).
       
        Returns
        -------
        int
            number of SNPs


        If needed, will cause a fast line-count of the .bim file.

        Example
        -------

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.sid_count)
            4

        """
        return self._count("bim")

    def _property_filepath(self, suffix):
        if suffix == "fam":
            return self._fam_filepath
        else:
            assert suffix == "bim"  # real assert
            return self._bim_filepath

    def _count(self, suffix):
        count = self._counts[suffix]
        if count is None:
            count = _rawincount(self._property_filepath(suffix))
            self._counts[suffix] = count
        return count

    @staticmethod
    def _check_file(filepointer):
        mode = filepointer.read(2)
        if mode != b"l\x1b":
            raise ValueError("Not a valid .bed file")
        mode = filepointer.read(1)  # \x01 = SNP major \x00 = individual major
        if mode != b"\x01":
            raise ValueError("only SNP-major is implemented")

    def __del__(self):
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        pass

    def _get_num_threads(self):
        if platform.system() == "Darwin":
            return 1
        if self._num_threads is not None:
            return self._num_threads
        if "MKL_NUM_THREADS" in os.environ:
            return int(os.environ["MKL_NUM_THREADS"])
        return multiprocessing.cpu_count()

    @staticmethod
    def _array_properties_are_ok(val, order):

        if order == "F":
            return val.flags["F_CONTIGUOUS"]
        else:
            assert order == "C"  # real assert
            return val.flags["C_CONTIGUOUS"]

    @property
    def shape(self):
        """
        Number of individuals (samples) and SNPs (variants).
       
        Returns
        -------
        (int, int)
            number of individuals, number of SNPs


        If needed, will cause a fast line-count of the .fam and .bim files.

        Example
        -------

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.shape)
            (3, 4)

        """
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
    def _write_fam_or_bim(base_filepath, properties, suffix, property_filepath):
        assert suffix in {"fam", "bim"}, "real assert"

        filepath = (
            Path(property_filepath)
            if property_filepath is not None
            else base_filepath.parent / (base_filepath.stem + "." + suffix)
        )

        fam_bim_list = []
        for key, mm in _meta_meta.items():
            if mm.suffix == suffix:
                assert len(fam_bim_list) == mm.column, "real assert"
                fam_bim_list.append(properties[key])

        sep = " " if suffix == "fam" else "\t"

        with open(filepath, "w") as filepointer:
            for index in range(len(fam_bim_list[0])):
                filepointer.write(
                    sep.join(str(seq[index]) for seq in fam_bim_list) + "\n"
                )

    @staticmethod
    def _fix_up_properties_array(input, dtype, missing_value, key):
        if input is None:
            return None
        if len(input) == 0:
            return np.zeros([0], dtype=dtype)

        if not isinstance(input, np.ndarray):
            return open_bed._fix_up_properties_array(
                np.array(input), dtype, missing_value, key
            )

        if len(input.shape) != 1:
            raise ValueError(f"{key} should be one dimensional")

        if not np.issubdtype(input.dtype, dtype):
            output = np.array(input, dtype=dtype)
        else:
            output = input

        # Change NaN in input to correct missing value
        if np.issubdtype(input.dtype, np.floating):
            output[input != input] = missing_value

        return output

    @staticmethod
    def _fix_up_properties(properties, iid_count, sid_count, use_fill_sequence):
        for key in properties:
            if key not in _meta_meta:
                raise KeyError(f"properties key '{key}' not known")

        count_dict = {"fam": iid_count, "bim": sid_count}
        properties_dict = {}
        for key, mm in _meta_meta.items():
            count = count_dict[mm.suffix]

            if key not in properties or (use_fill_sequence and properties[key] is None):
                if use_fill_sequence:
                    output = mm.fill_sequence(key, count, mm.missing_value, mm.dtype)
                else:
                    continue  # Test coverage reaches this, but doesn't report it.
            else:
                output = open_bed._fix_up_properties_array(
                    properties[key], mm.dtype, mm.missing_value, key
                )

            if output is not None:
                if count is None:
                    count_dict[mm.suffix] = len(output)
                else:
                    if count != len(output):
                        raise ValueError(
                            f"The length of override {key}, {len(output)}, should not be different from the current {_count_name[mm.suffix]}, {count}"
                        )
            properties_dict[key] = output
        return properties_dict, count_dict

    def _read_fam_or_bim(self, suffix):
        property_filepath = self._property_filepath(suffix)

        logging.info("Loading {0} file {1}".format(suffix, property_filepath))

        count = self._counts[suffix]

        delimiter = _delimiters[suffix]
        if delimiter in {r"\s+"}:
            delimiter = None
            delim_whitespace = True
        else:
            delim_whitespace = False

        usecolsdict = {}
        dtype_dict = {}
        for key, mm in _meta_meta.items():
            if mm.suffix is suffix and key not in self.properties_dict:
                usecolsdict[key] = mm.column
                dtype_dict[mm.column] = mm.dtype
        assert list(usecolsdict.values()) == sorted(usecolsdict.values())  # real assert
        assert len(usecolsdict) > 0  # real assert

        if os.path.getsize(property_filepath) == 0:
            fields = []
        else:
            fields = pd.read_csv(
                property_filepath,
                delimiter=delimiter,
                delim_whitespace=delim_whitespace,
                header=None,
                index_col=False,
                comment=None,
                dtype=dtype_dict,
                usecols=usecolsdict.values(),
            )

        if count is None:
            self._counts[suffix] = len(fields)
        else:
            if count != len(fields):
                raise ValueError(
                    f"The number of lines in the *.{suffix} file, {len(fields)}, should not be different from the current {_count_name[suffix]}, {count}"
                )
        for key in usecolsdict.keys():
            mm = _meta_meta[key]
            if len(fields) == 0:
                output = np.array([], dtype=mm.dtype)
            else:
                output = fields[mm.column].values
            self.properties_dict[key] = output


#if __name__ == "__main__":
    #logging.basicConfig(level=logging.INFO)
    #import os

    #if False:
    #    from bed_reader import open_bed, sample_file

    #    file_name = sample_file("small.bed")
    #    with open_bed(file_name) as bed:
    #        print(bed.iid)
    #        print(bed.sid)
    #        print(bed.read())

    # if False:
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

    #import pytest

    #pytest.main(["--doctest-modules", __file__])
