import logging
import multiprocessing
import os
import re
from dataclasses import dataclass
from io import BytesIO
from itertools import repeat, takewhile
from pathlib import Path
from typing import Any, List, Mapping, Optional, Union
from urllib.parse import ParseResult as UrlParseResult
from urllib.parse import urlparse

import numpy as np

try:
    from scipy import sparse
except ImportError:
    sparse = None

from .bed_reader import (  # type: ignore
    check_file_cloud,
    read_cloud_f32,
    read_cloud_f64,
    read_cloud_i8,
    read_f32,
    read_f64,
    read_i8,
    url_to_bytes,
)


# https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python
def _rawincount(f):
    f.seek(0)
    bufgen = takewhile(lambda x: x, (f.read(1024 * 1024) for _ in repeat(None)))
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
        (f"{key}{i + 1}" for i in range(length)), dtype=dtype, count=length
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


def get_num_threads(num_threads=None):
    if num_threads is not None:
        return num_threads
    if "PST_NUM_THREADS" in os.environ:
        return int(os.environ["PST_NUM_THREADS"])
    if "NUM_THREADS" in os.environ:
        return int(os.environ["NUM_THREADS"])
    if "MKL_NUM_THREADS" in os.environ:
        return int(os.environ["MKL_NUM_THREADS"])
    return multiprocessing.cpu_count()


def get_max_concurrent_requests(max_concurrent_requests=None):
    if max_concurrent_requests is not None:
        return max_concurrent_requests
    return 10


def get_max_chunk_bytes(max_chunk_bytes=None):
    if max_chunk_bytes is not None:
        return max_chunk_bytes
    return 8_000_000


class open_bed:
    """
    Open a PLINK .bed file, local or cloud, for reading.

    Parameters
    ----------
    location: pathlib.Path or str
        File path or URL to the .bed file. See :doc:`cloud_urls` for details on cloud URLs.
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
          of the appropriate dtype, if necessary. Any :data:`numpy.nan` values
          will converted to the appropriate missing value. The PLINK `.fam specification
          <https://www.cog-genomics.org/plink2/formats#fam>`_
          and `.bim specification <https://www.cog-genomics.org/plink2/formats#bim>`_
          lists the dtypes and missing values for each property.

    count_A1: bool, optional
        True (default) to count the number of A1 alleles (the PLINK standard).
        False to count the number of A2 alleles.
    num_threads: None or int, optional
        The number of threads with which to read data. Defaults to all available
        processors.
        Can also be set with these environment variables (listed in priority order):
        'PST_NUM_THREADS', 'NUM_THREADS', 'MKL_NUM_THREADS'.
    skip_format_check: bool, optional
        False (default) to immediately check for expected starting bytes in
        the .bed file. True to delay the check until (and if) data is read.
    fam_location: pathlib.Path or str or URL, optional
        Path to the file containing information about each individual (sample).
        Defaults to replacing the .bed file’s suffix with .fam.
    bim_location: pathlib.Path or str URL, optional
        Path to the file containing information about each SNP (variant).
        Defaults to replacing the .bed file’s suffix with .bim.
    cloud_options: dict, optional
        A dictionary of options for reading from cloud storage. The default is an empty.
    max_concurrent_requests: None or int, optional
        The maximum number of concurrent requests to make to the cloud storage service.
        Defaults to 10.
    max_chunk_bytes: None or int, optional
        The maximum number of bytes to read in a single request to the cloud storage
        service. Defaults to 8MB.
    filepath: same as location
        Deprecated. Use location instead.
    fam_filepath: same as fam_location
        Deprecated. Use fam_location instead.
    bim_filepath: same as bim_location
        Deprecated. Use bim_location instead.

    Returns
    -------
    open_bed
        an open_bed object

    Examples
    --------

    Open a local file and list individual (sample) :attr:`iid` and SNP (variant) :attr:`sid`. Then, :meth:`read`
    the whole file.

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

    Open a cloud file with a non-default timeout.
    Then, read the data for one SNP (variant)
    at index position 2.

    See :doc:`cloud_urls` for details on reading files from cloud storage.

    .. doctest::

        >>> import numpy as np
        >>> with open_bed("https://raw.githubusercontent.com/fastlmm/bed-sample-files/main/small.bed",
        ...               cloud_options={"timeout": "10s"}) as bed:
        ...     print(bed.read(np.s_[:,2]))
        [[nan]
         [nan]
         [ 2.]]

    With the local file, replace :attr:`iid`.


        >>> bed = open_bed(file_name, properties={"iid":["sample1","sample2","sample3"]})
        >>> print(bed.iid) # replaced
        ['sample1' 'sample2' 'sample3']
        >>> print(bed.sid) # same as before
        ['sid1' 'sid2' 'sid3' 'sid4']

    Give the number of individuals (samples) and SNPs (variants) so that the .fam and
    .bim files need never be opened.

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
        location: Union[str, Path, UrlParseResult],
        iid_count: Optional[int] = None,
        sid_count: Optional[int] = None,
        properties: Mapping[str, List[Any]] = {},
        count_A1: bool = True,
        num_threads: Optional[int] = None,
        skip_format_check: bool = False,
        fam_location: Union[str, Path, UrlParseResult] = None,
        bim_location: Union[str, Path, UrlParseResult] = None,
        cloud_options: Mapping[str, str] = {},
        max_concurrent_requests: Optional[int] = None,
        max_chunk_bytes: Optional[int] = None,
        # accept old keywords
        filepath: Union[str, Path] = None,
        fam_filepath: Union[str, Path] = None,
        bim_filepath: Union[str, Path] = None,
    ):
        location = self._combined(location, filepath, "location", "filepath")
        fam_location = self._combined(
            fam_location, fam_filepath, "fam_location", "fam_filepath"
        )
        bim_location = self._combined(
            bim_location, bim_filepath, "bim_location", "bim_filepath"
        )

        self.location = self._path_or_url(location)
        self.cloud_options = cloud_options
        self.count_A1 = count_A1
        self._num_threads = num_threads
        self._max_concurrent_requests = max_concurrent_requests
        self._max_chunk_bytes = max_chunk_bytes
        self.skip_format_check = skip_format_check
        self._fam_location = (
            self._path_or_url(fam_location)
            if fam_location is not None
            else self._replace_extension(self.location, "fam")
        )
        self._bim_location = (
            self._path_or_url(bim_location)
            if bim_location is not None
            else self._replace_extension(self.location, "bim")
        )

        self.properties_dict, self._counts = open_bed._fix_up_properties(
            properties, iid_count, sid_count, use_fill_sequence=False
        )
        self._iid_range = None
        self._sid_range = None
        if not self.skip_format_check:
            if self._is_url(self.location):
                check_file_cloud(self.location.geturl(), self.cloud_options)
            else:
                with open(self.location, "rb") as filepointer:
                    self._mode = self._check_file(filepointer)

    # # its an error to set both location and filepath
    # location = self._combined(location, filepath, "location", "filepath")
    # fam_location = self._combined(fam_location, fam_filepath, "fam_location", "fam_filepath")
    # bim_location = self._combined(bim_location, bim_filepath, "bim_location", "bim_filepath")
    @staticmethod
    def _combined(location, filepath, location_name, filepath_name):
        if location is not None and filepath is not None:
            raise ValueError(f"Cannot set both {location_name} and {filepath_name}")
        # None, None is ok for now
        return location if location is not None else filepath

    @staticmethod
    def _replace_extension(location, extension):
        if open_bed._is_url(location):
            # Split the path and change the extension
            path, _ = os.path.splitext(location.path)
            new_path = f"{path}.{extension}"

            # Create a new ParseResult with the updated path
            new_parse_result = UrlParseResult(
                scheme=location.scheme,
                netloc=location.netloc,
                path=new_path,
                params=location.params,
                query=location.query,
                fragment=location.fragment,
            )
            return new_parse_result
        else:
            assert isinstance(location, Path)  # real assert
            return location.parent / (location.stem + "." + extension)

    @staticmethod
    def _is_url(location):
        return isinstance(location, UrlParseResult)

    @staticmethod
    def _path_or_url(input):
        if isinstance(input, Path):
            return input
        if isinstance(input, UrlParseResult):
            return input
        assert isinstance(
            input, str
        ), "Expected a string or Path object or UrlParseResult"
        parsed = urlparse(input)
        if parsed.scheme and "://" in input:
            return parsed
        else:
            return Path(input)

    def read(
        self,
        index: Optional[Any] = None,
        dtype: Optional[Union[type, str]] = "float32",
        order: Optional[str] = "F",
        force_python_only: Optional[bool] = False,
        num_threads=None,
        max_concurrent_requests=None,
        max_chunk_bytes=None,
    ) -> np.ndarray:
        """
        Read genotype information.

        Parameters
        ----------
        index:
            An optional expression specifying the individuals (samples) and SNPs
            (variants) to read. (See examples, below).
            Defaults to ``None``, meaning read all.

            (If index is a tuple, the first component indexes the individuals and the
            second indexes
            the SNPs. If it is not a tuple and not None, it indexes SNPs.)

        dtype: {'float32' (default), 'float64', 'int8'}, optional
            The desired data-type for the returned array.
        order : {'F','C'}, optional
            The desired memory layout for the returned array.
            Defaults to ``F`` (Fortran order, which is SNP-major).
        force_python_only: bool, optional
            If False (default), uses the faster Rust code; otherwise it uses the slower
            pure Python code.

        num_threads: None or int, optional
            The number of threads with which to read data. Defaults to all available
            processors.
            Can also be set with :class:`open_bed` or these
            environment variables (listed in priority order):
            'PST_NUM_THREADS', 'NUM_THREADS', 'MKL_NUM_THREADS'.

        max_concurrent_requests: None or int, optional
            The maximum number of concurrent requests to make to the cloud storage
            service. Defaults to 10.

        max_chunk_bytes: None or int, optional
            The maximum number of bytes to read in a single request to the cloud
            storage service. Defaults to 8MB.

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

        To read selected individuals (samples) and/or SNPs (variants), set each part of
        a :data:`numpy.s_` to an `int`, a list of `int`, a slice expression, or
        a list of `bool`.
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
            >>> # read SNPs from 1 (inclusive) to 4 (exclusive)
            >>> print(bed.read(np.s_[:,1:4]))
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

        # Later happy with _iid_range and _sid_range or could it be done with
        # allocation them?
        if self._iid_range is None:
            self._iid_range = np.arange(self.iid_count, dtype="intp")
        if self._sid_range is None:
            self._sid_range = np.arange(self.sid_count, dtype="intp")

        iid_index = np.ascontiguousarray(
            self._iid_range[iid_index_or_slice_etc],
            dtype="intp",
        )
        sid_index = np.ascontiguousarray(
            self._sid_range[sid_index_or_slice_etc], dtype="intp"
        )

        if not force_python_only or open_bed._is_url(self.location):
            num_threads = get_num_threads(
                self._num_threads if num_threads is None else num_threads
            )
            max_concurrent_requests = get_max_concurrent_requests(
                self._max_concurrent_requests
                if max_concurrent_requests is None
                else max_concurrent_requests
            )
            max_chunk_bytes = get_max_chunk_bytes(
                self._max_chunk_bytes if max_chunk_bytes is None else max_chunk_bytes
            )

            val = np.zeros((len(iid_index), len(sid_index)), order=order, dtype=dtype)

            if self.iid_count > 0 and self.sid_count > 0:
                reader, location_str, is_cloud = self._pick_reader(dtype)

                if not is_cloud:
                    reader(
                        location_str,
                        self.cloud_options,
                        iid_count=self.iid_count,
                        sid_count=self.sid_count,
                        is_a1_counted=self.count_A1,
                        iid_index=iid_index,
                        sid_index=sid_index,
                        val=val,
                        num_threads=num_threads,
                    )
                else:
                    reader(
                        location_str,
                        self.cloud_options,
                        iid_count=self.iid_count,
                        sid_count=self.sid_count,
                        is_a1_counted=self.count_A1,
                        iid_index=iid_index,
                        sid_index=sid_index,
                        val=val,
                        num_threads=num_threads,
                        max_concurrent_requests=max_concurrent_requests,
                        max_chunk_bytes=max_chunk_bytes,
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
            # An earlier version of this code had a way to read consecutive SNPs of code
            # in one read. May want
            # to add that ability back to the code.
            # Also, note that reading with python will often result in
            # non-contiguous memory
            # logging.warn("using pure python plink parser (might be much slower!!)")

            if self.major == "SNP":
                minor_count = self.iid_count
                minor_index = iid_index
                major_index = sid_index
            else:
                minor_count = self.sid_count
                minor_index = sid_index
                major_index = iid_index

            val = np.zeros(
                ((int(np.ceil(0.25 * minor_count)) * 4), len(major_index)),
                order=order,
                dtype=dtype,
            )  # allocate it a little big

            nbyte = int(np.ceil(0.25 * minor_count))
            with open(self.location, "rb") as filepointer:
                for major_index_value, major_index_index in enumerate(major_index):
                    startbit = int(np.ceil(0.25 * minor_count) * major_index_index + 3)
                    filepointer.seek(startbit)
                    bytes = np.array(bytearray(filepointer.read(nbyte))).reshape(
                        (int(np.ceil(0.25 * minor_count)), 1), order="F"
                    )

                    val[3::4, major_index_value : major_index_value + 1] = byteZero
                    val[3::4, major_index_value : major_index_value + 1][
                        bytes >= 64
                    ] = missing
                    val[3::4, major_index_value : major_index_value + 1][
                        bytes >= 128
                    ] = 1
                    val[3::4, major_index_value : major_index_value + 1][
                        bytes >= 192
                    ] = byteThree
                    bytes = np.mod(bytes, 64)
                    val[2::4, major_index_value : major_index_value + 1] = byteZero
                    val[2::4, major_index_value : major_index_value + 1][
                        bytes >= 16
                    ] = missing
                    val[2::4, major_index_value : major_index_value + 1][
                        bytes >= 32
                    ] = 1
                    val[2::4, major_index_value : major_index_value + 1][
                        bytes >= 48
                    ] = byteThree
                    bytes = np.mod(bytes, 16)
                    val[1::4, major_index_value : major_index_value + 1] = byteZero
                    val[1::4, major_index_value : major_index_value + 1][
                        bytes >= 4
                    ] = missing
                    val[1::4, major_index_value : major_index_value + 1][bytes >= 8] = 1
                    val[1::4, major_index_value : major_index_value + 1][
                        bytes >= 12
                    ] = byteThree
                    bytes = np.mod(bytes, 4)
                    val[0::4, major_index_value : major_index_value + 1] = byteZero
                    val[0::4, major_index_value : major_index_value + 1][
                        bytes >= 1
                    ] = missing
                    val[0::4, major_index_value : major_index_value + 1][bytes >= 2] = 1
                    val[0::4, major_index_value : major_index_value + 1][
                        bytes >= 3
                    ] = byteThree
                val = val[minor_index, :]  # reorder or trim any extra allocation
                assert val.dtype == np.dtype(dtype)  # real assert
                if not open_bed._array_properties_are_ok(val, order):
                    val = val.copy(order=order)
            # if in force python mode, and individual-major mode, then we need to transpose
            if self.major == "individual":
                val = val.T
        return val

    def _pick_reader(self, dtype):
        if dtype == np.int8:
            file_reader = read_i8
            cloud_reader = read_cloud_i8
        elif dtype == np.float64:
            file_reader = read_f64
            cloud_reader = read_cloud_f64
        elif dtype == np.float32:
            file_reader = read_f32
            cloud_reader = read_cloud_f32
        else:
            raise ValueError(
                f"dtype '{dtype}' not known, only "
                + "'int8', 'float32', and 'float64' are allowed."
            )

        if open_bed._is_url(self.location):
            reader = cloud_reader
            location_str = self.location.geturl()
            is_cloud = True
        else:
            reader = file_reader
            location_str = str(self.location.as_posix())
            is_cloud = False
        return reader, location_str, is_cloud

    def __str__(self) -> str:
        return f"{self.__class__.__name__}('{self.location}',...)"

    @property
    def major(self) -> str:
        """
        Major mode of a local .bed file.

        Returns
        -------
        str
            'SNP' or 'individual'


        Almost all PLINK 1.9 .bed files are 'SNP' major. This makes
        reading the data by SNP(s) fast.

        Errors
        ------
        ValueError
            If the file is a cloud file.

        Example
        -------

        .. doctest::

            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("small.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.major)
            SNP

        """
        if self._is_url(self.location):
            raise ValueError("Cannot determine major mode for cloud files")
        # if self._mode is not set, set it
        if not hasattr(self, "mode"):
            with open(self.location, "rb") as filepointer:
                self._mode = self._check_file(filepointer)

        return "individual" if self._mode == b"\x00" else "SNP"

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

    def _property_location(self, suffix):
        if suffix == "fam":
            return self._fam_location
        else:
            assert suffix == "bim"  # real assert
            return self._bim_location

    def _count(self, suffix):
        count = self._counts[suffix]
        if count is None:
            location = self._property_location(suffix)
            if open_bed._is_url(location):
                # should not download twice from cloud
                if suffix == "fam":
                    if self.property_item("iid") is None:
                        # ... unless user doesn't want iid
                        file_bytes = bytes(
                            url_to_bytes(location.geturl(), self.cloud_options)
                        )
                        count = _rawincount(BytesIO(file_bytes))
                    else:
                        count = len(self.iid)
                elif suffix == "bim":
                    if self.property_item("sid") is None:
                        # ... unless user doesn't want sid
                        file_bytes = bytes(
                            url_to_bytes(location.geturl(), self.cloud_options)
                        )
                        count = _rawincount(BytesIO(file_bytes))
                    else:
                        count = len(self.sid)
                else:
                    raise ValueError("real assert")
            else:
                count = _rawincount(open(location, "rb"))
                self._counts[suffix] = count
        return count

    @staticmethod
    def _check_file(filepointer):
        magic_number = filepointer.read(2)
        if magic_number != b"l\x1b":
            raise ValueError("Not a valid .bed file")
        mode = filepointer.read(1)
        # Check if mode is either individual-major or SNP-major
        if mode not in (b"\x00", b"\x01"):
            raise ValueError("Not a valid .bed file")
        return mode

    def __del__(self):
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        pass

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
        return (self.iid_count, self.sid_count)

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

        do_missing_values = True
        if np.issubdtype(input.dtype, np.floating) and np.issubdtype(dtype, int):
            input[input != input] = missing_value
            old_settings = np.seterr(invalid="warn")
            try:
                output = np.array(input, dtype=dtype)
            finally:
                np.seterr(**old_settings)
        elif not np.issubdtype(input.dtype, dtype):
            # This will convert, for example, numerical sids to string sids or
            # floats that happen to be integers into ints,
            # but there will be a warning generated.
            old_settings = np.seterr(invalid="warn")
            try:
                output = np.array(input, dtype=dtype)
            finally:
                np.seterr(**old_settings)
        else:
            output = input

        # Change NaN in input to correct missing value
        if do_missing_values and np.issubdtype(input.dtype, np.floating):
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
                            f"The length of override {key}, {len(output)}, should not "
                            + "be different from the current "
                            + f"{_count_name[mm.suffix]}, {count}"
                        )
            properties_dict[key] = output
        return properties_dict, count_dict

    def _read_fam_or_bim(self, suffix):
        property_location = self._property_location(suffix)

        logging.info("Loading {0} file {1}".format(suffix, property_location))

        count = self._counts[suffix]

        delimiter = _delimiters[suffix]
        if delimiter in {r"\s+"}:
            delimiter = None

        usecolsdict = {}
        dtype_dict = {}
        for key, mm in _meta_meta.items():
            if mm.suffix is suffix and key not in self.properties_dict:
                usecolsdict[key] = mm.column
                dtype_dict[mm.column] = mm.dtype
        assert list(usecolsdict.values()) == sorted(usecolsdict.values())  # real assert
        assert len(usecolsdict) > 0  # real assert

        if self._is_url(property_location):
            file_bytes = bytes(
                url_to_bytes(property_location.geturl(), self.cloud_options)
            )
            if len(file_bytes) == 0:
                columns, row_count = [], 0
            else:  # note similar code below
                columns, row_count = _read_csv(
                    BytesIO(file_bytes),
                    delimiter=delimiter,
                    dtype=dtype_dict,
                    usecols=usecolsdict.values(),
                )
        else:
            if os.path.getsize(property_location) == 0:
                columns, row_count = [], 0
            else:
                columns, row_count = _read_csv(
                    property_location,
                    delimiter=delimiter,
                    dtype=dtype_dict,
                    usecols=usecolsdict.values(),
                )

        if count is None:
            self._counts[suffix] = row_count
        else:
            if count != row_count:
                raise ValueError(
                    f"The number of lines in the *.{suffix} file, {row_count}, "
                    + "should not be different from the current "
                    + "f{_count_name[suffix]}, {count}"
                )
        for i, key in enumerate(usecolsdict.keys()):
            mm = _meta_meta[key]
            if row_count == 0:
                output = np.array([], dtype=mm.dtype)
            else:
                output = columns[i]
                if not np.issubdtype(output.dtype, mm.dtype):
                    output = np.array(output, dtype=mm.dtype)
            self.properties_dict[key] = output

    def read_sparse(
        self,
        index: Optional[Any] = None,
        dtype: Optional[Union[type, str]] = "float32",
        batch_size: Optional[int] = None,
        format: Optional[str] = "csc",
        num_threads=None,
        max_concurrent_requests=None,
        max_chunk_bytes=None,
    ) -> (Union[sparse.csc_matrix, sparse.csr_matrix]) if sparse is not None else None:  # type: ignore
        """
        Read genotype information into a :mod:`scipy.sparse` matrix. Sparse matrices
        may be useful when the data is mostly zeros.

        .. note::
            This method requires :mod:`scipy`. Install `scipy` with:

            .. code-block:: bash

                pip install --upgrade bed-reader[sparse]

        Parameters
        ----------
        index:
            An optional expression specifying the individuals (samples) and SNPs
            (variants) to read. (See examples, below).
            Defaults to ``None``, meaning read all.

            (If index is a tuple, the first component indexes the individuals and the
            second indexes
            the SNPs. If it is not a tuple and not None, it indexes SNPs.)

        dtype: {'float32' (default), 'float64', 'int8'}, optional
            The desired data-type for the returned array.
        batch_size: None or int, optional
            Number of dense columns or rows to read at a time, internally.
            Defaults to round(sqrt(total-number-of-columns-or-rows-to-read)).
        format : {'csc','csr'}, optional
            The desired format of the sparse matrix.
            Defaults to ``csc`` (Compressed Sparse Column, which is SNP-major).
        num_threads: None or int, optional
            The number of threads with which to read data. Defaults to all available
            processors.
            Can also be set with :class:`open_bed` or these
            environment variables (listed in priority order):
            'PST_NUM_THREADS', 'NUM_THREADS', 'MKL_NUM_THREADS'.
        max_concurrent_requests: None or int, optional
            The maximum number of concurrent requests to make to the cloud storage
            service. Defaults to 10.
        max_chunk_bytes: None or int, optional
            The maximum number of bytes to read in a single request to the cloud
            storage service. Defaults to 8MB.


        Returns
        -------
        a :class:`scipy.sparse.csc_matrix` (default) or :class:`scipy.sparse.csr_matrix`

        Rows represent individuals (samples). Columns represent SNPs (variants).

        For ``dtype`` 'float32' and 'float64', NaN indicates missing values.
        For 'int8', -127 indicates missing values.


        The memory used by the final sparse matrix is approximately:

           # of non-zero values * (4 bytes + 1 byte (for int8))

        For example, consider reading 1000 individuals (samples) x 50,000 SNPs (variants)
        into csc format where the data is 97% sparse.
        The memory used will be about 7.5 MB (1000 x 50,000 x 3% x 5 bytes).
        This is 15% of the 50 MB needed by a dense matrix.

        Internally, the function reads the data via small dense matrices.
        For this example, by default, the function will read 1000 individuals x 224 SNPs
        (because 224 * 224 is about 50,000).
        The memory used by the small dense matrix is 1000 x 244 x 1 byte (for int8) = 0.224 MB.

        You can set `batch_size`. Larger values will be faster.
        Smaller values will use less memory.

        For this example, we might want to set the `batch_size` to 5000. Then,
        the memory used by the small dense matrix
        would be 1000 x 5000 x 1 byte (for int8) = 5 MB,
        similar to the 7.5 MB needed for the final sparse matrix.

        Examples
        --------

        Read all data in a .bed file into a :class:`scipy.sparse.csc_matrix`.
        The file has 10 individuals (samples) by 20 SNPs (variants).
        All but eight values are 0.

        .. doctest::

            >>> # pip install bed-reader[samples,sparse]  # if needed
            >>> from bed_reader import open_bed, sample_file
            >>>
            >>> file_name = sample_file("sparse.bed")
            >>> with open_bed(file_name) as bed:
            ...     print(bed.shape) # doctest:+NORMALIZE_WHITESPACE +ELLIPSIS
            ...     val_sparse = bed.read_sparse(dtype="int8")
            (10, 20)
            >>> print("Nonzero Values", val_sparse.data) # doctest:+NORMALIZE_WHITESPACE +ELLIPSIS
            Nonzero Values [1 2 2 1 1 1 1 1]

        To read selected individuals (samples) and/or SNPs (variants), set each part of
        a :data:`numpy.s_` to an `int`, a list of `int`, a slice expression, or
        a list of `bool`.
        Negative integers count from the end of the list.

        .. doctest::

            >>> import numpy as np
            >>> bed = open_bed(file_name)
            >>> print("Nonzero Values", bed.read_sparse(np.s_[:,5], dtype="int8").data)  # read the SNPs indexed by 5.
            Nonzero Values [2]
            >>> # read the SNPs indexed by 5, 4, and 0
            >>> print("Nonzero Values", bed.read_sparse(np.s_[:,[5,4,0]], dtype="int8").data)
            Nonzero Values [2 1]
            >>> # read SNPs from 1 (inclusive) to 11 (exclusive)
            >>> print("Nonzero Values", bed.read_sparse(np.s_[:,1:11], dtype="int8").data)
            Nonzero Values [1 2 2 1 1]
            >>> print(np.unique(bed.chromosome)) # print unique chrom values
            ['1' '5' 'Y']
            >>> # read all SNPs in chrom 5
            >>> print("Nonzero Values", bed.read_sparse(np.s_[:,bed.chromosome=='5'], dtype="int8").data)
            Nonzero Values [1 2 2 1 1 1 1 1]
            >>> # Read 1st individual (across all SNPs)
            >>> print("Nonzero Values", bed.read_sparse(np.s_[0,:], dtype="int8").data)
            Nonzero Values [2]
            >>> print("Nonzero Values", bed.read_sparse(np.s_[::2,:], dtype="int8").data) # Read every 2nd individual
            Nonzero Values [1 2 2 1 1]
            >>> # read last and 2nd-to-last individuals and the 15th-from-the-last SNP
            >>> print("Nonzero Values", bed.read_sparse(np.s_[[-1,-2],-15], dtype="int8").data)
            Nonzero Values [2]
        """
        if sparse is None:
            raise ImportError(
                "The function read_sparse() requires scipy. "
                + "Install it with 'pip install --upgrade bed-reader[sparse]'."
            )
        iid_index_or_slice_etc, sid_index_or_slice_etc = self._split_index(index)

        dtype = np.dtype(dtype)

        # Similar code in read().
        # Later happy with _iid_range and _sid_range or could it be done with
        # allocation them?
        if self._iid_range is None:
            self._iid_range = np.arange(self.iid_count, dtype="intp")
        if self._sid_range is None:
            self._sid_range = np.arange(self.sid_count, dtype="intp")

        iid_index = np.ascontiguousarray(
            self._iid_range[iid_index_or_slice_etc],
            dtype="intp",
        )
        sid_index = np.ascontiguousarray(
            self._sid_range[sid_index_or_slice_etc], dtype="intp"
        )

        if (
            len(iid_index) > np.iinfo(np.int32).max
            or len(sid_index) > np.iinfo(np.int32).max
        ):
            raise ValueError(
                "Too (many Individuals or SNPs (variants) requested. Maximum is {np.iinfo(np.int32).max}."
            )

        if batch_size is None:
            batch_size = round(np.sqrt(len(sid_index)))

        num_threads = get_num_threads(
            self._num_threads if num_threads is None else num_threads
        )
        max_concurrent_requests = get_max_concurrent_requests(
            self._max_concurrent_requests
            if max_concurrent_requests is None
            else max_concurrent_requests
        )
        max_chunk_bytes = get_max_chunk_bytes(
            self._max_chunk_bytes if max_chunk_bytes is None else max_chunk_bytes
        )

        if format == "csc":
            order = "F"
            indptr = np.zeros(len(sid_index) + 1, dtype=np.int32)
        elif format == "csr":
            order = "C"
            indptr = np.zeros(len(iid_index) + 1, dtype=np.int32)
        else:
            raise ValueError(f"format '{format}' not known. Expected 'csc' or 'csr'.")

        # We init data and indices with zero element arrays to set their dtype.
        data = [np.empty(0, dtype=dtype)]
        indices = [np.empty(0, dtype=np.int32)]

        if self.iid_count > 0 and self.sid_count > 0:
            reader, location_str, is_cloud = self._pick_reader(dtype)

            if format == "csc":
                val = np.zeros((len(iid_index), batch_size), order=order, dtype=dtype)
                for batch_start in range(0, len(sid_index), batch_size):
                    batch_end = batch_start + batch_size
                    if batch_end > len(sid_index):
                        batch_end = len(sid_index)
                        del val
                        val = np.zeros(
                            (len(iid_index), batch_end - batch_start),
                            order=order,
                            dtype=dtype,
                        )
                    batch_slice = np.s_[batch_start:batch_end]
                    batch_index = sid_index[batch_slice]

                    if not is_cloud:
                        reader(
                            location_str,
                            self.cloud_options,
                            iid_count=self.iid_count,
                            sid_count=self.sid_count,
                            is_a1_counted=self.count_A1,
                            iid_index=iid_index,
                            sid_index=batch_index,
                            val=val,
                            num_threads=num_threads,
                        )
                    else:
                        reader(
                            location_str,
                            self.cloud_options,
                            iid_count=self.iid_count,
                            sid_count=self.sid_count,
                            is_a1_counted=self.count_A1,
                            iid_index=iid_index,
                            sid_index=batch_index,
                            val=val,
                            num_threads=num_threads,
                            max_concurrent_requests=max_concurrent_requests,
                            max_chunk_bytes=max_chunk_bytes,
                        )

                    self.sparsify(
                        val, order, iid_index, batch_slice, data, indices, indptr
                    )
            else:
                assert format == "csr"  # real assert
                val = np.zeros((batch_size, len(sid_index)), order=order, dtype=dtype)
                for batch_start in range(0, len(iid_index), batch_size):
                    batch_end = batch_start + batch_size
                    if batch_end > len(iid_index):
                        batch_end = len(iid_index)
                        del val
                        val = np.zeros(
                            (batch_end - batch_start, len(sid_index)),
                            order=order,
                            dtype=dtype,
                        )

                    batch_slice = np.s_[batch_start:batch_end]
                    batch_index = iid_index[batch_slice]

                    if not is_cloud:
                        reader(
                            location_str,
                            self.cloud_options,
                            iid_count=self.iid_count,
                            sid_count=self.sid_count,
                            is_a1_counted=self.count_A1,
                            iid_index=batch_index,
                            sid_index=sid_index,
                            val=val,
                            num_threads=num_threads,
                        )
                    else:
                        reader(
                            location_str,
                            self.cloud_options,
                            iid_count=self.iid_count,
                            sid_count=self.sid_count,
                            is_a1_counted=self.count_A1,
                            iid_index=batch_index,
                            sid_index=sid_index,
                            val=val,
                            num_threads=num_threads,
                            max_concurrent_requests=max_concurrent_requests,
                            max_chunk_bytes=max_chunk_bytes,
                        )

                    self.sparsify(
                        val, order, sid_index, batch_slice, data, indices, indptr
                    )

        data = np.concatenate(data)
        indices = np.concatenate(indices)

        if format == "csc":
            return sparse.csc_matrix(
                (data, indices, indptr), (len(iid_index), len(sid_index))
            )
        else:
            assert format == "csr"  # real assert
            return sparse.csr_matrix(
                (data, indices, indptr), (len(iid_index), len(sid_index))
            )

    def sparsify(self, val, order, minor_index, batch_slice, data, indices, indptr):
        flatten = np.ravel(val, order=order)
        nz_indices = np.flatnonzero(flatten).astype(np.int32)
        column_indexes = nz_indices // len(minor_index)
        counts = np.bincount(
            column_indexes, minlength=batch_slice.stop - batch_slice.start
        ).astype(np.int32)
        counts_with_initial = np.r_[
            indptr[batch_slice.start : batch_slice.start + 1], counts
        ]

        data.append(flatten[nz_indices])
        indices.append(np.mod(nz_indices, len(minor_index)))
        indptr[1:][batch_slice] = np.cumsum(counts_with_initial)[1:]


def _read_csv(filepath, delimiter=None, dtype=None, usecols=None):
    pattern = re.compile(r"^np\.\w+\((.+?)\)$")

    # Prepare the usecols by ensuring it is a list of indices
    usecols_indices = list(usecols)
    transposed = np.loadtxt(
        filepath,
        dtype=np.str_,
        delimiter=delimiter,
        usecols=usecols_indices,
        unpack=True,
    )
    if transposed.ndim == 1:
        transposed = transposed.reshape(-1, 1)
    row_count = transposed.shape[1]  # because unpack=True

    # Convert column lists to numpy arrays with the specified dtype
    columns = []
    for output_index, input_index in enumerate(usecols_indices):
        col = transposed[output_index]

        # work around numpy/python bug
        if len(col) > 0 and pattern.fullmatch(col[0]):
            col = np.array([pattern.fullmatch(x).group(1) for x in col])

        # Find the dtype for this column
        col_dtype = dtype.get(input_index, np.str_)
        # Convert the column list to a numpy array with the specified dtype
        columns.append(_convert_to_dtype(col, col_dtype))

    return columns, row_count


def _convert_to_dtype(str_arr, dtype):
    assert dtype in [np.str_, np.float32, np.int32]  # real assert

    if dtype == np.str_:
        return str_arr

    try:
        new_arr = str_arr.astype(dtype)
    except ValueError as e:
        if dtype == np.float32:
            raise e
        # for backwards compatibility, see if intermediate float helps int conversion
        try:
            assert dtype == np.int32  # real assert
            float_arr = str_arr.astype(np.float32)
        except ValueError:
            raise e
        new_arr = float_arr.astype(np.int32)
        if not np.array_equal(new_arr, float_arr):
            raise ValueError(
                f"invalid literal for int: '{str_arr[np.where(new_arr != float_arr)][:1]}')"
            )
    return new_arr


if __name__ == "__main__":
    import pytest

    logging.basicConfig(level=logging.INFO)

    pytest.main(["--doctest-modules", __file__])
