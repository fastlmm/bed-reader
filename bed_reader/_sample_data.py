import tempfile
from pathlib import Path, PurePath
from typing import Union

try:
    import pooch

    """
    Load sample data.
    """

    POOCH = pooch.create(
        # Use the default cache folder for the OS
        path=pooch.os_cache("bed_reader"),
        # The remote data is on Github
        base_url="https://raw.githubusercontent.com/"
        + "fastlmm/bed-sample-files/main/",
        # If this is a development version, get the data from the master branch
        version_dev="main",
        # The registry specifies the files that can be fetched
        env="BED_READER_DATA_DIR",
    )

    # Get registry file from package_data
    registry_file = Path(__file__).parent / "tests/registry.txt"
    # Load this registry file
    POOCH.load_registry(registry_file)
except ImportError:
    pooch = None


def sample_file(filepath: Union[str, Path]) -> str:
    """
    Retrieve a sample .bed file. (Also retrieves associated .fam and .bim files).

    Parameters
    ----------
    filepath
        Name of the sample .bed file.

    Returns
    -------
    str
        Local name of sample .bed file.


    .. note::
        This function requires the :mod:`pooch` package. Install `pooch` with:

        .. code-block:: bash

            pip install --upgrade bed-reader[samples]


    By default this function puts files under the user's cache directory.
    Override this by setting
    the `BED_READER_DATA_DIR` environment variable.

    Example
    --------

    .. doctest::

        >>> # pip install bed-reader[samples]  # if needed
        >>> from bed_reader import sample_file
        >>>
        >>> file_name = sample_file("small.bed")
        >>> print(f"The local file name is '{file_name}'")
        The local file name is '...small.bed'
    """
    if pooch is None:
        raise ImportError(
            "The function sample_file() requires pooch. "
            + "Install it with 'pip install --upgrade bed-reader[samples]'."
        )

    filepath = Path(filepath)
    file_string = str(filepath)
    if file_string.lower().endswith(".bed"):
        POOCH.fetch(file_string[:-4] + ".fam")
        POOCH.fetch(file_string[:-4] + ".bim")
    return POOCH.fetch(file_string)


def sample_url(filepath: Union[str, Path]) -> str:
    """
    Retrieve a URL to a sample .bed file. (Also makes ready associated .fam and .bim files).

    Parameters
    ----------
    filepath
        Name of the sample .bed file.

    Returns
    -------
    str
        URL to sample .bed file.


    .. note::
        This function requires the :mod:`pooch` package. Install `pooch` with:

        .. code-block:: bash

            pip install --upgrade bed-reader[samples]


    By default this function puts files under the user's cache directory.
    Override this by setting
    the `BED_READER_DATA_DIR` environment variable.

    Example
    --------

    .. doctest::

        >>> # pip install bed-reader[samples]  # if needed
        >>> from bed_reader import sample_url
        >>>
        >>> url = sample_url("small.bed")
        >>> print(f"The url is '{url}'")
        The url is 'file:///.../small.bed'
    """
    file_name = sample_file(filepath)
    url = PurePath(file_name).as_uri()
    return url


def tmp_path() -> Path:
    """
    Return a :class:`pathlib.Path` to a temporary directory.

    Returns
    -------
    pathlib.Path
        a temporary directory

    Example
    --------

    .. doctest::

        >>> from bed_reader import to_bed, tmp_path
        >>>
        >>> output_file = tmp_path() / "small3.bed"
        >>> val = [[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]]
        >>> to_bed(output_file, val)

    """
    temp_dir = tempfile.gettempdir()
    path = Path(temp_dir) / "bed_reader_tmp_path"
    path.mkdir(parents=True, exist_ok=True)
    return path


# if __name__ == "__main__":
#    logging.basicConfig(level=logging.INFO)

#    import pytest

#    pytest.main(["--doctest-modules", __file__])
