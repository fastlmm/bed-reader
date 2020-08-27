from typing import Any, List, Optional, Tuple, Union
import pooch
from pathlib import Path

"""
Load sample data.
"""

POOCH = pooch.create(
    # Use the default cache folder for the OS
    path=pooch.os_cache("bed_reader"),
    # The remote data is on Github
    base_url="https://raw.githubusercontent.com/fastlmm/bed-reader/master/bed_reader/tests/data/",
    # If this is a development version, get the data from the master branch
    version_dev="master",
    # The registry specifies the files that can be fetched
    env="BED_READER_DATA_DIR",  #!!!cmk document this
)

# Get registry file from package_data
registry_file = Path(__file__).parent / "tests/registry.txt"
# Load this registry file
POOCH.load_registry(registry_file)


def sample_file(filepath: Union[str, Path]) -> str:  #!!!cmk doc
    """
    Retrieve a bed_reader sample BED file.

    Parameters
    ----------
    filepath
        Name of the sample file.

    Returns
    -------
    string
        Local name of sample file.

    Example
    --------

    .. doctest::

        >>> from bed_reader import sample_file
        >>>
        >>> file_name = sample_file("small.bed")
        >>> print(f"The local file name is '{file_name}'")
        The local file name is '...small.bed'
    """
    filepath = Path(filepath)
    file_string = str(filepath)
    if file_string.lower().endswith(".bed"):
        POOCH.fetch(file_string[:-4] + ".fam")
        POOCH.fetch(file_string[:-4] + ".bim")
    return POOCH.fetch(file_string)
