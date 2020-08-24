"""
Load sample data.
"""
import pooch
from pathlib import Path

POOCH = pooch.create(
    # Use the default cache folder for the OS
    path=pooch.os_cache("bed-reader"),
    # The remote data is on Github
    base_url="https://github.com/fastlmm/bed-reader/tree/master/bed_reader/tests/data/",
    # If this is a development version, get the data from the master branch
    version_dev="master",
    # The registry specifies the files that can be fetched
    env="BED_READER_DATA_DIR",
)

# Get registry file from package_data
registry_file = Path(__file__).parent / "tests/registry.txt"
# Load this registry file
POOCH.load_registry(registry_file)

def sample_file(file_name):
    return POOCH.fetch(file_name)
