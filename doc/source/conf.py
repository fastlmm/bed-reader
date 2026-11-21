import os
import time
import urllib.request
from pathlib import Path

from recommonmark.parser import CommonMarkParser

doctest_global_setup = "import bed_reader"


def read(filepath):
    with open(filepath, encoding="utf-8") as fp:
        return fp.read()


def find_version(filepath):
    import re

    version_file = read(filepath)
    version_match = re.search(r"^version = ['\"]([^'\"]*)['\"]", version_file, re.MULTILINE)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


project = "Bed Reader"
copyright = "2024, Carl Kadie"
author = "Carl Kadie"

version = find_version(Path(__file__).parents[0] / Path("../../pyproject.toml"))
release = version
today = time.strftime("%B %d, %Y")

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
]

autodoc_default_flags = ["members"]
autodoc_mock_imports = ["_tkinter"]

# Wrap long signatures so trailing parameters remain visible without scrolling.
python_maximum_signature_line_length = 88
autosummary_generate = True
napoleon_numpy_docstring = True
templates_path = ["_templates"]

source_suffix = ".rst"

master_doc = "index"
man_pages = [(master_doc, project, f"{project} documentation", [author], 1)]
language = "en"

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "conf.py"]

pygments_style = "default"

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "logo_only": False,
    # "display_version": True,
    "style_external_links": True,
    "collapse_navigation": False,
}
htmlhelp_basename = f"{project}doc"


def _intersphinx_inventory(name, url):
    """Where Sphinx looks for `name`'s objects.inv.

    When INTERSPHINX_CACHE is set (as CI does), download it there first, retrying, so one slow
    docs site does not fail the -W build. Sphinx tries the local copy, then the site, and warns
    only if both fail.
    """
    cache = os.environ.get("INTERSPHINX_CACHE")
    if not cache:
        return None
    path = Path(cache) / f"{name}.inv"
    path.parent.mkdir(parents=True, exist_ok=True)
    attempts = 3
    for attempt in range(attempts):
        try:
            # numpy.org refuses urllib's default User-Agent with HTTP 403.
            request = urllib.request.Request(f"{url}objects.inv", headers={"User-Agent": "bed-reader-docs-build"})
            with urllib.request.urlopen(request, timeout=30) as response:
                path.write_bytes(response.read())
            break
        except OSError as error:
            print(f"Downloading the {name} inventory failed (attempt {attempt + 1} of {attempts}): {error}")
            if attempt + 1 < attempts:
                time.sleep(10)
    return (str(path), None) if path.exists() else None


# Without a timeout, Sphinx would wait forever on a docs site that never answers.
intersphinx_timeout = 30
intersphinx_mapping = {
    name: (url, _intersphinx_inventory(name, url))
    for name, url in {
        "python": "https://docs.python.org/3/",
        "numpy": "https://numpy.org/doc/stable/",
        "pandas": "https://pandas.pydata.org/docs/",
        "scipy": "https://docs.scipy.org/doc/scipy/",
    }.items()
}

source_parsers = {
    ".md": CommonMarkParser,
}

source_suffix = [".rst", ".md"]
