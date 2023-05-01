import time
from pathlib import Path


def read(filepath):
    import codecs

    with codecs.open(filepath, "r") as fp:
        return fp.read()


def find_version(filepath):
    import re

    version_file = read(filepath)
    version_match = re.search(r"^version = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


project = "Bed Reader"
copyright = "2023, Carl Kadie"
author = "Carl Kadie"

version = find_version(Path(__file__).parents[0] / Path("../../Cargo.toml"))
release = version
today = time.strftime("%B %d, %Y")

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.autosectionlabel",
]

autodoc_default_flags = ["members"]
autodoc_mock_imports = ["_tkinter"]
autosummary_generate = True
napoleon_numpy_docstring = True
templates_path = ["_templates"]
autosectionlabel_prefix_document = True

source_suffix = ".rst"

master_doc = "index"
man_pages = [(master_doc, project, "{} documentation".format(project), [author], 1)]
language = 'en'

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "conf.py"]

pygments_style = "default"

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "logo_only": False,
    "display_version": True,
    "style_external_links": True,
    "collapse_navigation": False,
}
htmlhelp_basename = "{}doc".format(project)

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
}
