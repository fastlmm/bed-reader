import distutils.sysconfig
import os
import platform
import re
import shutil
import sys
from distutils.command.clean import clean as Clean
from pathlib import Path

import numpy
from setuptools import Extension, setup


def find_version(filepath):
    import re

    version_file = read(filepath)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


def read(filepath):
    import codecs

    with codecs.open(filepath, "r") as fp:
        return fp.read()


# Version number
version = find_version(Path(__file__).parent / "bed_reader/__init__.py")


def readme():
    with open(Path(__file__).parent / "README.md") as f:
        return f.read()


try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

# use_cython=False


class CleanCommand(Clean):
    description = "Remove build directories, and compiled files (including .pyc)"

    def run(self):
        Clean.run(self)
        if os.path.exists("build"):
            shutil.rmtree("build")
        for dirpath, _, filenames in os.walk("."):
            for filename in filenames:
                if (
                    filename.endswith(".so")
                    or filename.endswith(".pyd")
                    # or filename.find("wrap_plink_parser.cpp") != -1 # remove automatically generated source file
                    # or filename.find("wrap_matrix_subset.cpp") != -1 # remove automatically generated source file
                    or filename.endswith(".pyc")
                ):
                    tmp_fn = os.path.join(dirpath, filename)
                    print("removing", tmp_fn)
                    os.unlink(tmp_fn)


# set up macro
if platform.system() == "Darwin":
    macros = [("__APPLE__", "1")]
    openmp_root = os.path.join(os.path.dirname(__file__), "external/intel/linux")
    mp5lib = "iomp5"
    openmp_compiler_args = ["-fopenmp", "-std=c++11"]
    library_list = [openmp_root + "/compiler/lib/intel64"]
    runtime_library_dirs = library_list

elif platform.system() == "Windows":
    macros = [("_WIN32", "1"), ("_CRT_SECURE_NO_WARNINGS", "1")]
    openmp_root = os.path.join(os.path.dirname(__file__), "external/llvm/windows")
    mp5lib = "libomp"
    openmp_compiler_args = ["/EHsc", "/openmp"]
    library_list = [openmp_root + "/lib"]
    runtime_library_dirs = None

else:
    macros = [("_UNIX", "1")]
    openmp_root = os.path.join(os.path.dirname(__file__), "external/intel/linux")
    mp5lib = "iomp5"
    openmp_compiler_args = ["-fopenmp", "-std=c++11"]
    library_list = [openmp_root + "/compiler/lib/intel64"]
    runtime_library_dirs = library_list


# see http://stackoverflow.com/questions/4505747/how-should-i-structure-a-python-package-that-contains-cython-code
ext_modules = [
    Extension(
        name="bed_reader.wrap_plink_parser_onep",
        language="c++",
        sources=[
            "bed_reader/wrap_plink_parser_onep.pyx"
            if use_cython
            else "bed_reader/wrap_plink_parser_onep.cpp",
            "bed_reader/CPlinkBedFile.cpp",
        ],
        include_dirs=[numpy.get_include()],
        define_macros=macros,
    ),
    Extension(
        name="bed_reader.wrap_matrix_subset",
        language="c++",
        sources=[
            "bed_reader/wrap_matrix_subset.pyx"
            if use_cython
            else "bed_reader/wrap_matrix_subset.cpp",
            "bed_reader/MatrixSubset.cpp",
        ],
        include_dirs=[numpy.get_include()],
        define_macros=macros,
    ),
]
if platform.system() != "Darwin":
    ext_modules.append(
        Extension(
            name="bed_reader.wrap_plink_parser_openmp",
            language="c++",
            sources=[
                "bed_reader/wrap_plink_parser_openmp.pyx"
                if use_cython
                else "bed_reader/wrap_plink_parser_openmp.cpp",
                "bed_reader/CPlinkBedFile.cpp",
            ],
            libraries=[mp5lib],
            library_dirs=library_list,
            runtime_library_dirs=runtime_library_dirs,
            include_dirs=library_list + [numpy.get_include()],
            extra_compile_args=openmp_compiler_args,
            define_macros=macros,
        )
    )

if use_cython:
    cmdclass = {"build_ext": build_ext, "clean": CleanCommand}
else:
    cmdclass = {}

install_requires = ["numpy>=1.11.3"]

win_data = (
    [("lib/site-packages/bed_reader", ["external/llvm/windows/bin/libomp.dll"])]
    if platform.system() == "Windows"
    else []
)
setup(
    name="bed-reader",
    version=version,
    description="Bed Reader",
    long_description=readme(),
    long_description_content_type="text/markdown",
    keywords="bioinformatics plink",
    url="https://fastlmm.github.io/",
    author="FaST-LMM Team",
    author_email="fastlmm-dev@python.org",
    license="Apache 2.0",
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python",
    ],
    packages=[
        "bed_reader",
        "bed_reader/tests",
    ],  # basically everything with a __init__.py
    data_files=win_data,
    package_data={"bed_reader/tests": ["registry.txt"],},
    install_requires=install_requires,
    # extensions
    cmdclass=cmdclass,
    ext_modules=ext_modules,
)
