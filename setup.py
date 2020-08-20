import platform
import os
import sys
import shutil
from setuptools import setup, Extension
from distutils.command.clean import clean as Clean
import numpy
import distutils.sysconfig

# Version number
version = "0.0.2a0"


def readme():
    with open("README.md") as f:
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
    intel_root = os.path.join(os.path.dirname(__file__),"external/intel/linux")
    mp5lib = 'iomp5'
    openmp_compiler_args = ["-fopenmp","-std=c++11"]

elif "win" in platform.system().lower():
    macros = [("_WIN32", "1"),("_CRT_SECURE_NO_WARNINGS","1")]
    intel_root = os.path.join(os.path.dirname(__file__),"external/intel/windows")
    mp5lib = 'libiomp5md'
    openmp_compiler_args = ["/EHsc", "/openmp"]
else:
    macros = [("_UNIX", "1")]
    intel_root = os.path.join(os.path.dirname(__file__),"external/intel/linux")
    mp5lib = 'iomp5'
    openmp_compiler_args = ["-fopenmp","-std=c++11"]

library_list = [intel_root+"/compiler/lib/intel64"]
runtime_library_dirs = None if "win" in platform.system().lower() else library_list


# see http://stackoverflow.com/questions/4505747/how-should-i-structure-a-python-package-that-contains-cython-code
if use_cython:
    ext_modules = [
        Extension(
            name="bed_reader.wrap_plink_parser_onep",
            language="c++",
            sources=[
                "bed_reader/wrap_plink_parser_onep.pyx",
                "bed_reader/CPlinkBedFile.cpp",
            ],
            include_dirs=[numpy.get_include()],
            define_macros=macros,
        ),
        Extension(
            name="bed_reader.wrap_matrix_subset",
            language="c++",
            sources=[
                "bed_reader/wrap_matrix_subset.pyx",
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
                "bed_reader/wrap_plink_parser_openmp.pyx",
                "bed_reader/CPlinkBedFile.cpp",
            ],
            libraries = [mp5lib],
            library_dirs = library_list,
            runtime_library_dirs = runtime_library_dirs,
            include_dirs=library_list+[numpy.get_include()],
            extra_compile_args=openmp_compiler_args,
            define_macros=macros+[("OPENMP","1")]))

    cmdclass = {"build_ext": build_ext, "clean": CleanCommand}
else:
    assert False, "need codecmk"
    ext_modules = [
        Extension(
            name="bed_reader.wrap_plink_parser",
            language="c++",
            sources=[
                "bed_reader/wrap_plink_parser.cpp",
                "bed_reader/CPlinkBedFile.cpp",
            ],
            libraries = [mp5lib],
            library_dirs = library_list,
            runtime_library_dirs = runtime_library_dirs,
            include_dirs=library_list+[numpy.get_include()],
            extra_compile_args=extra_compile_args,
            define_macros=macros,
        ),
        Extension(
            name="bed_reader.wrap_matrix_subset",
            language="c++",
            sources=[
                "bed_reader/wrap_matrix_subset.cpp",
                "bed_reader/MatrixSubset.cpp",
            ],
            include_dirs=[numpy.get_include()],
            define_macros=macros,
        ),
    ]
    cmdclass = {}

install_requires = ["numpy>=1.11.3"]

# !!!cmk see FIXUP's
setup(
    name='bed-reader',
    version=version,
    description='Bed Reader',
    long_description=readme(),
    long_description_content_type = 'text/markdown',
    keywords='bioinformatics plink',
    url="https://fastlmm.github.io/",
    author='FaST-LMM Team',
    author_email='fastlmm-dev@python.org',
    license='Apache 2.0',
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python",
    ],
    packages=["bed_reader","bed_reader/tests"],  # basically everything with a __init__.py
    data_files=[("lib/site-packages/bed_reader", ["external/intel/windows/compiler/lib/intel64/libiomp5md.dll"])] if "win" in platform.system().lower() else [],
    #package_data={"bed_reader":["external\intel\windows\compiler\lib\intel64\libiomp5md.dll"]} if "win" in platform.system().lower() else {},
    install_requires=install_requires,
    # extensions
    cmdclass=cmdclass,
    ext_modules=ext_modules,
)

