PySnpTools#!!!cmk update
====================

PySnpTools is a library for reading and manipulating genetic data.

Main Features:

* [SnpReader](http://fastlmm.github.io/PySnpTools): Efficiently read genetic PLINK formats including \*.bed/bim/fam files.
          Also, efficiently read parts of files, read kernel data, and standardize data. 
          New features include on-the-fly SNP generation, larger in-memory data, and
          cluster-ready BED data.

* [DistReader](https://fastlmm.github.io/PySnpTools/#module-pysnptools.distreader): Efficiently work with
         unphased BGEN format and other diploid, biallelic distribution data.
          Also, efficiently read parts of files. See [Distribution IPython Notebook](https://nbviewer.jupyter.org/github/fastlmm/PySnpTools/blob/master/doc/ipynb/Dist.ipynb).

* [util](https://fastlmm.github.io/PySnpTools/#module-pysnptools.util): In one line, intersect and re-order IIDs from snpreader and other sources.
          Also, efficiently extract a submatrix from an ndarray. 

* [IntRangeSet](https://fastlmm.github.io/PySnpTools/#util-intrangeset): Efficiently manipulate ranges of integers - for example, genetic position - with set operators including union, intersection, and set difference. 

* [mapreduce1](https://fastlmm.github.io/PySnpTools/#module-pysnptools.util.mapreduce1): Run loops locally, on multiple processors, or on any cluster. 

* [filecache](https://fastlmm.github.io/PySnpTools/#module-pysnptools.util.filecache):  Read and write files locally or from/to any remote storage.

Documentation
=================================

* [API Documentation](http://fastlmm.github.io/PySnpTools/) with examples. It includes links to tutorial slides, notebooks, and video.
* [Project Home and Full Annotated Bibliography](https://fastlmm.github.io/)

Code
=================================
* [PyPi](https://pypi.org/project/pysnptools/)
* [GitHub](https://github.com/fastlmm/PySnpTools)

Contacts
=================================

* Email the developers at fastlmm-dev@python.org.
* [Join](mailto:fastlmm-user-join@python.org?subject=Subscribe) the user discussion and announcement list (or use [web sign up](https://mail.python.org/mailman3/lists/fastlmm-user.python.org)).
* [Open an issue](https://github.com/fastlmm/PySnpTools/issues) on GitHub.


Quick install:
====================

If you have pip installed, installation is as easy as:

    pip install pysnptools


Detailed Package Install Instructions:
========================================

pysnptools has the following dependencies:

python 3.7 or 3.8 (on MacOS, just 3.7)

Packages:

* numpy
* scipy
* pandas
* cython
* psutil
* h5py
* dill


(1) Installation of dependent packages
-----------------------------------------

We recommend using a Python distribution such as 
[Anaconda](https://www.anaconda.com/distribution/).
This distribution can be used on Linux and Windows and is free.
It is the easiest way to get all the required package
dependencies.


(2) Installing from source
-----------------------------------------

Go to the directory where you copied the source code for pysnptools.

On Linux:

At the shell, type: 

    sudo python setup.py install


On Windows:

At the OS command prompt, type 

    python setup.py install



For developers (and also to run regression tests)
=========================================================

When working on the developer version, first add the src directory of the package to your PYTHONPATH 
environment variable.

For building C-extensions, first make sure all of the above dependencies are installed (including cython)

To build extension (from .\src dir), type the following at the OS prompt:

    python setup.py build_ext --inplace

Note, if this fails with a gcc permission denied error, then specifying the correct compiler will likely fix the problem, e.g.

    python setup.py build_ext --inplace --compiler=msvc


Don't forget to set your PYTHONPATH to point to the directory above the one named pysnptools in
the pysnptools source code. For e.g. if pysnptools is in the [somedir] directory, then
in the unix shell use:

    export PYTHONPATH=$PYTHONPATH:[somedir]

Or in the Windows DOS terminal,
one can use: 

    set PYTHONPATH=%PYTHONPATH%;[somedir]

(or use the Windows GUI for env variables).

Note for Windows: You must have Visual Studio installed. If you have VisualStudio2008 installed 
(which was used to build python2.7) you need to nothing more. Otherwise, follow these instructions:

If you have Visual Studio 2015 installed:

    SET VS90COMNTOOLS=%VS130COMNTOOLS%

or with Visual Studio 2017 installed:

    SET VS90COMNTOOLS=%VS140COMNTOOLS%

or with Visual Studio 2019 installed:

    SET VS90COMNTOOLS=%VS150COMNTOOLS%

Running regression tests
-----------------------------

From the directory tests at the top level, run:

    python test.py

This will run a
series of regression tests, reporting "." for each one that passes, "F" for each
one that does not match up, and "E" for any which produce a run-time error. After
they have all run, you should see the string "............" indicating that they 
all passed, or if they did not, something such as "....F...E......", after which
you can see the specific errors.

Note that you must use "python setup.py build_ext --inplace" to run the 
regression tests, and not "python setup.py install".
