bed-reader
====================

A simple and efficient PLINK \*.bed file format reader and writer.

Read and write genetic PLINK \*.bed/bim/fam files.
Also, efficiently read parts of files and access metadata.

Documentation
=================================

* [API Documentation](http://fastlmm.github.io/bed-reader/) with examples.
* [Project Home and Full Annotated Bibliography](https://fastlmm.github.io/)

Code
=================================
* [PyPi](https://pypi.org/project/bed-reader/)
* [GitHub](https://github.com/fastlmm/bed-reader)

Contacts
=================================

* Email the developers at fastlmm-dev@python.org.
* [Join](mailto:fastlmm-user-join@python.org?subject=Subscribe) the user discussion and announcement list (or use [web sign up](https://mail.python.org/mailman3/lists/fastlmm-user.python.org)).
* [Open an issue](https://github.com/fastlmm/PySnpTools/issues) on GitHub.


Quick install:
====================

If you have pip installed, installation is as easy as:

    pip install bed-reader


Detailed Package Install Instructions:
========================================

pysnptools has the following dependencies:

python 3.7 or 3.8

Package(s):

* numpy

(1) Installing from source
-----------------------------------------

Go to the directory where you copied the source code for bed-reader

On Linux:

At the shell, type: 

    sudo python setup.py install


On Windows:

At the OS command prompt, type 

    python setup.py install



For developers (and also to run regression tests)
=========================================================

pip install -r requirements-dev.txt

When working on the developer version, first add the src directory of the package to your PYTHONPATH 
environment variable.

To build extension (from .\src dir), type the following at the OS prompt:

    python setup.py build_ext --inplace

Don't forget to set your PYTHONPATH to point to the directory above the one named bed-reader in
the bed-reader source code. For e.g. if bed-reader is in the [somedir] directory, then
in the Unix shell use:

    export PYTHONPATH=$PYTHONPATH:[somedir]

Or in the Windows DOS terminal,
one can use: 

    set PYTHONPATH=%PYTHONPATH%;[somedir]

Note for Windows: You must have Visual Studio installed.

Running regression tests
-----------------------------

From the directory tests at the top level, run:

    pytest

This will run a series of regression tests.

Note that you must use "python setup.py build_ext --inplace" to run the 
regression tests, and not "python setup.py install".
