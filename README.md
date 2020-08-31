bed-reader
====================

A simple and efficient PLINK .bed file format reader and writer.

Also, efficiently reads slices of data and accesses individual (sample) and SNP (variant) properties.

Examples
========

Read everything

```python
>>> import numpy as np
>>> from bed_reader import open_bed, sample_file
>>> file_name = sample_file("small.bed")
>>> bed = open_bed(file_name)
>>> val = bed.read()
>>> print(val)
[[ 1.  0. nan  0.]
 [ 2.  0. nan  2.]
 [ 0.  1.  2.  0.]]
>>> del bed


Read every 2nd individual (sample) from
SNP (variant) index position 20 (inclusive)
to 30 (exclusive).

>>> file_name2 = sample_file("some_missing.bed")
>>> bed2 = open_bed(file_name2)
>>> val2 = bed2.read(index=np.s_[::2,20:30])
>>> print(val2.shape)
(50, 10)
>>> del bed2

List the first 5 individual (sample) ids, the
first 5 SNP (variant) ids, and every unique
chromosome. Then, read every chromosome 5 value.

>>> with open_bed(file_name2) as bed3:
...     print(bed3.iid[:5])
...     print(bed3.sid[:5])
...     print(np.unique(bed3.chromosome))
...     val3 = bed3.read(index=np.s_[:,bed3.chromosome=='5'])
...     print(val3.shape)
['iid_0' 'iid_1' 'iid_2' 'iid_3' 'iid_4']
['sid_0' 'sid_1' 'sid_2' 'sid_3' 'sid_4']
['1' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '2' '20' '21' '22'
 '3' '4' '5' '6' '7' '8' '9']
(100, 6)

```

#cmk Doc test python -m doctest -v README.md

Documentation
=================================

* [Documentation](http://fastlmm.github.io/bed-reader/) with examples.
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


Detailed Package Install Instructions: #cmk need this???
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
