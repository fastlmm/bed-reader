.. Comment: This must be kept in sync with the README.md file manually.

.. image:: https://badge.fury.io/py/bed-reader.svg
    :target: https://badge.fury.io/py/bed-reader
.. image:: https://github.com/fastlmm/bed-reader/actions/workflows/ci.yml/badge.svg?branch=master
    :target: https://github.com/fastlmm/bed-reader/actions/workflows/ci.yml
.. image:: https://img.shields.io/pypi/pyversions/bed-reader

################################
:mod:`bed_reader` Documentation
################################

.. currentmodule:: bed_reader

Read and write the PLINK BED format, simply and efficiently.

Features:

* Fast multi-threaded Rust engine.
* Supports all Python indexing methods. Slice data by individuals (samples) and/or SNPs (variants).
* Used by `PySnpTools <https://github.com/fastlmm/PySnpTools>`_, `FaST-LMM <https://github.com/fastlmm/FaST-LMM>`_, and `PyStatGen <https://github.com/pystatgen>`_.
* Supports `PLINK 1.9 <https://www.cog-genomics.org/plink2/formats>`_.
* Read data locally or from the cloud, efficiently and directly.

Install
====================

**Full version**: With all optional dependencies:

.. code-block:: bash

    pip install bed-reader[samples,sparse]

**Minimal version**: Depends only on `numpy`:

.. code-block:: bash

    pip install bed-reader

Usage
========

Read genotype data from a .bed file.

::

    >>> import numpy as np
    >>> from bed_reader import open_bed, sample_file
    >>>
    >>> file_name = sample_file("small.bed")
    >>> bed = open_bed(file_name)
    >>> val = bed.read()
    >>> print(val)
    [[ 1.  0. nan  0.]
     [ 2.  0. nan  2.]
     [ 0.  1.  2.  0.]]
    >>> del bed

Read every second individual and SNPs (variants) from 20 to 30.

::

    >>> file_name2 = sample_file("some_missing.bed")
    >>> bed2 = open_bed(file_name2)
    >>> val2 = bed2.read(index=np.s_[::2,20:30])
    >>> print(val2.shape)
    (50, 10)
    >>> del bed2

List the first 5 individual (sample) ids, the
first 5 SNP (variant) ids, and every unique
chromosome. Then, read every value in chromosome 5.

::

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

From the cloud: open a file and read data for one SNP (variant)
at index position 2. (See :doc:`cloud_urls` for details on cloud URLs.)

::

    >>> with open_bed("https://raw.githubusercontent.com/fastlmm/bed-sample-files/main/small.bed") as bed:
    ...     val = bed.read(index=np.s_[:,2], dtype="float64")
    ...     print(val)
    [[nan]
     [nan]
     [ 2.]]


Project Links
==============

- `Documentation <http://fastlmm.github.io/bed-reader>`_
- Questions to `fastlmm-dev@python.org <mailto:fastlmm-dev@python.org>`_
- `Source code <https://github.com/fastlmm/bed-reader>`_
- `PyPI <https://pypi.org/project/bed-reader>`_
- `Bug reports <https://github.com/fastlmm/bed-reader/issues>`_
- `Mailing list <https://mail.python.org/mailman3/lists/fastlmm-user.python.org>`_
- `Project Website <https://fastlmm.github.io/>`_
- `Change Log <https://github.com/fastlmm/bed-reader/blob/master/CHANGELOG.md>`_


Summary
======================

Open, Read, and Write
----------------------

.. autosummary::

    open_bed
    open_bed.read
    open_bed.read_sparse
    to_bed

Properties of Individuals (samples) and SNPs (variants)
------------------------------------------------------------------

.. autosummary::

    open_bed.iid_count
    open_bed.sid_count
    open_bed.shape

    open_bed.fid
    open_bed.iid
    open_bed.father
    open_bed.mother
    open_bed.sex
    open_bed.pheno
    open_bed.sid
    open_bed.chromosome
    open_bed.cm_position
    open_bed.bp_position
    open_bed.allele_1
    open_bed.allele_2

    open_bed.properties
    open_bed.property_item

Utilities
----------------------

.. autosummary::

    sample_file
    tmp_path

Details
========

open_bed
----------

.. autoclass:: open_bed
   :members:
   :inherited-members:

to_bed
----------

.. autofunction:: to_bed

sample_file
-------------

.. autofunction:: sample_file

tmp_path
-------------

.. autofunction:: tmp_path

.. only:: html 


Environment Variables
======================

By default :meth:`sample_file` puts files under the user's cache directory. Override this by setting
the ``BED_READER_DATA_DIR`` environment variable.

By default, :class:`open_bed` uses all available processors. Override this with the ``num_threads``
parameter or by setting environment variable (listed in priority order):
'PST_NUM_THREAD

Cloud URL Examples
======================

.. toctree::
   :maxdepth: 1

   cloud_urls 

Indices and Tables
====================

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`