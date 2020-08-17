################################
:mod:`pysnptools` Documentation
################################

PySnpTools: A library for reading and manipulating genetic data.

See `PySnpTools's README.md <https://github.com/fastlmm/PySnpTools/blob/master/README.md>`_ for installation instructions, documentation, and code.

:Synopsis:

* :mod:`.snpreader`: Efficiently read genetic PLINK formats including \*.bed/bim/fam and phenotype files. Also, efficiently read *parts* of files and standardize data.

* :new: :mod:`.distreader`: Efficiently work with unphased BGEN format and other diploid, biallelic distribution data. Also, efficiently read *parts* of files.
    See `Distribution IPython Notebook <https://nbviewer.jupyter.org/github/fastlmm/PySnpTools/blob/master/doc/ipynb/Dist.ipynb>`_

* :mod:`.kernelreader`: Efficiently create, read, and manipulate kernel data.

* :mod:`.standardizer`: Specify standardizers for :mod:`.snpreader`.

* :mod:`.kernelstandardizer`: Specify standardizers for :mod:`.kernelreader`.

* :mod:`.pstreader`: Generalizes :mod:`.snpreader` and :mod:`.kernelreader` (provides the efficiency of numpy arrays with some of the flexibility of pandas)

* :mod:`.util`: In one line, intersect and re-order IIDs from :mod:`.snpreader`, :mod:`.kernelreader` and other sources. Also, efficiently extract a submatrix from an ndarray.

* :class:`.util.IntRangeSet`: Efficiently manipulate ranges of integers -- for example, genetic position -- with set operators including
  union, intersection, and set difference. 


* :mod:`.util.mapreduce1`: Run loops in parallel on multiple processes, threads, or clusters.

* :mod:`.util.filecache`: Automatically copy files to and from any remote storage.


:Tutorial:

* `White Paper <https://fastlmm.github.io/KadiePySnpTools.pdf>`_

*From PyData Conference, Seattle*

* `Slides <https://1drv.ms/p/s!AkoPP4cC5J64jOc15JXgDsXmD9EEaw>`_
* `Video <https://www.youtube.com/watch?v=KPI6479ctAQ>`_
* IPython Notebooks
   * `General Notebook <https://nbviewer.jupyter.org/github/fastlmm/PySnpTools/blob/master/doc/ipynb/tutorial.ipynb>`_
   * `Distribution Notebook <https://nbviewer.jupyter.org/github/fastlmm/PySnpTools/blob/master/doc/ipynb/Dist.ipynb>`_

:Code:

* `PyPi <https://pypi.org/project/pysnptools/>`_
* `GitHub <https://github.com/fastlmm/PySnpTools>`_

:Contacts:

* Email the developers at `fastlmm-dev@python.org <mailto:fastlmm-dev@python.org>`_.
* `Join <mailto:fastlmm-user-join@python.org?subject=Subscribe>`_ the user discussion
  and announcements list (or use `web sign up <https://mail.python.org/mailman3/lists/fastlmm-user.python.org>`_).
* `Open an issue <https://github.com/fastlmm/PySnpTools/issues>`_ on GitHub.
* `Project Home <https://fastlmm.github.io/>`_.



.. automodule:: pysnptools
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:

***********************
:mod:`snpreader` Module
***********************

.. automodule:: pysnptools.snpreader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:


:class:`snpreader.SnpReader`
+++++++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpReader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs


:class:`snpreader.Bed`
++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Bed
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row

:class:`snpreader.SnpData`
++++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpData
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row

:class:`snpreader.Pheno`
++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Pheno
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row

:class:`snpreader.Ped`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Ped
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row

:class:`snpreader.Dat`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Dat
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row

:class:`snpreader.Dense`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Dense
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row

:class:`snpreader.SnpHdf5`
+++++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpHdf5
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property

:class:`snpreader.SnpNpz`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpNpz
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


:class:`snpreader.SnpMemMap`
++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpMemMap
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property

:class:`snpreader.SnpGen`
++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpGen
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property

:class:`snpreader.DistributedBed`
++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.DistributedBed
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


****************************
:mod:`kernelreader` Module
****************************
.. automodule:: pysnptools.kernelreader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:


:class:`kernelreader.KernelReader`
+++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.KernelReader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs, col_property, row_property


:class:`kernelreader.KernelData`
++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.KernelData
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


:class:`kernelreader.KernelNpz`
++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.KernelNpz
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property

:class:`kernelreader.KernelHdf5`
++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.KernelHdf5
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


:class:`kernelreader.Identity`
+++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.Identity
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property

:class:`kernelreader.SnpKernel`
+++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.SnpKernel
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property




***************************
:mod:`standardizer` Module
***************************

.. automodule:: pysnptools.standardizer
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:


:class:`standardizer.Standardizer`
+++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.Standardizer
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs

:class:`standardizer.Unit`
+++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.Unit
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize

:class:`standardizer.Beta`
+++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.Beta
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize

:class:`standardizer.Identity`
+++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.Identity
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize,is_constant

:class:`standardizer.DiagKtoN`
+++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.DiagKtoN
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize

:class:`standardizer.UnitTrained`
+++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.UnitTrained
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize,is_constant

:class:`standardizer.DiagKtoNTrained`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.DiagKtoNTrained
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize,is_constant


******************************************************
:mod:`kernelstandardizer` Module
******************************************************

:class:`kernelstandardizer.KernelStandardizer`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelstandardizer.KernelStandardizer
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs

:class:`kernelstandardizer.DiagKtoN`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelstandardizer.DiagKtoN
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize

:class:`kernelstandardizer.DiagKtoNTrained`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelstandardizer.DiagKtoNTrained
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize,is_constant

:class:`kernelstandardizer.Identity`
+++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelstandardizer.Identity
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize

**************************
:mod:`distreader` Module
**************************

.. automodule:: pysnptools.distreader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:


:class:`distreader.DistReader`
++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.distreader.DistReader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs


:class:`distreader.Bgen`
+++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.distreader.Bgen
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row

.. autofunction:: pysnptools.distreader.bgen.default_iid_function

.. autofunction:: pysnptools.distreader.bgen.default_sid_function

.. autofunction:: pysnptools.distreader.bgen.default_sample_function

.. autofunction:: pysnptools.distreader.bgen.default_id_rsid_function

:class:`distreader.DistData`
+++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.distreader.DistData
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row


:class:`distreader.DistHdf5`
++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.distreader.DistHdf5
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property

:class:`distreader.DistNpz`
++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.distreader.DistNpz
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


:class:`distreader.DistMemMap`
++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.distreader.DistMemMap
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property

:class:`distreader.DistGen`
++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.distreader.DistGen
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


***********************
:mod:`pstreader` Module
***********************

.. automodule:: pysnptools.pstreader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:


:class:`pstreader.PstReader`
+++++++++++++++++++++++++++++
.. autoclass:: pysnptools.pstreader.PstReader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs


:class:`pstreader.PstData`
++++++++++++++++++++++++++
.. autoclass:: pysnptools.pstreader.PstData
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


:class:`pstreader.PstHdf5`
+++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.pstreader.PstHdf5
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


:class:`pstreader.PstNpz`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.pstreader.PstNpz
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property

:class:`pstreader.PstMemMap`
++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.pstreader.PstMemMap
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


***********************
:mod:`util` Module
***********************

:mod:`util`
++++++++++++++++++++++++++
.. automodule:: pysnptools.util
    :members:
    :undoc-members:
	:show-inheritance:

.. autofunction:: pysnptools.util.example_file

.. autofunction:: pysnptools.util.snp_gen

:class:`util.IntRangeSet`
++++++++++++++++++++++++++
.. autoclass:: pysnptools.util.IntRangeSet
  :members:
  :undoc-members:
  :show-inheritance:
  :special-members:
  :exclude-members: __and__, __weakref__,__module__,__dict__, __add__


:mod:`util.pheno`
+++++++++++++++++++++++++++++++++++++++++++++
.. automodule:: pysnptools.util.pheno
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:

*******************************************
:mod:`util.mapreduce1` Module
*******************************************
.. automodule:: pysnptools.util.mapreduce1
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:


:mod:`util.mapreduce1`
++++++++++++++++++++++++++

.. autofunction:: pysnptools.util.mapreduce1.map_reduce

:class:`util.mapreduce1.runner.Runner`
++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.util.mapreduce1.runner.Runner
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: run

:class:`util.mapreduce1.runner.Local`
++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.util.mapreduce1.runner.Local
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: run

:class:`util.mapreduce1.runner.LocalMultiProc`
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.util.mapreduce1.runner.LocalMultiProc
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: run

:class:`util.mapreduce1.runner.LocalMultiThread`
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.util.mapreduce1.runner.LocalMultiThread
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: run

:class:`util.mapreduce1.runner.LocalInParts`
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.util.mapreduce1.runner.LocalInParts
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: run


*******************************************
:mod:`util.filecache` Module
*******************************************
.. automodule:: pysnptools.util.filecache
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:

:class:`util.filecache.FileCache`
++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.util.filecache.FileCache
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members:

:class:`util.filecache.LocalCache`
++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.util.filecache.LocalCache
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members:

:class:`util.filecache.PeerToPeer`
++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.util.filecache.PeerToPeer
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members:

:class:`util.filecache.Hashdown`
++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.util.filecache.Hashdown
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members:


.. only:: html 

***********************
Indices and Tables
***********************

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
