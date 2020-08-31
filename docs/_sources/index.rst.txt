################################
:mod:`bed_reader` Documentation
################################

.. currentmodule:: bed_reader

cmk example here, too?

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   

***********************
Summary
***********************

Open, Read, and Write
======================

.. autosummary::

    open_bed
    open_bed.read
    to_bed

Properties of Individuals (samples) and SNPs (variants)
==========================================================

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
==========

.. autosummary::

    sample_file
    tmp_path

***********************
Details
***********************

open_bed
==========

.. autoclass:: open_bed
   :members:
   :inherited-members:

to_bed
==========

.. autofunction:: to_bed

sample_file
==============

.. autofunction:: sample_file

tmp_path
==========

.. autofunction:: tmp_path

.. only:: html 

***********************
Indices and Tables
***********************

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`

