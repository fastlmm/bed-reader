################################
:mod:`bed_reader` Documentation
################################

.. currentmodule:: bed_reader

Open, Read, and Write:

.. autosummary::

    open_bed
    open_bed.read
    to_bed

Properties of Individuals (samples) and SNPs (variants)

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

.. autosummary::

    sample_file
    tmp_path

.. autoclass:: open_bed
   :members:
   :inherited-members:

.. autofunction:: to_bed

.. autofunction:: sample_file

.. autofunction:: tmp_path

CMKTOC
