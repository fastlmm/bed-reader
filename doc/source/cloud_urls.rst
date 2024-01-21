Cloud URLs Examples
====================

*Table of Contents*:

.. cmk do these links work?

- `Http <#http>`_
- `local file <#local-file>`_
- `AWS S3 <#aws-s3>`_

.. note::

   The `bed-reader` package also supports Azure and GCP, but we don't have examples.

To specify a file in the cloud, you must specify URL string plus optional cloud options.

The exact details depend on the cloud service. We'll look at `http`, at `local files`, and at `AWS S3`.

Http
----

You can read \*.bed from web sites directly. For small files, access will be fast. For medium-size files, you may need to extend the default `timeout`.

Reading from large files can also be practical and even fast under these conditions:

- You need only some of the information
- (Optional, but helpful) You can provide some metadata about individuals (samples) and SNPs (variants) locally.

Let's first look at reading a small or medium-sized dataset.

*Example:*

Read an entire file and find the fraction of missing values.

.. code-block:: python

    >>> import numpy as np
    >>> from bed_reader import open_bed
    >>> with open_bed("https://raw.githubusercontent.com/fastlmm/bed-sample-files/main/small.bed") as bed:
    ...     val = bed.read()
    ...     missing_count = np.isnan(val).sum()
    ...     missing_fraction = missing_count / val.size
    ...     missing_fraction  # doctest: +ELLIPSIS
    0.1666...

When reading a medium-sized file, you may need to set a `timeout` in your options. With a `timeout`, you can give your code more than the default 30 seconds to read metadata from the \*.fam and \*.bim files (or genomic data from \*.bed). You may also wish to use `.skip_early_check()` to avoid a fast, early check of the \*.bed file's header.

Here we print the first five iids (individual or sample ids) and first find sids (SNP or variant ids). We then, print all unique chromosome values. Finally, we read all data from chromosome 5 and print its dimensions.

.. code-block:: python

    >>> import numpy as np
    >>> from bed_reader import open_bed
    >>> with open_bed(
    ...     "https://raw.githubusercontent.com/fastlmm/bed-sample-files/main/toydata.5chrom.bed",
    ...     cloud_options={"timeout": "100s"},
    ...     skip_format_check=True,
    ...     ) as bed:
    ...     bed.iid[:5]
    ...     bed.sid[:5]
    ...     np.unique(bed.chromosome)
    ...     val = bed.read(index=np.s_[:, bed.chromosome == "5"])
    ...     val.shape
    array(['per0', 'per1', 'per2', 'per3', 'per4'], dtype='<U11')
    array(['null_0', 'null_1', 'null_2', 'null_3', 'null_4'], dtype='<U9')
    array(['1', '2', '3', '4', '5'], dtype='<U9')
    (500, 440)

Now, let's read from a large file containing data from over 1 million individuals (samples) and over 300,000 SNPs (variants). The file size is 91 GB. In this example, we read data for just one SNP (variant). If we know the number of individuals (samples) and SNPs (variants) exactly, we can read this SNP quickly and with just one file access.

What is the mean value of the SNP (variant) at index position 100,000?

.. code-block:: python

    >>> import numpy as np
    >>> from bed_reader import open_bed
    >>> with open_bed(
    ...     "https://www.ebi.ac.uk/biostudies/files/S-BSST936/genotypes/synthetic_v1_chr-10.bed",
    ...     cloud_options={"timeout": "100s"},
    ...     skip_format_check=True,
    ...     iid_count=1_008_000,
    ...     sid_count=361_561,
    ...     ) as bed:
    ...     val = bed.read(index=np.s_[:, 100_000], dtype=np.float32)
    ...     np.mean(val) # doctest: +ELLIPSIS
    0.033913...

Local File
----------

We can specify a local file as if it is in the cloud. This is a great way to test cloud functions. For real work and better efficiency, however, use `Bed` instead of `BedCloud`.

Local File URL
++++++++++++++

The URL for a local file takes the form `file:///{encoded_file_name}`. No cloud options are needed, so we use `EMPTY_OPTIONS`.

*Example:*

.. code-block:: python

    >>> import numpy as np
    >>> from bed_reader import open_bed, sample_file
    >>> from urllib.parse import urljoin
    >>> from pathlib import Path
    >>> file_name = str(sample_file("small.bed"))
    >>> print(f"file name: {file_name}")   # doctest: +ELLIPSIS
    file name: ...small.bed
    >>> url = urljoin("file:", Path(file_name).as_uri())
    >>> print(f"url: {url}") # doctest: +ELLIPSIS
    url: file:///.../small.bed
    >>> with open_bed(url) as bed:
    ...     val = bed.read(index=np.s_[:, 2], dtype=np.float64)
    ...     print(val)
    [[nan]
     [nan]
     [ 2.]]

AWS S3
------

Let's look next at reading a file (or part of a file) from AWS S3.

The URL for an AWS S3 file takes the form `s3://{bucket_name}/{s3_path}`.

AWS forbids putting some needed information in the URL. Instead, that information must go into a string-to-string map of options. Specifically, we'll put `"aws_region"`, `"aws_access_key_id"`, and `"aws_secret_access_key"` in the options. For security, we pull the last two option values from a file rather than hard-coding them into the program.

*Example:*

.. note::

   I can run this, but others can't because of the authentication checks.

.. code-block:: python

    # Your Python example here
