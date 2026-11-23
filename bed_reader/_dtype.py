"""Resolve precision defaults for newly read arrays."""
import os

import numpy as np


def array_dtype(dtype, *, default):
    """Resolve an explicit dtype, then ARRAY_DTYPE, then the calling API's default.

    ARRAY_DTYPE accepts ``float32`` or ``float64``. An explicit dtype bypasses the
    environment and must be ``int8``, ``float32``, or ``float64``; anything else
    raises ``ValueError`` here, before any read path is chosen. This does not
    change existing arrays or NumPy defaults.

    Examples
    --------
    >>> from unittest.mock import patch
    >>> from bed_reader._dtype import array_dtype
    >>> with patch.dict("os.environ", {"ARRAY_DTYPE": "float32"}):
    ...     print(array_dtype(None, default="float32"), array_dtype("float64", default="float32"))
    float32 float64
    """
    if dtype is not None:
        dtype = np.dtype(dtype)
        # Validate here so the Python-only reader and reads of empty files reject
        # the same dtypes as the native reader (see open_bed._pick_reader).
        if dtype not in (np.int8, np.float32, np.float64):
            raise ValueError(
                f"dtype '{dtype}' not known, only "
                + "'int8', 'float32', and 'float64' are allowed.",
            )
        return dtype
    configured = os.environ.get("ARRAY_DTYPE")
    if configured is None:
        return np.dtype(default)
    if configured not in ("float32", "float64"):
        raise ValueError(f"ARRAY_DTYPE must be 'float32' or 'float64', not {configured!r}")
    return np.dtype(configured)
