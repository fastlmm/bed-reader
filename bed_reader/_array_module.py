"""Resolve precision defaults for newly read arrays."""
import os

import numpy as np


def array_dtype(dtype=None, *, default="float64"):
    """Resolve explicit dtype, ARRAY_DTYPE, then the API's original default.

    ARRAY_DTYPE accepts ``float32`` or ``float64``. Explicit dtypes bypass the
    environment; the calling API validates its supported types, including int8
    where applicable. This does not change existing arrays or NumPy defaults.

    Examples
    --------
    >>> from unittest.mock import patch
    >>> from bed_reader._array_module import array_dtype
    >>> with patch.dict("os.environ", {"ARRAY_DTYPE": "float32"}):
    ...     print(array_dtype(), array_dtype("float64"))
    float32 float64
    """
    if dtype is not None:
        return np.dtype(dtype)
    configured = os.environ.get("ARRAY_DTYPE")
    if configured is None:
        return np.dtype(default)
    if configured not in ("float32", "float64"):
        raise ValueError(f"ARRAY_DTYPE must be 'float32' or 'float64', not {configured!r}")
    return np.dtype(configured)
