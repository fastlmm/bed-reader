"""Resolve the destination of newly created arrays."""
import importlib
import os


def array_module(xp=None):
    if xp is None:
        xp = os.environ.get("ARRAY_MODULE", "numpy")
    name = xp if isinstance(xp, str) else getattr(xp, "__name__", None)
    if name not in ("numpy", "cupy"):
        raise ValueError(f"Unknown array module: {xp!r}")
    return importlib.import_module(name)
