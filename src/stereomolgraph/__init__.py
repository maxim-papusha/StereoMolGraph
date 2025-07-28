# ruff: noqa: F401
"""Simple access of key classes MolGraph, StereoMolGraph,
CondensedReactionGraph and StereoCondensedReactionGraph"""

from __future__ import annotations

from importlib import import_module
from importlib.metadata import version

def __getattr__(name):
    """Lazy import submodules on first access (Python 3.7+)"""
    if name in {'sub1', 'sub2'}:
        return import_module(f'.{name}', __name__)
    if name == "__version__":
        return version("stereomolgraph")

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
        



