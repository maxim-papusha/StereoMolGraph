# ruff: noqa: F401
"""Simple access of key classes MolGraph, StereoMolGraph,
CondensedReactionGraph and StereoCondensedReactionGraph"""

from __future__ import annotations

from stereomolgraph.graphs.crg import CondensedReactionGraph
from stereomolgraph.graphs.mg import AtomId, Bond, MolGraph
from stereomolgraph.graphs.scrg import StereoCondensedReactionGraph
from stereomolgraph.graphs.smg import StereoMolGraph


def __getattr__(name: str):
    match name:
        case "__version__":
            from importlib.metadata import version

            return version("stereomolgraph")

        case name:
            from importlib import import_module

            import_module(f".{name}", __name__)
