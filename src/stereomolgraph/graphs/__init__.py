# ruff: noqa: F401
"""Simple access of key classes MolGraph, StereoMolGraph,
CondensedReactionGraph and StereoCondensedReactionGraph"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from stereomolgraph.graphs.crg import CondensedReactionGraph
    from stereomolgraph.graphs.mg import AtomId, Bond, MolGraph
    from stereomolgraph.graphs.scrg import StereoCondensedReactionGraph
    from stereomolgraph.graphs.smg import StereoMolGraph


def __getattr__(name: str):
    match name:
        case "AtomId":
            from stereomolgraph.graphs.mg import AtomId

            return AtomId
        case "Bond":
            from stereomolgraph.graphs.mg import Bond

            return Bond
        case "MolGraph":
            from stereomolgraph.graphs.mg import MolGraph

            return MolGraph
        case "StereoMolGraph":
            from stereomolgraph.graphs.smg import StereoMolGraph

            return StereoMolGraph
        case "CondensedReactionGraph":
            from stereomolgraph.graphs.crg import CondensedReactionGraph

            return CondensedReactionGraph
        case "StereoCondensedReactionGraph":
            from stereomolgraph.graphs.scrg import StereoCondensedReactionGraph

            return StereoCondensedReactionGraph

        case _:
            raise AttributeError(f"module has no attribute {name}")
