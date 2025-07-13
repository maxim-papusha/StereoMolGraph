# ruff: noqa: F401
"""Simple access of key classes MolGraph, StereoMolGraph,
CondensedReactionGraph and StereoCondensedReactionGraph"""

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from stereomolgraph.graphs.mg import AtomId, Bond, MolGraph
    from stereomolgraph.graphs.smg import StereoMolGraph
    from stereomolgraph.graphs.crg import CondensedReactionGraph
    from stereomolgraph.graphs.scrg import StereoCondensedReactionGraph
    from stereomolgraph.periodictable import (
        Element,
        PERIODIC_TABLE,
        COVALENT_RADII,
    )


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

        case "Element":
            from stereomolgraph.periodictable import Element

            return Element
        case "PERIODIC_TABLE":
            from stereomolgraph.periodictable import PERIODIC_TABLE

            return PERIODIC_TABLE
        case "COVALENT_RADII":
            from stereomolgraph.periodictable import COVALENT_RADII

            return COVALENT_RADII
        case "__version__":
            from importlib.metadata import version

            return version("stereomolgraph")

        case _:
            raise AttributeError(f"module has no attribute {name}")
