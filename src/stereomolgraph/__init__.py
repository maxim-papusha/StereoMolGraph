"""Simple access of key classes MolGraph, StereoMolGraph,
CondensedReactionGraph and StereoCondensedReactionGraph"""

__all__ = [
    "MolGraph",
    "StereoMolGraph",
    "CondensedReactionGraph",
    "StereoCondensedReactionGraph",
    "Element",
    "PERIODIC_TABLE",
    "COVALENT_RADII",
]

def __getattr__(name):
    match name:
        case "MolGraph":
            from stereomolgraph.graph import MolGraph
            return MolGraph
        case "StereoMolGraph":
            from stereomolgraph.graph import StereoMolGraph
            return StereoMolGraph
        case "CondensedReactionGraph":
            from stereomolgraph.graph import CondensedReactionGraph
            return CondensedReactionGraph
        case "StereoCondensedReactionGraph":
            from stereomolgraph.graph import StereoCondensedReactionGraph
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