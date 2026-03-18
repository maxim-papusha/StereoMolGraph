from __future__ import annotations

from collections import deque

from stereomolgraph import StereoMolGraph
from stereomolgraph.algorithms.color_refine import color_refine_smg
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms


def topological_symmetry_number(graph: StereoMolGraph) -> int:
    """
    Calculated from the number of graph isomorphisms which conserve the
    stereo information.
    symmetry_number = internal_symmetry_number * rotational_symmetry_number
    TODO: add paper reference
    """

    if any(stereo.parity is None for stereo in graph.stereo.values()):
        raise NotImplementedError(
            "all stereocenters have to be defined"
            " to calculate the symmetry number"
        )
    colorings = color_refine_smg(graph)
    mappings = vf2pp_all_isomorphisms(
        graph,
        graph,
        atom_labels=(colorings, colorings),  # type: ignore[arg-type]
        stereo=True,
    )  # type: ignore[arg-type]
    return deque(enumerate(mappings, 1), maxlen=1)[0][0]
