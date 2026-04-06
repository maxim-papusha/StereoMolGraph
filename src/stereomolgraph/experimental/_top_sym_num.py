from collections import deque
from collections.abc import Iterable

from stereomolgraph import StereoMolGraph
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms


def _top_sym_num_from_mappings(mappings: Iterable[dict[int, int]]) -> int:
    try:
        return len(mappings)
    except Exception:
        return deque(enumerate(mappings, 1), maxlen=1)[0][0]


def topological_symmetry_number(graph: StereoMolGraph) -> int:
    """
    Calculated from the number of graph isomorphisms which conserve the
    stereo information.
    symmetry_number = internal_symmetry_number * rotational_symmetry_number
    """

    if any(stereo.parity is None for stereo in graph.stereo.values()):
        raise NotImplementedError(
            "all stereocenters have to be defined to calculate the symmetry number"
        )
    mappings = vf2pp_all_isomorphisms(
        graph,
        graph,
        stereo=True,
    )
    return _top_sym_num_from_mappings(mappings)
