"""Focused smoke tests for the experimental stereoisomer draft."""

from stereomolgraph import StereoMolGraph
from stereomolgraph.experimental.stereoisomer_draft import (
    detect_stereo_elements,
    enumerate_stereo_assignments,
)
from stereomolgraph.stereodescriptors import Tetrahedral


def _tetrahedral_graph():
    graph = StereoMolGraph()
    for atom, element in enumerate(("C", "F", "Cl", "Br", "I")):
        graph.add_atom(atom, element)
    for ligand in (1, 2, 3, 4):
        graph.add_bond(0, ligand)
    return graph


def test_explicit_tetrahedral_candidate_and_enumeration():
    graph = _tetrahedral_graph()
    graph.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), None))
    analysis = detect_stereo_elements(graph)

    assert len(analysis.sites) == 1
    assert analysis.sites[0].classification == "distinct_colors"
    assert len(analysis.sites[0].states) == 2
    assert len(list(enumerate_stereo_assignments(graph, analysis))) == 2


def test_infer_tetrahedral_and_limit():
    graph = _tetrahedral_graph()
    assert len(detect_stereo_elements(graph).sites) == 0
    analysis = detect_stereo_elements(graph, infer_tetrahedral=True)

    assert len(analysis.sites) == 1
    assert analysis.sites[0].key == ("atom", 0)
    assert len(list(enumerate_stereo_assignments(graph, analysis, max_isomers=1))) == 1
    assert len(graph.atom_stereo) == 0  # Input is not mutated.


def test_structural_dependency_is_detected_but_not_used_for_pruning():
    graph = StereoMolGraph()
    types = {
        0: "C", 1: "C", 2: "C", 3: "F", 4: "Cl",
        5: "H", 6: "H", 7: "H", 8: "H", 9: "H", 10: "H",
    }
    for a, element in types.items():
        graph.add_atom(a, element)
    for end in (1, 2, 3, 4):
        graph.add_bond(0, end)
    for end in (5, 6, 7):
        graph.add_bond(1, end)
    for end in (8, 9, 10):
        graph.add_bond(2, end)

    analysis = detect_stereo_elements(graph, infer_tetrahedral=True)
    sites = {site.key: site for site in analysis.sites}
    assert len(sites) == 3
    assert sites[("atom", 0)].classification == "unresolved"
    assert ("atom", 1) in analysis.dependencies[("atom", 0)]
    assert ("atom", 2) in analysis.dependencies[("atom", 0)]
