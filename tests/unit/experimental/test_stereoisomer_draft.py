"""Small tests for the input-descriptor-only experimental stereoisomer draft."""

from stereomolgraph import StereoMolGraph
from stereomolgraph.experimental.stereoisomer_draft import (
    detect_stereo_elements,
    enumerate_stereo_assignments,
    locally_stereogenic_at_assignment,
)
import pytest

from stereomolgraph.stereodescriptors import SquarePlanar, Tetrahedral


def _tetrahedral_graph(ligands=("F", "Cl", "Br", "I")):
    graph = StereoMolGraph()
    graph.add_atom(0, "C")
    for atom, element in enumerate(ligands, 1):
        graph.add_atom(atom, element)
        graph.add_bond(0, atom)
    return graph


def test_only_declared_stereo_sites_are_analyzed():
    graph = _tetrahedral_graph()
    assert detect_stereo_elements(graph).sites == ()
    isomers = list(enumerate_stereo_assignments(graph))
    assert len(isomers) == 1
    assert len(isomers[0].stereo) == 0


def test_existing_tetrahedral_candidate_and_enumeration():
    graph = _tetrahedral_graph()
    graph.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), None))
    analysis = detect_stereo_elements(graph)

    assert len(analysis.sites) == 1
    assert analysis.sites[0].classification == "distinct_colors"
    assert len(analysis.sites[0].states) == 2
    assert len(list(enumerate_stereo_assignments(graph, analysis))) == 2
    assert len(list(enumerate_stereo_assignments(
        graph, analysis, max_isomers=1
    ))) == 1
    assert graph.get_atom_stereo(0).parity is None


def test_exact_local_stereogenicity_on_assigned_tetrahedra():
    graph = _tetrahedral_graph()
    graph.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), 1))
    assert locally_stereogenic_at_assignment(graph, ("atom", 0))

    symmetric = _tetrahedral_graph(ligands=("F", "F", "Br", "I"))
    symmetric.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), 1))
    assert not locally_stereogenic_at_assignment(symmetric, ("atom", 0))


def test_dependency_estimate_only_uses_supplied_sites():
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
    for a in (0, 1, 2):
        graph.set_atom_stereo(Tetrahedral((a, *sorted(graph.bonded_to(a))), None))

    analysis = detect_stereo_elements(graph)
    sites = {site.key: site for site in analysis.sites}
    assert len(sites) == 3
    assert sites[("atom", 0)].classification == "unresolved"
    assert ("atom", 1) in analysis.dependencies[("atom", 0)]
    assert ("atom", 2) in analysis.dependencies[("atom", 0)]


def test_existing_square_planar_descriptor_has_three_states():
    graph = StereoMolGraph()
    graph.add_atom(0, "Pt")
    for atom, element in enumerate(("F", "Cl", "Br", "I"), 1):
        graph.add_atom(atom, element)
        graph.add_bond(0, atom)
    graph.set_atom_stereo(SquarePlanar((0, 1, 2, 3, 4), None))

    analysis = detect_stereo_elements(graph)
    assert len(analysis.sites) == 1
    assert len(analysis.sites[0].states) == 3
    assert len(list(enumerate_stereo_assignments(graph, analysis))) == 3
    assert graph.get_atom_stereo(0).parity is None


def test_stale_analysis_is_rejected():
    graph = _tetrahedral_graph()
    graph.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), None))
    analysis = detect_stereo_elements(graph)

    # Same site key but a different ordered descriptor template.
    graph.set_atom_stereo(Tetrahedral((0, 2, 1, 3, 4), None))
    with pytest.raises(ValueError, match="analysis descriptors"):
        list(enumerate_stereo_assignments(graph, analysis))


def test_missing_ligand_slot_is_supported():
    graph = StereoMolGraph()
    graph.add_atom(0, "C")
    for idx, element in enumerate(("F", "Cl", "Br"), 1):
        graph.add_atom(idx, element)
        graph.add_bond(0, idx)
    graph.set_atom_stereo(Tetrahedral((0, 1, 2, 3, None), None))

    analysis = detect_stereo_elements(graph)
    assert len(analysis.sites) == 1
    assert len(analysis.sites[0].states) == 2
    results = list(enumerate_stereo_assignments(graph, analysis))
    assert len(results) == 2
    assert len(graph.atom_stereo) == 1
    assert graph.get_atom_stereo(0).parity is None
