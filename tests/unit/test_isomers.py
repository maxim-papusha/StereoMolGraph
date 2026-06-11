from __future__ import annotations

from stereomolgraph import StereoMolGraph
from stereomolgraph.experimental import generate_stereoisomers
from stereomolgraph.stereodescriptors import PlanarBond, Tetrahedral


def _single_tetrahedral_graph() -> StereoMolGraph:
    graph = StereoMolGraph()
    graph.add_atom(0, "C")
    graph.add_atom(1, "H")
    graph.add_atom(2, "F")
    graph.add_atom(3, "Cl")
    graph.add_atom(4, "Br")
    graph.add_bond(0, 1)
    graph.add_bond(0, 2)
    graph.add_bond(0, 3)
    graph.add_bond(0, 4)
    graph.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), None))
    return graph


def _meso_candidate_graph() -> StereoMolGraph:
    graph = StereoMolGraph()
    graph.add_atom(0, "C")
    graph.add_atom(1, "C")
    graph.add_atom(2, "H")
    graph.add_atom(3, "F")
    graph.add_atom(4, "Cl")
    graph.add_atom(5, "H")
    graph.add_atom(6, "F")
    graph.add_atom(7, "Cl")

    graph.add_bond(0, 1)
    graph.add_bond(0, 2)
    graph.add_bond(0, 3)
    graph.add_bond(0, 4)
    graph.add_bond(1, 5)
    graph.add_bond(1, 6)
    graph.add_bond(1, 7)

    graph.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), None))
    graph.set_atom_stereo(Tetrahedral((1, 0, 5, 6, 7), None))
    return graph


def _planar_bond_graph() -> StereoMolGraph:
    graph = StereoMolGraph()
    graph.add_atom(0, "F")
    graph.add_atom(1, "H")
    graph.add_atom(2, "C")
    graph.add_atom(3, "C")
    graph.add_atom(4, "F")
    graph.add_atom(5, "H")

    graph.add_bond(0, 2)
    graph.add_bond(1, 2)
    graph.add_bond(2, 3)
    graph.add_bond(3, 4)
    graph.add_bond(3, 5)

    graph.set_bond_stereo(PlanarBond((0, 1, 2, 3, 4, 5), None))
    return graph


def _generated_set(
    graph: StereoMolGraph,
    strategy: str,
    **kwargs,
) -> set[StereoMolGraph]:
    return set(generate_stereoisomers(graph, strategy=strategy, **kwargs))


def test_staged_matches_bruteforce_for_single_tetrahedral_center():
    graph = _single_tetrahedral_graph()

    brute_force = _generated_set(graph, "bruteforce")
    staged = _generated_set(graph, "staged")

    assert staged == brute_force
    assert len(staged) == 2


def test_staged_matches_bruteforce_for_meso_candidate():
    graph = _meso_candidate_graph()

    brute_force = _generated_set(graph, "bruteforce")
    staged = _generated_set(graph, "staged")

    assert staged == brute_force
    assert len(staged) == 3


def test_staged_can_collapse_enantiomers_after_staged_enumeration():
    graph = _meso_candidate_graph()

    staged = list(generate_stereoisomers(graph, strategy="staged", enantiomers=False))

    assert len(staged) == 2


def test_staged_matches_bruteforce_for_planar_bond():
    graph = _planar_bond_graph()

    brute_force = _generated_set(graph, "bruteforce")
    staged = _generated_set(graph, "staged")

    assert staged == brute_force
    assert len(staged) == 2
