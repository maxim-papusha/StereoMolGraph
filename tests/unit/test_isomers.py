from stereomolgraph import StereoCondensedReactionGraph, StereoMolGraph
from stereomolgraph.experimental import (
    generate_fleeting_stereoisomers,
    generate_stereoisomers,
)
from stereomolgraph.graphs.crg import Change
from stereomolgraph.stereodescriptors import Tetrahedral


def _undefined_tetrahedral_graph() -> StereoMolGraph:
    graph = StereoMolGraph()
    graph.add_atom(0, "C")
    graph.add_atom(1, "Cl")
    graph.add_atom(2, "Br")
    graph.add_atom(3, "F")
    graph.add_atom(4, "I")

    for atom in (1, 2, 3, 4):
        graph.add_bond(0, atom)

    graph.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), None))
    return graph


def _undefined_fleeting_tetrahedral_graph() -> StereoCondensedReactionGraph:
    graph = StereoCondensedReactionGraph()
    graph.add_atom(0, "C")
    graph.add_atom(1, "Cl")
    graph.add_atom(2, "Br")
    graph.add_atom(3, "F")
    graph.add_atom(4, "I")

    for atom in (1, 2, 3, 4):
        graph.add_bond(0, atom)

    graph.set_atom_stereo_change(fleeting=Tetrahedral((0, 1, 2, 3, 4), None))
    return graph


def test_generate_stereoisomers_enantiomers_false_no_duplicates():
    graph = _undefined_tetrahedral_graph()

    stereoisomers = list(generate_stereoisomers(graph, enantiomers=False))
    unique_isomers = {g.copy(frozen=True) for g in stereoisomers}

    assert len(stereoisomers) == 2
    assert len(unique_isomers) == 2
    assert {g.get_atom_stereo(0).parity for g in stereoisomers} == {1, -1}


def test_generate_stereoisomers_enantiomers_true_one_representative():
    graph = _undefined_tetrahedral_graph()

    stereoisomers = list(generate_stereoisomers(graph, enantiomers=True))
    unique_isomers = {g.copy(frozen=True) for g in stereoisomers}

    assert len(stereoisomers) == 1
    assert len(unique_isomers) == 1


def test_generate_fleeting_stereoisomers_enantiomers_true_no_duplicates():
    graph = _undefined_fleeting_tetrahedral_graph()

    stereoisomers = list(generate_fleeting_stereoisomers(graph, enantiomers=True))
    unique_isomers = {g.copy(frozen=True) for g in stereoisomers}

    assert len(stereoisomers) == 2
    assert len(unique_isomers) == 2
    assert {
        g.get_atom_stereo_change(0)[Change.FLEETING].parity for g in stereoisomers
    } == {1, -1}


def test_generate_fleeting_stereoisomers_enantiomers_false_one_representative():
    graph = _undefined_fleeting_tetrahedral_graph()

    stereoisomers = list(generate_fleeting_stereoisomers(graph, enantiomers=False))
    unique_isomers = {g.copy(frozen=True) for g in stereoisomers}

    assert len(stereoisomers) == 1
    assert len(unique_isomers) == 1
