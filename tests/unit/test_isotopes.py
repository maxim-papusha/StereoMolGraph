"""Regression tests for isotope-aware molecular graphs."""
import numpy as np
import pytest
from rdkit import Chem

from stereomolgraph import (
    CondensedReactionGraph,
    MolGraph,
    StereoCondensedReactionGraph,
    StereoMolGraph,
)
from stereomolgraph.algorithms.circular import label_hash
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.graph2rdmol import mol_graph_to_rdmol
from stereomolgraph.stereodescriptors import Tetrahedral


@pytest.mark.parametrize(
    "graph_type",
    [MolGraph, StereoMolGraph, CondensedReactionGraph, StereoCondensedReactionGraph],
)
def test_isotope_aware_equality_and_hash(graph_type):
    labeled = graph_type()
    labeled.add_atom(0, "C", isotope=13)
    labeled.add_atom(1, "H", isotope=2)
    labeled.add_bond(0, 1)

    renumbered = graph_type()
    renumbered.add_atom(11, "H", isotope=2)
    renumbered.add_atom(12, "C", isotope=13)
    renumbered.add_bond(12, 11)
    assert labeled.isotopes == (13, 2)
    assert labeled == renumbered
    assert hash(labeled.freeze()) == hash(renumbered.freeze())

    ordinary = graph_type()
    ordinary.add_atom(0, "C")
    ordinary.add_atom(1, "H")
    ordinary.add_bond(0, 1)
    assert labeled != ordinary


def test_unspecified_isotope_and_validation():
    graph = MolGraph()
    graph.add_atom(0, "C")
    assert graph.get_atom_isotope(0) == 0
    assert "isotope" not in graph.get_atom_attributes(0)

    graph.set_atom_attribute(0, "isotope", 13)
    assert graph.get_atom_isotope(0) == 13
    graph.set_atom_attribute(0, "isotope", 0)
    assert graph.get_atom_isotope(0) == 0
    assert "isotope" not in graph.get_atom_attributes(0)

    for invalid in (-1, None, True, 1.5, "13"):
        with pytest.raises(ValueError):
            graph.set_atom_attribute(0, "isotope", invalid)


def test_initial_isotope_labels():
    graph = MolGraph()
    graph.add_atom(0, "C")
    graph.add_atom(1, "C", isotope=12)
    graph.add_atom(2, "C", isotope=13)
    assert label_hash(graph).tolist() == [6, 6 + 12 * 256, 6 + 13 * 256]
    assert label_hash(graph, ("atom_type",)).tolist() == [6, 6, 6]


def test_isomorphism_requires_exact_isotopes_with_custom_colors():
    a = MolGraph()
    b = MolGraph()
    a.add_atom(0, "H", isotope=1)
    b.add_atom(0, "H", isotope=2)
    labels = (np.zeros(1, dtype=np.int64), np.zeros(1, dtype=np.int64))
    assert not list(vf2pp_all_isomorphisms(a, b, atom_labels=labels))


def test_methane_deuteration_reduces_automorphisms():
    graph = MolGraph()
    graph.add_atom(0, "C")
    for atom in range(1, 5):
        graph.add_atom(atom, "H")
        graph.add_bond(0, atom)
    assert len(list(vf2pp_all_isomorphisms(graph, graph))) == 24
    graph.set_atom_attribute(1, "isotope", 2)
    assert len(list(vf2pp_all_isomorphisms(graph, graph))) == 6


def test_isotopic_tetrahedral_enantiomers():
    graph = StereoMolGraph()
    graph.add_atom(0, "C")
    for atom, element, isotope in (
        (1, "H", 1),
        (2, "H", 2),
        (3, "H", 3),
        (4, "F", 0),
    ):
        graph.add_atom(atom, element, isotope=isotope)
        graph.add_bond(0, atom)
    graph.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), 1))
    assert graph != graph.enantiomer()


@pytest.mark.parametrize("graph_type", [MolGraph, StereoMolGraph])
def test_rdkit_isotope_roundtrip(graph_type):
    rd_mol = Chem.MolFromSmiles("[2H][13CH3]")
    assert rd_mol is not None
    if graph_type is StereoMolGraph:
        graph = graph_type.from_rdmol(rd_mol, resonance=False)
    else:
        graph = graph_type.from_rdmol(rd_mol)
    assert sorted(graph.isotopes) == [2, 13]
    converted, _ = mol_graph_to_rdmol(graph)
    assert sorted(atom.GetIsotope() for atom in converted.GetAtoms()) == [2, 13]


@pytest.mark.parametrize("graph_type", [CondensedReactionGraph, StereoCondensedReactionGraph])
def test_reaction_graph_isotopes(graph_type):
    base_type = StereoMolGraph if graph_type is StereoCondensedReactionGraph else MolGraph
    reactant = base_type()
    reactant.add_atom(0, "C", isotope=13)
    reactant.add_atom(1, "H", isotope=2)
    product = reactant.copy()
    product.add_bond(0, 1)

    reaction = graph_type.from_graphs(reactant, product)
    assert reaction.isotopes == (13, 2)
    assert reaction.reactant(keep_attributes=False).isotopes == (13, 2)
    assert reaction.product(keep_attributes=False).isotopes == (13, 2)
