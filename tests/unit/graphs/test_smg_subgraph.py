"""Stereo subgraphs retain lone-pair slots but require every referenced atom."""

import pytest

from stereomolgraph import StereoCondensedReactionGraph, StereoMolGraph
from stereomolgraph.stereodescriptors import PlanarBond, Tetrahedral


@pytest.fixture(params=[StereoMolGraph, StereoCondensedReactionGraph])
def graph_type(request):
    return request.param


@pytest.fixture(
    params=[
        Tetrahedral((0, 1, 2, 3, None), 1),
        Tetrahedral((0, 1, 2, 3, None), -1),
        Tetrahedral((0, 1, 2, 3, None), None),
        Tetrahedral((0, 1, 2, 3, 4), 1),
        PlanarBond((0, None, 1, 2, 3, None), 0),
        PlanarBond((0, None, 1, 2, 3, None), None),
        PlanarBond((0, 4, 1, 2, 3, 5), 0),
    ],
    ids=[
        "atom-positive-lone-pair",
        "atom-negative-lone-pair",
        "atom-unspecified-lone-pair",
        "atom-no-placeholder",
        "bond-lone-pairs",
        "bond-unspecified-lone-pairs",
        "bond-no-placeholder",
    ],
)
def stereo_graph(graph_type, request):
    stereo = request.param
    graph = graph_type()
    for atom in range(7):
        graph.add_atom(atom, "C")
    if isinstance(stereo, Tetrahedral):
        for atom in stereo.atoms[1:]:
            if atom is not None:
                graph.add_bond(0, atom)
        graph.set_atom_stereo(stereo)
    else:
        graph.add_bond(1, 2)
        for center, ligands in ((1, stereo.atoms[:2]), (2, stereo.atoms[4:])):
            for atom in ligands:
                if atom is not None:
                    graph.add_bond(center, atom)
        graph.set_bond_stereo(stereo)
    return graph, stereo


@pytest.mark.parametrize(
    "selection_type", [list, set, iter], ids=["list", "set", "iterator"]
)
def test_subgraph_preserves_complete_stereo(stereo_graph, selection_type):
    graph, stereo = stereo_graph
    selected = {atom for atom in stereo.atoms if atom is not None}
    subgraph = graph.subgraph(selection_type(sorted(selected)))

    assert set(subgraph.atoms) == selected
    assert 6 not in subgraph.atoms  # Unrelated atom is excluded.
    assert len(subgraph.stereo) == 1
    actual = next(iter(subgraph.stereo.values()))
    # Stereo equality treats unspecified parity as a wildcard; inspect fields.
    assert type(actual) is type(stereo)
    assert actual.atoms == stereo.atoms
    assert actual.parity == stereo.parity
    assert len(graph.atoms) == 7
    assert next(iter(graph.stereo.values())).atoms == stereo.atoms


@pytest.mark.parametrize(
    "remove_site", [False, True], ids=["missing-ligand", "missing-site"]
)
def test_subgraph_drops_stereo_referencing_missing_real_atoms(
    stereo_graph, remove_site
):
    graph, stereo = stereo_graph
    selected = {atom for atom in stereo.atoms if atom is not None}
    if isinstance(stereo, Tetrahedral):
        selected.remove(0 if remove_site else 1)
    else:
        selected.remove(1 if remove_site else 0)
    subgraph = graph.subgraph(selected)

    assert set(subgraph.atoms) == selected
    assert not subgraph.stereo
    assert len(graph.stereo) == 1
