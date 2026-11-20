"""Allene atom-stereo descriptors and RDKit round trips."""

import pytest
from rdkit import Chem

from stereomolgraph import StereoMolGraph
from stereomolgraph.algorithms.circular import circular_stereo_generator
from stereomolgraph.graph2rdmol import stereo_mol_graph_to_rdmol
from stereomolgraph.stereodescriptors import Allene


def _allene_mol(smiles: str = "FC(Cl)=C=C(Br)I", permutation: int = 1):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    centers = [
        atom
        for atom in mol.GetAtoms()
        if len(atom.GetNeighbors()) == 2
        and all(b.GetBondType() == Chem.BondType.DOUBLE for b in atom.GetBonds())
    ]
    assert len(centers) == 1
    center = centers[0]
    center.SetChiralTag(Chem.ChiralType.CHI_ALLENE)
    center.SetUnsignedProp("_chiralPermutation", permutation)
    return mol, center.GetIdx()


def _graph(permutation: int) -> tuple[StereoMolGraph, int]:
    mol, center = _allene_mol(permutation=permutation)
    return StereoMolGraph.from_rdmol(mol, resonance=False), center


def test_allene_descriptor_permutations():
    stereo = Allene((0, 1, 2, 3, 4, 5, 6), 1)
    assert stereo == Allene((0, 2, 1, 5, 6, 3, 4), 1)
    assert stereo == Allene((0, 1, 2, 4, 3, 6, 5), 1)
    assert stereo == Allene((0, 1, 2, 4, 3, 5, 6), -1)
    assert stereo != stereo.invert()
    assert hash(stereo) == hash(Allene((0, 2, 1, 5, 6, 3, 4), 1))
    assert len(stereo.get_isomers()) == 2


@pytest.mark.parametrize("permutation,parity", [(1, 1), (2, -1), (0, None)])
def test_rdkit_allene_roundtrip(permutation, parity):
    mol, center = _allene_mol(permutation=permutation)
    graph = StereoMolGraph.from_rdmol(mol, resonance=False)
    stereo = graph.get_atom_stereo(center)
    assert isinstance(stereo, Allene)
    assert stereo.central_atom == center
    assert stereo.parity == parity
    assert len(stereo.atoms) == 7

    restored, rdidx_to_graphid = stereo_mol_graph_to_rdmol(graph)
    restored_center = next(
        rdidx for rdidx, graph_id in rdidx_to_graphid.items() if graph_id == center
    )
    rd_atom = restored.GetAtomWithIdx(restored_center)
    assert rd_atom.GetChiralTag() == Chem.ChiralType.CHI_ALLENE
    assert rd_atom.GetUnsignedProp("_chiralPermutation") == permutation
    assert all(b.GetBondType() == Chem.BondType.DOUBLE for b in rd_atom.GetBonds())

    roundtrip = StereoMolGraph.from_rdmol(restored, resonance=False)
    assert graph == roundtrip


def test_allene_enantiomers_appear_only_at_radius_two():
    first, center = _graph(1)
    second, _ = _graph(2)
    it1 = circular_stereo_generator(first)
    it2 = circular_stereo_generator(second)

    assert tuple(next(it1)) == tuple(next(it2))  # radius 0
    assert tuple(next(it1)) == tuple(next(it2))  # radius 1
    at_radius_two_first = next(it1)
    at_radius_two_second = next(it2)
    atom_position = list(first.atoms).index(center)
    assert at_radius_two_first[atom_position] != at_radius_two_second[atom_position]
    assert first != second
    first.freeze()
    second.freeze()
    assert hash(first) != hash(second)


def test_allene_relabeling_and_enantiomer():
    original, center = _graph(1)
    new_ids = {atom: atom + 100 for atom in original.atoms}
    relabeled = original.relabel_atoms(new_ids)
    assert relabeled == original
    assert isinstance(relabeled.get_atom_stereo(new_ids[center]), Allene)

    converted, _ = stereo_mol_graph_to_rdmol(relabeled)
    restored = StereoMolGraph.from_rdmol(converted, resonance=False)
    assert relabeled == restored

    opposite, _ = _graph(2)
    assert original.enantiomer() == opposite


def test_implicit_terminal_hydrogen():
    mol, center = _allene_mol("FC=C=C(Br)I", permutation=1)
    graph = StereoMolGraph.from_rdmol(mol, resonance=False)
    stereo = graph.get_atom_stereo(center)
    assert isinstance(stereo, Allene)
    assert stereo.atoms.count(None) == 1

    output, _ = stereo_mol_graph_to_rdmol(graph)
    Chem.SanitizeMol(output)
    result = StereoMolGraph.from_rdmol(output, resonance=False)
    assert result == graph


def test_allene_without_permutation_is_unspecified():
    mol, center = _allene_mol()
    atom = mol.GetAtomWithIdx(center)
    atom.ClearProp("_chiralPermutation")
    graph = StereoMolGraph.from_rdmol(mol, resonance=False)
    assert graph.get_atom_stereo(center).parity is None
