import numpy as np
import pytest

from stereomolgraph import (
    Bond,
    MolGraph,
)
from stereomolgraph.coords import Geometry
from stereomolgraph.graph2rdmol import mol_graph_to_rdmol
from stereomolgraph.periodic_table import PERIODIC_TABLE as PTOE

REACTION_SMILES_4007 = (
    "[CH2:1]1[CH:2]=[CH:3][N:4]([CH3:5])[CH:6]=[C:7]1[C:8]([NH2:9])=[O:10]."
    "[CH3:11][O:12][C:13](=[O:14])[C:15]([CH3:16])=[O+:17][O-:18]>>"
    "[CH2:1]1[C@H:2]2[C@@H:3]([N:4]([CH3:5])[CH:6]=[C:7]1[C:8]([NH2:9])=[O:10])"
    "[O:18][O:17][C:15]2([C:13]([O:12][CH3:11])=[O:14])[CH3:16]"
)


class TestMolGraph:
    _TestClass: type[MolGraph] = MolGraph

    @pytest.fixture
    def enantiomer_graph1(self, enantiomer_geos):
        return self._TestClass.from_geometry(enantiomer_geos[0])

    @pytest.fixture
    def enantiomer_graph2(self, enantiomer_geos):
        return self._TestClass.from_geometry(enantiomer_geos[1])

    @pytest.fixture
    def water_graph(self, water_geo):
        return self._TestClass.from_geometry(water_geo)

    @pytest.fixture
    def empty_mol_graph(self):
        return self._TestClass()

    @pytest.fixture
    def mol_graph(self):
        mol_graph = self._TestClass()
        mol_graph.add_atom(0, atom_type="C")
        mol_graph.add_atom(1, atom_type="H")
        mol_graph.add_atom(2, atom_type="O")
        mol_graph.add_bond(0, 1, bond_order=1)
        return mol_graph

    @pytest.fixture
    def chiral_product_geo1(self, data_path):
        filepath = data_path / "disrot_reaction" / "(Z)-(4S)-3,4-Dichlor-2-pentene.xyz"

        return Geometry.from_xyz_file(filepath)

    @pytest.fixture
    def chiral_product_graph1(self, chiral_product_geo1):
        return self._TestClass.from_geometry(chiral_product_geo1)

    @pytest.fixture
    def chiral_product_geo2(self, data_path):
        filepath = data_path / "conrot_reaction/(Z)-(4S)-3,4-Dichlor-2-pentene.xyz"
        return Geometry.from_xyz_file(filepath)

    @pytest.fixture
    def chiral_product_graph2(self, chiral_product_geo2):
        return self._TestClass.from_geometry(chiral_product_geo2)

    @pytest.fixture
    def chiral_reactant_geo(self, data_path):
        filepath = (
            data_path
            / "conrot_reaction/(2S,3S)-1,1-Dichlor-2,3-dimethylcyclopropane.xyz"
        )

        return Geometry.from_xyz_file(filepath)

    @pytest.fixture
    def chiral_reactant_graph(self, chiral_reactant_geo):
        return self._TestClass.from_geometry(chiral_reactant_geo)

    def test_len(self, enantiomer_graph1):
        assert enantiomer_graph1 is not None
        assert len(enantiomer_graph1) == 8

    def test_add_atom(self, empty_mol_graph, *args, **kwargs):
        empty_mol_graph.add_atom(0, atom_type="H", *args, **kwargs)
        empty_mol_graph.add_atom(1, atom_type=PTOE["C"], *args, **kwargs)
        assert empty_mol_graph

    def test_remove_atom(self, enantiomer_graph1):
        enantiomer_graph1.remove_atom(0)
        assert 0 not in enantiomer_graph1.atoms

    def test_add_bond(self, enantiomer_graph1):
        enantiomer_graph1.add_bond(0, 6)
        assert Bond({0, 6}) in enantiomer_graph1.bonds

    def test_remove_bond(self, water_graph):
        water_graph.remove_bond(0, 1)
        assert Bond({0, 1}) not in water_graph.bonds

    def test_atoms(self, enantiomer_graph1):
        assert set(enantiomer_graph1.atoms) == {0, 1, 2, 3, 4, 5, 6, 7}

    def test_bonds(self, mol_graph):
        assert Bond((0, 1)) in mol_graph.bonds

    def test_get_atom_attribute(self, mol_graph):
        assert mol_graph.get_atom_attribute(1, attr="atom_type") == PTOE["H"]

    def test_get_bond_attribute(self, mol_graph):
        assert mol_graph.get_bond_attribute(0, 1, attr="bond_order") == 1

    def test_get_atom_attributes(self, mol_graph):
        assert mol_graph.get_atom_attributes(1) == {"atom_type": PTOE["H"]}

    def test_get_atoms_with_attributes(self, mol_graph):
        assert mol_graph.atoms_with_attributes == {
            0: {"atom_type": PTOE["C"]},
            1: {"atom_type": PTOE["H"]},
            2: {"atom_type": PTOE["O"]},
        }

    def test_set_atom_attribute(self, mol_graph):
        mol_graph.set_atom_attribute(1, attr="test_attr", value="test")
        assert mol_graph.get_atom_attribute(1, attr="test_attr") == "test"
        mol_graph.set_atom_attribute(1, attr="atom_type", value="He")
        assert mol_graph.get_atom_attribute(1, attr="atom_type") == PTOE["He"]
        with pytest.raises(ValueError):
            mol_graph.set_atom_attribute(1, attr="atom_type", value="test")

    def test_set_bond_attribute(self, mol_graph):
        mol_graph.set_bond_attribute(0, 1, attr="lengh", value="very_long")
        assert mol_graph.get_bond_attribute(0, 1, attr="lengh") == "very_long"
        mol_graph.set_bond_attribute(0, 1, attr="bond_order", value=13)
        assert mol_graph.get_bond_attribute(0, 1, attr="bond_order") == 13

    def test_delete_atom_attribute(self, mol_graph):
        mol_graph.set_atom_attribute(1, attr="test_attr", value="test")
        assert mol_graph.get_atom_attribute(1, attr="test_attr") is not None
        mol_graph.delete_atom_attribute(1, attr="test_attr")
        assert mol_graph.get_atom_attribute(1, attr="test_attr") is None

    def test_delete_bond_attribute(self, mol_graph):
        mol_graph.delete_bond_attribute(0, 1, attr="bond_order")
        assert mol_graph.get_bond_attribute(0, 1, attr="bond_oder") is None

    def test_get_bond_attributes(self, mol_graph):
        assert mol_graph.get_bond_attributes(0, 1) == {"bond_order": 1}

    def test_get_bonds_with_attributes(self, mol_graph):
        assert mol_graph.bonds_with_attributes == {Bond({0, 1}): {"bond_order": 1}}

    def test_connectivity_matrix(self, mol_graph):
        assert np.array_equal(
            mol_graph.connectivity_matrix(),
            np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=int),
        )

    def test_relabel_atoms(self, enantiomer_graph1):
        mapping = {0: 15, 5: 10, 3: 99}
        enantiomer_graph1.relabel_atoms(mapping, copy=False)
        assert all(atom in enantiomer_graph1.atoms for atom in mapping.values())
        assert all(atom not in enantiomer_graph1.atoms for atom in mapping.keys())

    def test_relabel_atoms_copy(self, enantiomer_graph1):
        mapping = {0: 15, 5: 10, 3: 99}
        new_graph = enantiomer_graph1.relabel_atoms(mapping, copy=True)
        assert all(atom in new_graph.atoms for atom in mapping.values())
        assert all(atom not in new_graph.atoms for atom in mapping.keys())

    def test_connected_components(self, mol_graph):
        assert [i for i in mol_graph.connected_components()] == [{0, 1}, {2}]

    def test_subgraph(self, mol_graph):
        subgraph1 = mol_graph.subgraph([0, 1])
        assert (
            1 in subgraph1.atoms
            and 0 in subgraph1.atoms
            and Bond({0, 1}) in subgraph1.bonds
            and 2 not in subgraph1.atoms
        )
        subgraph2 = mol_graph.subgraph([1, 2])
        assert (
            1 in subgraph2.atoms
            and 2 in subgraph2.atoms
            and Bond({1, 2}) not in subgraph2.bonds
            and 0 not in subgraph2.atoms
        )

    def test_copy(self, enantiomer_graph1):
        copied_mol_graph = enantiomer_graph1.copy()
        assert id(copied_mol_graph) != id(enantiomer_graph1)

    def test_compose(self, water_graph, mol_graph, empty_mol_graph):
        comp_graph = self._TestClass.compose([water_graph, empty_mol_graph])
        assert comp_graph.atom_types == water_graph.atom_types
        assert comp_graph.atoms == water_graph.atoms
        assert comp_graph.bonds == water_graph.bonds

        comp_graph = self._TestClass.compose([water_graph, mol_graph])
        assert comp_graph.atom_types == mol_graph.atom_types
        assert Bond((0, 2)) in comp_graph.bonds
        assert Bond((0, 1)) in comp_graph.bonds
        assert comp_graph.get_bond_attribute(0, 1, "bond_order") == 1

    def test_from_composed_chiral_molgraphs(
        self, chiral_product_graph1, chiral_product_graph2
    ):
        relabel_mapping = {
            atom: atom + chiral_product_graph1.n_atoms
            for atom in chiral_product_graph2.atoms
        }
        chiral_product_graph2.relabel_atoms(relabel_mapping, copy=False)

        combined = self._TestClass.compose(
            [chiral_product_graph1, chiral_product_graph2]
        )

        assert (
            combined.atoms_with_attributes
            == chiral_product_graph1.atoms_with_attributes
            | chiral_product_graph2.atoms_with_attributes
        )
        assert (
            combined.bonds_with_attributes
            == chiral_product_graph1.bonds_with_attributes
            | chiral_product_graph2.bonds_with_attributes
        )

    def test_to_rdmol(self, water_graph):
        rdmol, _ = mol_graph_to_rdmol(water_graph)
        assert (
            tuple([Atom.GetAtomicNum() for Atom in rdmol.GetAtoms()])
            == water_graph.atom_types
        )
        assert {
            Bond((rd_b.GetBeginAtomIdx(), rd_b.GetEndAtomIdx()))
            for rd_b in rdmol.GetBonds()
        } == {Bond(b) for b in water_graph.bonds}

    @pytest.mark.parametrize(
        "inchi",
        [
            (r"InChI=1S/C3H8O/c1-3(2)4/h3-4H,1-2H3"),
            (r"InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3"),
        ],
        ids=["isopropanol", "caffeine"],
    )
    def test_from_rdmol_to_rdmol_not_chiral(self, inchi):
        # RDKit construction tests have been moved to `test_rdkit_conversion.py`.
        pytest.skip("RDKit construction tests moved to test_rdkit_conversion.py")

    def test_equality_relabeled_water(self, water_graph):
        assert water_graph == water_graph.copy()
        assert water_graph == water_graph.relabel_atoms({0: 1, 1: 0, 2: 13})

    def test_equality(
        self,
        chiral_product_graph1,
        chiral_product_graph2,
        chiral_reactant_graph,
    ):
        assert chiral_product_graph1 == chiral_product_graph2
        assert chiral_product_graph1 != chiral_reactant_graph != chiral_product_graph2

    def test_hash_enantiomers(self, enantiomer_graph1, enantiomer_graph2):
        assert hash(enantiomer_graph1.copy(frozen=True)) == hash(
            enantiomer_graph2.copy(frozen=True)
        )

    def test_hash_relabel(self, water_graph):
        relabel_water = water_graph.relabel_atoms({0: 1, 1: 0, 2: 13}, copy=True)
        assert hash(water_graph.copy(frozen=True)) == hash(
            relabel_water.copy(frozen=True)
        )
