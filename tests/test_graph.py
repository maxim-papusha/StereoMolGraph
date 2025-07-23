from collections import defaultdict
from copy import deepcopy
from itertools import permutations

import numpy as np
import pytest
import rdkit.Chem  # type: ignore

from stereomolgraph import PERIODIC_TABLE as PTOE
from stereomolgraph import (Bond,
                            MolGraph,
                            StereoMolGraph,
                            CondensedReactionGraph,
                            StereoCondensedReactionGraph)
from stereomolgraph.coords import Geometry, are_planar
from stereomolgraph.graphs.crg import Change
from stereomolgraph.graphs.scrg import ChangeDict

from stereomolgraph.stereodescriptors import (
    AtropBond,
    Octahedral,
    PlanarBond,
    Tetrahedral,
    TrigonalBipyramidal,
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
        filepath = data_path /"conrot_reaction/(2S,3S)-1,1-Dichlor-2,3-dimethylcyclopropane.xyz"
        
        return Geometry.from_xyz_file(filepath)

    @pytest.fixture
    def chiral_reactant_graph(self, chiral_reactant_geo):
        return self._TestClass.from_geometry(chiral_reactant_geo)

    def test_len(self, enantiomer_graph1):
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
        assert (
            mol_graph.get_atom_attribute(1, attr="test_attr")
            is not None
        )
        mol_graph.delete_atom_attribute(1, attr="test_attr")
        assert (
            mol_graph.get_atom_attribute(1, attr="test_attr")
            is None
        )

    def test_delete_bond_attribute(self, mol_graph):
        mol_graph.delete_bond_attribute(0, 1, attr="bond_order")
        assert (
            mol_graph.get_bond_attribute(0, 1, attr="bond_oder")
            is None
        )

    def test_get_bond_attributes(self, mol_graph):
        assert mol_graph.get_bond_attributes(0, 1) == {"bond_order": 1}

    def test_get_bonds_with_attributes(self, mol_graph):
        assert mol_graph.bonds_with_attributes == {
            Bond({0, 1}): {"bond_order": 1}
        }

    def test_connectivity_matrix(self, mol_graph):
        assert np.array_equal(
            mol_graph.connectivity_matrix(),
            np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=int),
        )

    def test_relabel_atoms(self, enantiomer_graph1):
        mapping = {0: 15, 5: 10, 3: 99}
        enantiomer_graph1.relabel_atoms(mapping, copy=False)
        assert all(
            atom in enantiomer_graph1.atoms for atom in mapping.values()
        )
        assert all(
            atom not in enantiomer_graph1.atoms for atom in mapping.keys()
        )

    def test_relabel_atoms_copy(self, enantiomer_graph1):
        mapping = {0: 15, 5: 10, 3: 99}
        new_graph = enantiomer_graph1.relabel_atoms(mapping, copy=True)
        assert all(atom in new_graph.atoms for atom in mapping.values())
        assert all(atom not in new_graph.atoms for atom in mapping.keys())

    #@pytest.mark.skip("Not Implemented")
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

    def test_get_isomorphic_mappings(self, water_graph, mol_graph):
        assert [] == [
            i for i in mol_graph.get_isomorphic_mappings(water_graph)
        ]
        
        assert all(
            mapping in ({0: 0, 1: 1, 2: 2}, {0: 0, 2: 1, 1: 2})
            for mapping in (
                i for i in water_graph.get_isomorphic_mappings(water_graph)
            )
        )

    def test_get_isomorphic_mappings_of_enantiomers(
        self, enantiomer_graph1, enantiomer_graph2
    ):
        assert 2 == len(
            list(enantiomer_graph1.get_isomorphic_mappings(enantiomer_graph2))
        )


    def test_get_automorphic_mappings(self, water_graph):
        assert all(
            mapping in ({0: 0, 1: 1, 2: 2}, {0: 0, 2: 1, 1: 2})
            for mapping in (i for i in water_graph.get_isomorphic_mappings(water_graph))
        )

    def test_compose(
        self, water_graph, mol_graph, empty_mol_graph
    ):

        comp_graph = self._TestClass.compose(
            [water_graph, empty_mol_graph]
        )
        assert comp_graph.atom_types == water_graph.atom_types
        assert comp_graph.atoms == water_graph.atoms
        assert comp_graph.bonds == water_graph.bonds

        comp_graph = self._TestClass.compose(
            [water_graph, mol_graph]
        )
        assert comp_graph.atom_types == mol_graph.atom_types
        assert Bond((0, 2)) in comp_graph.bonds
        assert Bond((0, 1)) in comp_graph.bonds
        assert (
            comp_graph.get_bond_attribute(0, 1, "bond_order") == 1
        )

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
        rdmol, _ = water_graph._to_rdmol()
        assert [Atom.GetAtomicNum() for Atom in rdmol.GetAtoms()] == [
            atom_type.atomic_nr for atom_type in water_graph.atom_types
        ]
        assert ({Bond(
                (rd_b.GetBeginAtomIdx(), rd_b.GetEndAtomIdx()))
                for rd_b in rdmol.GetBonds()}
                == {Bond(b) for b in water_graph.bonds})
        
    def test_from_rdmol(self):

        rdmol = rdkit.Chem.MolFromSmiles("[C:1]([O:2][C:33]([C:4]([O:5][H:13])"
                                         "([H:111])[H:12])([H:9])[H:10])([H:6])"
                                         "([H:77])[H:8]", sanitize=False)
        mol_graph = self._TestClass.from_rdmol(rdmol, use_atom_map_number=True)

        assert set(mol_graph.atoms) == set((1,2,33,4,5,6,77,8,9,10,111,12,13))
        assert mol_graph.get_atom_attribute(1, "atom_type") == PTOE["C"]
        assert mol_graph.get_atom_attribute(2, "atom_type") == PTOE["O"]
        assert mol_graph.get_atom_attribute(33, "atom_type") == PTOE["C"]
        assert mol_graph.get_atom_attribute(4, "atom_type") == PTOE["C"]
        assert mol_graph.get_atom_attribute(5, "atom_type") == PTOE["O"]
        assert mol_graph.get_atom_attribute(6, "atom_type") == PTOE["H"]
        assert mol_graph.get_atom_attribute(77, "atom_type") == PTOE["H"]
        assert mol_graph.get_atom_attribute(8, "atom_type") == PTOE["H"]
        assert mol_graph.get_atom_attribute(9, "atom_type") == PTOE["H"]
        assert mol_graph.get_atom_attribute(10, "atom_type") == PTOE["H"]
        assert mol_graph.get_atom_attribute(111, "atom_type") == PTOE["H"]
        assert mol_graph.get_atom_attribute(12, "atom_type") == PTOE["H"]
        assert mol_graph.get_atom_attribute(13, "atom_type") == PTOE["H"]

        assert mol_graph.has_bond(1, 2)
        assert mol_graph.has_bond(2, 33)
        assert mol_graph.has_bond(33, 4)
        assert mol_graph.has_bond(4, 5)
        assert mol_graph.has_bond(1, 6)
        assert mol_graph.has_bond(1, 77)
        assert mol_graph.has_bond(1, 8)
        assert mol_graph.has_bond(33, 9)
        assert mol_graph.has_bond(33, 10)
        assert mol_graph.has_bond(4, 111)
        assert mol_graph.has_bond(4, 12)
        assert mol_graph.has_bond(5, 13)

    @pytest.mark.parametrize("inchi", [
    (r"InChI=1S/C3H8O/c1-3(2)4/h3-4H,1-2H3"),
    (r"InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3")
    ], ids = ["isopropanol",
              "caffeine"],)
    def test_from_rdmol_to_rdmol_not_chiral(self, inchi):
        rdmol = rdkit.Chem.MolFromInchi(inchi, sanitize=False, removeHs=False)
        assert inchi == rdkit.Chem.MolToInchi(rdmol, treatWarningAsError=True) # type: ignore
        
        molgraph = self._TestClass.from_rdmol(rdmol)
        rdmol2, _ = molgraph._to_rdmol(generate_bond_orders=True) 
        # TODO: check why bond orders are needed

        assert inchi == rdkit.Chem.MolToInchi(rdmol2, treatWarningAsError=True) # type: ignore

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
        assert (
            chiral_product_graph1
            != chiral_reactant_graph
            != chiral_product_graph2
        )

    def test_hash_enantiomers(self, enantiomer_graph1, enantiomer_graph2):
        assert (hash(enantiomer_graph1)
                == hash(enantiomer_graph2))

    def test_hash_relabel(self, water_graph):
        relabel_water = water_graph.relabel_atoms(
            {0: 1, 1: 0, 2: 13}, copy=True
        )
        assert (hash(water_graph)
                == hash(relabel_water))


class TestCondensedReactionGraph(TestMolGraph):
    _TestClass: type[CondensedReactionGraph] = CondensedReactionGraph

    @pytest.fixture
    def crg(self):
        crg = self._TestClass()
        crg.add_atom(0, atom_type="C", atom_label="atom_label_C1")
        crg.add_atom(1, atom_type="H")
        crg.add_atom(2, atom_type="O")
        crg.add_atom(3, atom_type="H")
        crg.add_bond(0, 1, reaction=Change.FORMED)
        crg.add_bond(0, 2, bond_label="bond_label_C1O3")
        crg.add_bond(2, 3, reaction=Change.BROKEN)
        return crg

    def test_bonds(self, chiral_reactant_graph, chiral_product_graph1,
                            chiral_reaction_scrg1):
        unchanged_b = {bond for bond
                           in chiral_reaction_scrg1.bonds
                           if bond
                           not in chiral_reaction_scrg1.get_broken_bonds()
                           and bond
                           not in chiral_reaction_scrg1.get_formed_bonds()}

        formed_b = chiral_reaction_scrg1.get_formed_bonds()
        broken_b = chiral_reaction_scrg1.get_broken_bonds()

        assert set(chiral_reactant_graph.bonds) == unchanged_b | broken_b
        assert set(chiral_product_graph1.bonds) == unchanged_b | formed_b

    def test_add_bond_error(self, crg):
        with pytest.raises(TypeError):
            crg.add_bond(0, 1, reaction="test")

    def test_add_bond_with_reaction_attr(self, crg):
        crg.add_bond(0, 1, reaction=Change.FORMED, bond_order=1)
        assert (
            crg.get_bond_attribute(0, 1, attr="reaction") == Change.FORMED
        )
        assert crg.get_bond_attribute(0, 1, attr="bond_order") == 1

    def test_set_bond_attribute_reaction_exception(self, crg):
        with pytest.raises(ValueError):
            crg.set_bond_attribute(0, 1, attr="reaction", value="test")

    def test_set_bond_attribute_reaction(self, crg):
        crg.set_bond_attribute(0, 1, attr="reaction", value=Change.FORMED)
        assert (
            crg.get_bond_attribute(0, 1, attr="reaction") == Change.FORMED
        )

    def test_set_bond_attribute(self, crg):
        crg.set_bond_attribute(0, 1, attr="bond_label", value="test")
        assert crg.get_bond_attribute(0, 1, attr="bond_label") == "test"

    def test_add_formed_bond(self, crg):
        crg.add_formed_bond(1, 2)
        assert (
            crg.get_bond_attribute(1, 2, attr="reaction") == Change.FORMED
        )

    def test_add_broken_bond(self, crg):
        crg.add_broken_bond(1, 2)
        assert (
            crg.get_bond_attribute(2, 1, attr="reaction") == Change.BROKEN
        )

    def test_get_formed_bonds(self, crg):
        assert set(crg.get_formed_bonds()) == {Bond((0, 1)),}

    def test_get_broken_bonds(self, crg):
        assert set(crg.get_broken_bonds()) == {Bond((2, 3)),}

    def test_active_atoms_crg(self, crg):
        assert set(crg.active_atoms(additional_layer = 0)) == {0, 1, 2, 3}
        assert set(crg.active_atoms(additional_layer = 1)) == {0, 1, 2, 3}

        crg.add_atom(4, "O")
        crg.add_bond(4, 2)
        crg.add_atom(5, "H")
        crg.add_bond(5, 4)

        assert set(crg.active_atoms(additional_layer = 0)) == {0, 1, 2, 3}
        assert set(crg.active_atoms(additional_layer = 1)) == {0, 1, 2, 3, 4}

    def test_active_atoms(self):
        crg = self._TestClass()
        crg.add_atom(1, "C")
        crg.add_atom(0, "C")
        crg.add_broken_bond(0,1)
        for i in range(2, 10):
            crg.add_atom(i, "C")
            crg.add_bond(i, i-1)
        for j in range(0, 5):
            assert set(crg.active_atoms( additional_layer=j)) == {*range(j+2)}

    def test_reactant_with_attributes(self, crg):
        crg_copy = crg.copy()
        for bond in crg.get_formed_bonds():
            crg_copy.remove_bond(*bond)
        for bond in crg.get_broken_bonds():
            crg_copy.delete_bond_attribute(*bond, attr="reaction")

        expected_result = self.__class__.__bases__[0]._TestClass(crg_copy)
        assert (
            expected_result.bonds_with_attributes
            == crg.reactant(keep_attributes=True).bonds_with_attributes
        )

    def test_reactant_without_attributes(self, crg):
        crg_copy = crg.copy()
        for bond in crg.get_formed_bonds():
            crg_copy.remove_bond(*bond)
        for bond in crg.get_broken_bonds():
            crg_copy.delete_bond_attribute(*bond, attr="reaction")
        for bond, attr_dict in crg_copy.bonds_with_attributes.items():
            to_delete = [bond_attr for bond_attr in attr_dict.keys()]
            for bond_attr in to_delete:
                crg_copy.delete_bond_attribute(*bond, attr=bond_attr)
        for atom, attr_dict in crg_copy.atoms_with_attributes.items():
            to_delete = [
                atom_attr
                for atom_attr in attr_dict.keys()
                if atom_attr != "atom_type"
            ]
            for atom_attr in to_delete:
                crg_copy.delete_atom_attribute(atom, atom_attr)
        expected_result = self.__class__.__bases__[0]._TestClass(crg_copy)
        assert (
            expected_result.bonds_with_attributes
            == crg.reactant(keep_attributes=False).bonds_with_attributes
        )

    def test_product_with_attributes(self, crg):
        crg_copy = crg.copy()
        for bond in crg.get_broken_bonds():
            crg_copy.remove_bond(*bond)
        for bond in crg.get_formed_bonds():
            crg_copy.delete_bond_attribute(*bond, attr="reaction")

        expected_result = self.__class__.__bases__[0]._TestClass(crg_copy)
        assert (
            expected_result.bonds_with_attributes
            == crg.product(keep_attributes=True).bonds_with_attributes
        )

    def test_product_without_attributes(self, crg):
        crg_copy = crg.copy()
        for bond in crg.get_broken_bonds():
            crg_copy.remove_bond(*bond)
        for bond in crg.get_formed_bonds():
            crg_copy.delete_bond_attribute(*bond, attr="reaction")
        for bond, attr_dict in crg_copy.bonds_with_attributes.items():
            to_delete = [bond_attr for bond_attr in attr_dict.keys()]
            for bond_attr in to_delete:
                crg_copy.delete_bond_attribute(*bond, attr=bond_attr)
        for atom, attr_dict in crg_copy.atoms_with_attributes.items():
            to_delete = [
                atom_attr
                for atom_attr in attr_dict.keys()
                if atom_attr != "atom_type"
            ]
            for atom_attr in to_delete:
                crg_copy.delete_atom_attribute(atom, atom_attr)
        expected_result = self.__class__.__bases__[0]._TestClass(crg_copy)
        assert (
            expected_result.bonds_with_attributes
            == crg.product(keep_attributes=False).bonds_with_attributes
        )

    @pytest.fixture
    def chiral_reaction_scrg1(self, chiral_reactant_geo, chiral_product_geo1):
        return self._TestClass.from_reactant_and_product_geometry(
            chiral_reactant_geo, chiral_product_geo1
        )

    @pytest.fixture
    def chiral_reaction_scrg2(self, chiral_reactant_geo, chiral_product_geo2):
        return self._TestClass.from_reactant_and_product_geometry(
            chiral_reactant_geo, chiral_product_geo2
        )


    def test_reverse_reaction(self, chiral_reaction_scrg1):
        reversed_reaction = chiral_reaction_scrg1.reverse_reaction()
        assert (
            reversed_reaction.get_broken_bonds()
            == chiral_reaction_scrg1.get_formed_bonds()
        )
        assert (
            reversed_reaction.get_formed_bonds()
            == chiral_reaction_scrg1.get_broken_bonds()
        )
        
        double_reverset_reaction = reversed_reaction.reverse_reaction()
        assert double_reverset_reaction == chiral_reaction_scrg1

    def test_isomorphism_same_reactant_and_product_without_ts(
        self, chiral_reaction_scrg1, chiral_reaction_scrg2
    ):
        assert chiral_reaction_scrg1.product().is_isomorphic(
            chiral_reaction_scrg2.product()
        )
        assert chiral_reaction_scrg1.reactant().is_isomorphic(
            chiral_reaction_scrg2.reactant()
        )
        assert chiral_reaction_scrg1.is_isomorphic(chiral_reaction_scrg2)


class TestPlanar:
    def test_are_planar_true(self):
        coords = np.array([[0.0, 0.0, 0.0],[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[1.0, 1.0, 0.1]], dtype=np.float64)
        assert are_planar(coords, threshold=0.5)

    def test_are_planar_false(self):
        coords = np.array([[0.0, 0.0, 0.0],[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[1.0, 1.0, 1.0]], dtype=np.float64)
        assert not are_planar(coords, threshold=0.5)


class TestTetrahedral:
    def test_from_coords(self, enantiomer_geos):
        coords1 = enantiomer_geos[0].coords
        coords2 = enantiomer_geos[1].coords
        atoms = (3, 0, 1, 2, 4)

        stereo1 = Tetrahedral.from_coords(atoms, coords1.take(atoms, axis=0))
        stereo2 = Tetrahedral.from_coords(atoms, coords2.take(atoms, axis=0))
        assert stereo1.parity == 1
        assert stereo2.parity == -1

    def test_from_permuted_coords(self, enantiomer_geos):
        coords = enantiomer_geos[0].coords
        different_perms = {Tetrahedral.from_coords(atoms:=(3, *perm),
                                                   coords.take(atoms, axis=0))
                    for perm in permutations((0, 1, 2, 4))}
        assert len(different_perms) == 1

    def test_equality(self):
        stereo1 = Tetrahedral((6 ,0, 1, 2, 3), 1)
        stereo2 = Tetrahedral((6, 1, 2, 0, 3), 1)
        stereo3 = Tetrahedral((6, 0, 2, 1, 3), -1)
        stereo4 = Tetrahedral((6, 0, 2, 1, 3), 1)
        assert stereo1 == stereo2 == stereo3 != stereo4
        assert hash(stereo1) == hash(stereo2) == hash(stereo3) != hash(stereo4)

    def test_equality_with_none(self):
        stereo1 = Tetrahedral((6, 0, 1, 2, 3), None)
        stereo2 = Tetrahedral((6, 1, 2, 3, 0), None)
        assert stereo1 == stereo2

    def test_permutations(self):
        stereo1 = Tetrahedral((6, 0, 1, 2, 3), 1)
        stereo2 = Tetrahedral((6, 1, 2, 0, 3), 1)
        assert set(stereo1._perm_atoms()) == set(stereo2._perm_atoms())


class TestTrigonalBipyramidal:
    def test_from_coords(self, data_path):
        pcl5 = Geometry.from_xyz_file(data_path / "PCl5.xyz")
        result = TrigonalBipyramidal.from_coords(
            (0, 1, 2, 3, 4, 5), pcl5.coords.take((0, 1, 2, 3, 4, 5), axis=0)
        )
        assert result.parity == 1
        assert set(result.atoms) == {3, 4, 0, 1, 2, 5}

    def test_equality(self):
        stereo1 = TrigonalBipyramidal(
            (6, 0, 1, 2, 3, 4), 1
        )
        stereo2 = TrigonalBipyramidal(
            (6, 1, 0, 2, 3, 4), -1
        )
        stereo3 = TrigonalBipyramidal(
            (6, 1, 0, 3, 4, 2), -1
        )
        stereo4 = TrigonalBipyramidal(
            (6, 1, 0, 3, 4, 2), 1
        )
        stereo5 = deepcopy(stereo4)
        assert stereo1 == stereo2 == stereo3 != stereo4
        assert (
            hash(stereo1)
            == hash(stereo2)
            == hash(stereo3)
            != hash(stereo4)
            == hash(stereo5)
        )

    def test_equality_with_none(self):
        stereo1 = TrigonalBipyramidal(
            (0, 1, 2, 3, 4, 5), None
        )
        stereo2 = TrigonalBipyramidal(
            (1, 0, 2, 3, 4, 5), None
        )
        assert stereo1 == stereo2


class TestPlanarBond:
    def test_equality(self):
        stereo1 = PlanarBond((5, 4, 3, 2, 1, 0), 0)
        stereo2 = PlanarBond((4, 5, 3, 2, 0, 1), 0)
        stereo3 = PlanarBond((1, 0, 2, 3, 5, 4), 0)
        stereo4 = PlanarBond((4, 5, 3, 2, 1, 0), 0)
        stereo5 = deepcopy(stereo4)
        assert stereo1 == stereo2 == stereo3 != stereo4 == stereo5
        assert (
            hash(stereo1)
            == hash(stereo2)
            == hash(stereo3)
            != hash(stereo4)
            == hash(stereo5)
        )

    def test_equality_with_none(self):
        stereo3 = PlanarBond((1, 0, 2, 3, 5, 4), None)
        stereo4 = PlanarBond((5, 4, 3, 2, 1, 0), None)
        assert stereo3 == stereo4


class TestStereoMolGraph(TestMolGraph):
    _TestClass: type[StereoMolGraph] = StereoMolGraph

    def test_from_geometries1(self, chiral_reactant_graph):
        graph = chiral_reactant_graph
        expected_atom_stereo = {
            1: Tetrahedral((1, 0, 2, 3, 9), 1),
            0: Tetrahedral((0, 1, 2, 13, 14), -1),
            9: Tetrahedral((9, 1, 10, 11, 12), -1),
            5: Tetrahedral((5, 2, 6, 7, 8), -1),
        }
        expected_atom_types = [
            PTOE[atom]
            for atom in (
                "C",
                "C",
                "C",
                "H",
                "H",
                "C",
                "H",
                "H",
                "H",
                "C",
                "H",
                "H",
                "H",
                "Cl",
                "Cl",
            )
        ]
        expected_bonds = {Bond(pair) for pair in
            [(0, 1),
            (0, 2),
            (0, 13),
            (0, 14),
            (1, 2),
            (1, 3),
            (1, 9),
            (2, 4),
            (2, 5),
            (5, 6),
            (5, 7),
            (5, 8),
            (9, 10),
            (9, 11),
            (9, 12),]
        }
        assert graph.atom_types == tuple(expected_atom_types)
        assert set(graph.bonds) == expected_bonds
        assert all(
            graph.get_atom_stereo(key) == value
            for key, value in expected_atom_stereo.items()
        )

    def test_atom_stereo(self, chiral_product_graph1):
        expected = {
            1: Tetrahedral((1, 0, 3, 9, 13), 1),
            5: Tetrahedral((5, 2, 6, 7, 8), -1),
            9: Tetrahedral((9, 1, 10, 11, 12), -1),
        }
        expected2 = {
            Bond((0, 2)): PlanarBond(
                (1, 14, 0, 2, 4, 5), 0
            )
        }
        assert all(
            key in chiral_product_graph1._atom_stereo
            for key in set(expected.keys())
        )
        assert all(
            key in expected
            for key in set(chiral_product_graph1._atom_stereo.keys())
        )
        assert all(
            expected[key] == value
            for key, value in chiral_product_graph1._atom_stereo.items()
        )
        assert expected2 == chiral_product_graph1._bond_stereo

    def test_to_rdmol_double_bond(self):
        g1 = self._TestClass()
        g1.add_atom(0, atom_type="F")
        g1.add_atom(1, atom_type="H")
        g1.add_atom(2, atom_type="C")
        g1.add_atom(3, atom_type="C")
        g1.add_atom(4, atom_type="F")
        g1.add_atom(5, atom_type="H")

        g1.add_bond(0, 2)
        g1.add_bond(1, 2)
        g1.add_bond(2, 3)
        g1.add_bond(3, 4)
        g1.add_bond(3, 5)

        g2 = g1.copy()
        g3 = g1.copy()
        g1.set_bond_stereo(PlanarBond((0, 1, 2, 3, 4, 5),
                                                 0))
        g2.set_bond_stereo(PlanarBond((1, 0, 2, 3, 4, 5),
                                                    0))
        g3.set_bond_stereo(PlanarBond((0, 1, 2, 3, 4, 5),
                                                    None))
        rdmol_g1, idx_atom_map_dict_g1 = g1._to_rdmol()
        rdmol_g2, idx_atom_map_dict_g2 = g2._to_rdmol()
        rdmol_g3, idx_atom_map_dict_g3 = g3._to_rdmol()

        db1 = rdmol_g1.GetBondBetweenAtoms(2,3)
        stereo_atoms1 = {idx_atom_map_dict_g1[i] for i in db1.GetStereoAtoms()}
        assert stereo_atoms1 == {0, 4} or stereo_atoms1 == {1, 5}
        assert db1.GetStereo() == rdkit.Chem.rdchem.BondStereo.STEREOZ # type: ignore

        db2 = rdmol_g2.GetBondBetweenAtoms(2,3)
        stereo_atoms2 = {idx_atom_map_dict_g2[i] for i in db2.GetStereoAtoms()}
        assert stereo_atoms2 == {1, 4} or stereo_atoms2 == {0, 5}
        assert db2.GetStereo() == rdkit.Chem.rdchem.BondStereo.STEREOZ # type: ignore

        db3 = rdmol_g3.GetBondBetweenAtoms(2,3)

        assert db3.GetStereo() == rdkit.Chem.rdchem.BondStereo.STEREONONE # type: ignore

    def test_to_rdmol_tetrahedral(self):
        g = self._TestClass()
        g.add_atom(0, atom_type="C")
        g.add_atom(1, atom_type="H")
        g.add_atom(2, atom_type="F")
        g.add_atom(3, atom_type="Cl")
        g.add_atom(4, atom_type="Br")
        g.add_bond(0, 1)
        g.add_bond(0, 2)
        g.add_bond(0, 3)
        g.add_bond(0, 4)
        g.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), 1))

        mol, _ = g._to_rdmol()
        chiral_tag = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW # type: ignore
        assert mol.GetAtomWithIdx(0).GetChiralTag() == chiral_tag

        g.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), -1))
        mol, _ = g._to_rdmol()
        chiral_tag = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW # type: ignore
        assert mol.GetAtomWithIdx(0).GetChiralTag() == chiral_tag

    @pytest.mark.parametrize("inchi", [
        (r"InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m0/s1"),
        (r"InChI=1S/C2H2Cl2/c3-1-2-4/h1-2H/b2-1+"),
        (r"InChI=1S/C2H2Cl2/c3-1-2-4/h1-2H/b2-1-"),
        (r"InChI=1S/C4H6/c1-3-4-2/h3-4H,1-2H2")
    ], ids = ["(R)-Bromochlorofluoromethane",
              "Trans-1,2-Dichloroethylene",
              "Cis-1,2-Dichloroethylene",
              "Butadiene"],)
    def test_from_rdmol_to_rdmol_stereo(self, inchi):
        rdmol = rdkit.Chem.MolFromInchi(inchi, sanitize=False)
        rdmol = rdkit.Chem.AddHs(rdmol, explicitOnly=True)
        molgraph = self._TestClass.from_rdmol(rdmol)
        rdmol2, _ = molgraph._to_rdmol(generate_bond_orders=True)
        assert inchi == rdkit.Chem.MolToInchi(rdmol2, treatWarningAsError=True) # type: ignore

    @pytest.mark.parametrize("smiles", [
        "[H][Pt@SP1](F)(Cl)Br",
        "Cl[Pt@SP1](Cl)([NH3])[NH3]",
        "Cl[Pt@SP2](Cl)([NH3])[NH3]"
    ], ids = ["Hydridofluorochlorobromoplatinum(II)",
              "(SP-4-2)-diamminedichloroplatinum",
              "(SP-4-1)-diamminedichloroplatinum"],)
    def test_from_rdmol_to_rdmol_square_planar(self, smiles):
        rdmol = rdkit.Chem.MolFromSmiles(smiles, sanitize=True)
        rdmol = rdkit.Chem.AddHs(rdmol, explicitOnly=True)
        molgraph = self._TestClass.from_rdmol(rdmol)
        #raise Exception (molgraph.stereo)
        rdmol2, _ = molgraph._to_rdmol(generate_bond_orders=True,
                                       allow_charged_fragments=True)
        rdkit.Chem.SanitizeMol(rdmol2, sanitizeOps=rdkit.Chem.SanitizeFlags.SANITIZE_ALL)
        for atom in rdmol2.GetAtoms():
            atom.SetAtomMapNum(0)
        assert rdkit.Chem.MolToSmiles(rdmol) == rdkit.Chem.MolToSmiles(rdmol2)

    def test_from_rdmol_square_planar_different(self):
        smiles = ["Cl[Pt@SP1](Cl)([NH3])[NH3]", "Cl[Pt@SP2](Cl)([NH3])[NH3]"]
        rdmols = [rdkit.Chem.MolFromSmiles(s, sanitize=False) for s in smiles]
        rdmol = [rdkit.Chem.AddHs(mol, explicitOnly=True) for mol in rdmols]
        molgraphs = [self._TestClass.from_rdmol(mol) for mol in rdmol]
        assert molgraphs[0] != molgraphs[1]

    def test_from_rdmol_square_planar(self):
        smiles = ("C[Pt@SP1](F)(Cl)[H]", "C[Pt@SP2](Cl)(F)[H]",
                  "C[Pt@SP3](F)([H])Cl")
        mols = [rdkit.Chem.MolFromSmiles(i, sanitize=False) for i in smiles]
        mols = [rdkit.Chem.AddHs(i, explicitOnly=True) for i in mols]
        molgraphs = [self._TestClass.from_rdmol(i) for i in mols]
        assert all(molgraph == molgraphs[0] for molgraph in molgraphs)

    def test_from_rdmol_trigonal_bipyramidal(self):
        smiles = ("S[As@TB1](F)(Cl)(Br)N", "S[As@TB2](F)(Br)(Cl)N",
                  "S[As@TB3](F)(Cl)(N)Br", "S[As@TB4](F)(Br)(N)Cl",
                  "S[As@TB5](F)(N)(Cl)Br", "S[As@TB6](F)(N)(Br)Cl",
                  "S[As@TB7](N)(F)(Cl)Br", "S[As@TB8](N)(F)(Br)Cl",
                  "F[As@TB9](S)(Cl)(Br)N", "F[As@TB11](S)(Br)(Cl)N",
                  "F[As@TB10](S)(Cl)(N)Br", "F[As@TB12](S)(Br)(N)Cl",
                  "F[As@TB13](S)(N)(Cl)Br", "F[As@TB14](S)(N)(Br)Cl",
                  "F[As@TB15](Cl)(S)(Br)N", "F[As@TB20](Br)(S)(Cl)N",
                  "F[As@TB16](Cl)(S)(N)Br", "F[As@TB19](Br)(S)(N)Cl",
                  "F[As@TB17](Cl)(Br)(S)N", "F[As@TB18](Br)(Cl)(S)N")
        mols = [rdkit.Chem.MolFromSmiles(i, sanitize=False) for i in smiles]
        mols = [rdkit.Chem.AddHs(i, explicitOnly=True) for i in mols]
        molgraphs = [self._TestClass.from_rdmol(i) for i in mols]
        assert all(any( isinstance(stereo, TrigonalBipyramidal) for
                       stereo in molgraph.stereo.values())
                   for molgraph in molgraphs)
        assert all(molgraph == molgraphs[0] for molgraph in molgraphs)


    def test_from_rdmol_octahedral(self):
        smiles = ("O[Co@OH1](Cl)(C)(N)(F)P", "O[Co@OH2](Cl)(F)(N)(C)P",
                  "O[Co@OH3](Cl)(C)(N)(P)F", "O[Co@OH16](Cl)(F)(N)(P)C",
                  "O[Co@OH6](Cl)(C)(P)(N)F", "O[Co@OH18](Cl)(F)(P)(N)C",
                  "O[Co@OH19](Cl)(P)(C)(N)F", "O[Co@OH24](Cl)(P)(F)(N)C",
                  "O[Co@OH25](P)(Cl)(C)(N)F", "O[Co@OH30](P)(Cl)(F)(N)C",
                  "O[Co@OH4](Cl)(C)(F)(N)P", "O[Co@OH14](Cl)(F)(C)(N)P",
                  "O[Co@OH5](Cl)(C)(F)(P)N", "O[Co@OH15](Cl)(F)(C)(P)N",
                  "O[Co@OH7](Cl)(C)(P)(F)N", "O[Co@OH17](Cl)(F)(P)(C)N",
                  "O[Co@OH20](Cl)(P)(C)(F)N", "O[Co@OH23](Cl)(P)(F)(C)N",
                  "O[Co@OH26](P)(Cl)(C)(F)N", "O[Co@OH29](P)(Cl)(F)(C)N",
                  "O[Co@OH10](Cl)(N)(F)(C)P", "O[Co@OH8](Cl)(N)(C)(F)P",
                  "O[Co@OH11](Cl)(N)(F)(P)C", "O[Co@OH9](Cl)(N)(C)(P)F",
                  "O[Co@OH13](Cl)(N)(P)(F)C", "O[Co@OH12](Cl)(N)(P)(C)F",
                  "O[Co@OH22](Cl)(P)(N)(F)C", "O[Co@OH21](Cl)(P)(N)(C)F",
                  "O[Co@OH28](P)(Cl)(N)(F)C", "O[Co@OH27](P)(Cl)(N)(C)F")

        mols = [rdkit.Chem.MolFromSmiles(i, sanitize=False) for i in smiles]
        mols = [rdkit.Chem.AddHs(i, explicitOnly=True) for i in mols]
        molgraphs = [self._TestClass.from_rdmol(i) for i in mols]
        assert all(any( isinstance(stereo, Octahedral) for
                       stereo in molgraph.stereo.values())
                   for molgraph in molgraphs)
        assert all(molgraph == molgraphs[0] for molgraph in molgraphs)

    def test_from_rdmol_octahedral_compare(self):
        identical1 = ("Cl[Co@OH1](N)(N)(O)(Cl)Cl", "Cl[Co@OH2](N)(N)(O)(Cl)Cl",
                      "Cl[Co@OH3](N)(N)(O)(Cl)Cl", "Cl[Co@OH4](N)(N)(O)(Cl)Cl",
                      "Cl[Co@OH5](N)(N)(O)(Cl)Cl","Cl[Co@OH14](N)(N)(O)(Cl)Cl",
                      "Cl[Co@OH15](N)(N)(O)(Cl)Cl","Cl[Co@OH16](N)(N)(O)(Cl)Cl",
                      "Cl[Co@OH21](N)(N)(O)(Cl)Cl","Cl[Co@OH22](N)(N)(O)(Cl)Cl",
                      "Cl[Co@OH27](N)(N)(O)(Cl)Cl","Cl[Co@OH28](N)(N)(O)(Cl)Cl")
        identical2 = ("Cl[Co@OH6](N)(N)(O)(Cl)Cl","Cl[Co@OH7](N)(N)(O)(Cl)Cl",
                      "Cl[Co@OH17](N)(N)(O)(Cl)Cl","Cl[Co@OH18](N)(N)(O)(Cl)Cl",
                      "Cl[Co@OH19](N)(N)(O)(Cl)Cl","Cl[Co@OH20](N)(N)(O)(Cl)Cl",
                      "Cl[Co@OH23](N)(N)(O)(Cl)Cl","Cl[Co@OH24](N)(N)(O)(Cl)Cl",
                      "Cl[Co@OH25](N)(N)(O)(Cl)Cl","Cl[Co@OH26](N)(N)(O)(Cl)Cl",
                      "Cl[Co@OH29](N)(N)(O)(Cl)Cl","Cl[Co@OH30](N)(N)(O)(Cl)Cl")
        identical3 = ("Cl[Co@OH8](N)(N)(O)(Cl)Cl","Cl[Co@OH9](N)(N)(O)(Cl)Cl",
                      "Cl[Co@OH10](N)(N)(O)(Cl)Cl","Cl[Co@OH11](N)(N)(O)(Cl)Cl",
                      "Cl[Co@OH12](N)(N)(O)(Cl)Cl","Cl[Co@OH13](N)(N)(O)(Cl)Cl")

        mols1 = [rdkit.Chem.MolFromSmiles(i, sanitize=False)
                 for i in identical1]
        mols1 = [rdkit.Chem.AddHs(i, explicitOnly=True) for i in mols1]
        molgraphs1 = [self._TestClass.from_rdmol(i) for i in mols1]
        mols2 = [rdkit.Chem.MolFromSmiles(i, sanitize=False)
                 for i in identical2]
        mols2 = [rdkit.Chem.AddHs(i, explicitOnly=True) for i in mols2]
        molgraphs2 = [self._TestClass.from_rdmol(i) for i in mols2]
        mols3 = [rdkit.Chem.MolFromSmiles(i, sanitize=False)
                 for i in identical3]
        mols3 = [rdkit.Chem.AddHs(i, explicitOnly=True) for i in mols3]
        molgraphs3 = [self._TestClass.from_rdmol(i) for i in mols3]

        assert all(molgraph == molgraphs1[0] for molgraph in molgraphs1)
        assert all(molgraph == molgraphs2[0] for molgraph in molgraphs2)
        assert all(molgraph == molgraphs3[0] for molgraph in molgraphs3)
        assert molgraphs1[0] != molgraphs2[0]
        assert molgraphs1[0] != molgraphs3[0]
        assert molgraphs2[0] != molgraphs3[0]

        set1 = {rdkit.Chem.MolToSmiles(molgraph._to_rdmol()[0],
                                       canonical=True,
                                       ignoreAtomMapNumbers=True,
                                       isomericSmiles=True)
                   for molgraph in molgraphs1}
        set2 = {rdkit.Chem.MolToSmiles(molgraph._to_rdmol()[0],
                                       canonical=True,
                                       ignoreAtomMapNumbers=True,
                                       isomericSmiles=True)
                   for molgraph in molgraphs2}
        set3 = {rdkit.Chem.MolToSmiles(molgraph._to_rdmol()[0],
                                       canonical=True,
                                       ignoreAtomMapNumbers=True,
                                       isomericSmiles=True)
                   for molgraph in molgraphs3}

        assert set() == set1 & set2 == set1 & set3 == set2 & set3

    def test_from_atrop(self):
        self._TestClass()
        smg = StereoMolGraph()
        smg.add_atom(0, "H")
        smg.add_atom(1, "Cl")
        smg.add_atom(2, "C")
        smg.add_atom(3, "C")
        smg.add_atom(4, "F")
        smg.add_atom(5, "I")

        smg.add_bond(0,2)
        smg.add_bond(1,2)
        smg.add_bond(2,3)
        smg.add_bond(3,4)
        smg.add_bond(3,5)

        smg.set_bond_stereo(AtropBond(atoms=(0, 1, 2, 3, 4, 5), parity=1))

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
        assert (
            combined.stereo
            == chiral_product_graph1.stereo
            | chiral_product_graph2.stereo
        )

    def test_atom_stereo_is_isomorphic(self, chiral_product_graph1):
        isomorphic_graph = chiral_product_graph1.copy()
        isomorphic_graph.relabel_atoms(
            {0: 1, 1: 0, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}
        )
        assert chiral_product_graph1.is_isomorphic(isomorphic_graph)

    def test_get_atom_stereo_isomorphic_mappings(self, chiral_product_graph1):
        isomorphic_graph = chiral_product_graph1.copy()
        isomorphic_graph.relabel_atoms(
            {0: 1, 1: 0, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}
        )
        assert {
            frozenset(mapping.items()) for mapping in
            chiral_product_graph1.get_isomorphic_mappings(isomorphic_graph)
         } == {frozenset(mapping.items()) for mapping in (
            {
                0: 0,
                1: 1,
                2: 2,
                3: 3,
                4: 4,
                5: 5,
                6: 6,
                7: 7,
                8: 8,
                9: 9,
                10: 10,
                11: 11,
                12: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                1: 1,
                2: 2,
                3: 3,
                4: 4,
                5: 5,
                6: 6,
                7: 7,
                8: 8,
                9: 9,
                11: 10,
                12: 11,
                10: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                1: 1,
                2: 2,
                3: 3,
                4: 4,
                5: 5,
                6: 6,
                7: 7,
                8: 8,
                9: 9,
                12: 10,
                10: 11,
                11: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                1: 1,
                2: 2,
                3: 3,
                4: 4,
                5: 5,
                7: 6,
                8: 7,
                6: 8,
                9: 9,
                10: 10,
                11: 11,
                12: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                1: 1,
                2: 2,
                3: 3,
                4: 4,
                5: 5,
                7: 6,
                8: 7,
                6: 8,
                9: 9,
                11: 10,
                12: 11,
                10: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                1: 1,
                2: 2,
                3: 3,
                4: 4,
                5: 5,
                7: 6,
                8: 7,
                6: 8,
                9: 9,
                12: 10,
                10: 11,
                11: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                1: 1,
                2: 2,
                3: 3,
                4: 4,
                5: 5,
                8: 6,
                6: 7,
                7: 8,
                9: 9,
                10: 10,
                11: 11,
                12: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                1: 1,
                2: 2,
                3: 3,
                4: 4,
                5: 5,
                8: 6,
                6: 7,
                7: 8,
                9: 9,
                11: 10,
                12: 11,
                10: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                1: 1,
                2: 2,
                3: 3,
                4: 4,
                5: 5,
                8: 6,
                6: 7,
                7: 8,
                9: 9,
                12: 10,
                10: 11,
                11: 12,
                13: 13,
                14: 14,
            },
        )}

    def test_get_isomorphic_mappings_of_enantiomers(
        self, enantiomer_graph1, enantiomer_graph2
    ):
        assert not enantiomer_graph1.is_isomorphic(enantiomer_graph2)

    def test_get_isomorphic_mappings(
        self, chiral_product_graph1, chiral_product_graph2
    ):
        expected_isomorphisms = (
            {
                0: 0,
                2: 1,
                1: 2,
                4: 3,
                3: 4,
                9: 5,
                10: 6,
                11: 7,
                12: 8,
                5: 9,
                6: 10,
                7: 11,
                8: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                2: 1,
                1: 2,
                4: 3,
                3: 4,
                9: 5,
                10: 6,
                11: 7,
                12: 8,
                5: 9,
                7: 10,
                8: 11,
                6: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                2: 1,
                1: 2,
                4: 3,
                3: 4,
                9: 5,
                10: 6,
                11: 7,
                12: 8,
                5: 9,
                8: 10,
                6: 11,
                7: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                2: 1,
                1: 2,
                4: 3,
                3: 4,
                9: 5,
                11: 6,
                12: 7,
                10: 8,
                5: 9,
                6: 10,
                7: 11,
                8: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                2: 1,
                1: 2,
                4: 3,
                3: 4,
                9: 5,
                11: 6,
                12: 7,
                10: 8,
                5: 9,
                7: 10,
                8: 11,
                6: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                2: 1,
                1: 2,
                4: 3,
                3: 4,
                9: 5,
                11: 6,
                12: 7,
                10: 8,
                5: 9,
                8: 10,
                6: 11,
                7: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                2: 1,
                1: 2,
                4: 3,
                3: 4,
                9: 5,
                12: 6,
                10: 7,
                11: 8,
                5: 9,
                6: 10,
                7: 11,
                8: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                2: 1,
                1: 2,
                4: 3,
                3: 4,
                9: 5,
                12: 6,
                10: 7,
                11: 8,
                5: 9,
                7: 10,
                8: 11,
                6: 12,
                13: 13,
                14: 14,
            },
            {
                0: 0,
                2: 1,
                1: 2,
                4: 3,
                3: 4,
                9: 5,
                12: 6,
                10: 7,
                11: 8,
                5: 9,
                8: 10,
                6: 11,
                7: 12,
                13: 13,
                14: 14,
            },
        )

        assert chiral_product_graph1.is_isomorphic(chiral_product_graph2)
        assert chiral_product_graph1.is_isomorphic(chiral_product_graph2)
        assert all(
            i in expected_isomorphisms
            for i in chiral_product_graph1.get_isomorphic_mappings(
                chiral_product_graph2
            )
        )

    def test_enantiomer(self, enantiomer_graph1, enantiomer_graph2):
        assert enantiomer_graph1.enantiomer().is_isomorphic(enantiomer_graph2)
        assert not enantiomer_graph1.is_isomorphic(enantiomer_graph2)

    def test_hash_enantiomers(self, enantiomer_graph1, enantiomer_graph2):
        assert enantiomer_graph1._atom_stereo != enantiomer_graph2._atom_stereo
        #raise Exception(enantiomer_graph1._atom_stereo, enantiomer_graph2._atom_stereo)
        assert (hash(enantiomer_graph1)
                != hash(enantiomer_graph2))

    def test_valid_stereo(self, chiral_product_graph1):
        assert chiral_product_graph1.is_stereo_valid()


class TestStereoCondensedReactionGraph(
    TestStereoMolGraph, TestCondensedReactionGraph
):
    _TestClass: type[StereoCondensedReactionGraph] = (
        StereoCondensedReactionGraph
    )

    @pytest.fixture
    def chiral_ts_geo1(self, data_path):
        return Geometry.from_xyz_file(
            data_path/ "conrot_reaction/ts.xyz")

    @pytest.fixture
    def chiral_ts_geo2(self, data_path):
        return Geometry.from_xyz_file(
            data_path / "disrot_reaction/ts.xyz")
        

    def test_product_with_attributes(self, crg):
        super().test_product_with_attributes(crg)
        expected_mol_atom_stereo = {
            key: value
            for key, value in crg._atom_stereo.items()
            if value in (None, Change.FORMED)
        }
        assert (
            crg.product(keep_attributes=True)._atom_stereo
            == expected_mol_atom_stereo
        )

    def test_product_without_attributes(self, crg):
        super().test_product_without_attributes(crg)
        expected_mol_atom_stereo = {
            key: value
            for key, value in crg._atom_stereo.items()
            if value in (None, Change.FORMED)
        }
        assert (
            crg.product(keep_attributes=True)._atom_stereo
            == expected_mol_atom_stereo
        )

    def test_reactant_with_attributes(self, crg):
        super().test_reactant_with_attributes(crg)
        expected_mol_atom_stereo = {
            key: value
            for key, value in crg._atom_stereo.items()
            if value in (None, Change.BROKEN)
        }
        assert (
            crg.product(keep_attributes=True)._atom_stereo
            == expected_mol_atom_stereo
        )

    def test_reactant_without_attributes(self, crg):
        super().test_reactant_without_attributes(crg)
        expected_mol_atom_stereo = {
            key: value
            for key, value in crg._atom_stereo.items()
            if value in (None, Change.BROKEN)
        }
        assert (
            crg.product(keep_attributes=True)._atom_stereo
            == expected_mol_atom_stereo
        )

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
        assert (
            combined.stereo
            == chiral_product_graph1.stereo
            | chiral_product_graph2.stereo
        )
        assert (
            {**combined.atom_stereo_changes, **combined.bond_stereo_changes}
            == {**chiral_product_graph1.atom_stereo_changes,
                **chiral_product_graph1.bond_stereo_changes,
                **chiral_product_graph2.atom_stereo_changes,
                **chiral_product_graph2.bond_stereo_changes}
        )

    def test_from_chain_of_states_reaction(self, data_path):
        reactant_geo = Geometry.from_xyz_file(
            data_path / "methylamine_phosgenation_trans_r.xyz"
            
        )
        product_geo = Geometry.from_xyz_file(
            data_path / "methylamine_phosgenation_trans_p.xyz"
        )
        ts_geo =  Geometry.from_xyz_file(
            data_path / "methylamine_phosgenation_trans_ts.xyz"
        )
        
        scrg = self._TestClass.from_reactant_product_and_ts_geometry(
            reactant_geo, product_geo, ts_geo
        )

        atoms = {
            0: {"atom_type": PTOE["H"]},
            1: {"atom_type": PTOE["N"]},
            2: {"atom_type": PTOE["H"]},
            3: {"atom_type": PTOE["C"]},
            4: {"atom_type": PTOE["Cl"]},
            5: {"atom_type": PTOE["O"]},
            6: {"atom_type": PTOE["Cl"]},
            7: {"atom_type": PTOE["C"]},
            8: {"atom_type": PTOE["H"]},
            9: {"atom_type": PTOE["H"]},
            10: {"atom_type": PTOE["H"]},
        }
        bonds = {
            Bond((0, 1)): {},
            Bond((1, 2)): {"reaction": Change.BROKEN},
            Bond((1, 7)): {},
            Bond((1, 3)): {},
            Bond((2, 6)): {"reaction": Change.FORMED},
            Bond((3, 4)): {},
            Bond((3, 6)): {"reaction": Change.BROKEN},
            Bond((3, 5)): {},
            Bond((7, 10)): {},
            Bond((7, 9)): {},
            Bond((7, 8)): {},
        }
        stereo = {7: Tetrahedral((7, 1, 8, 9, 10), -1)}
        atom_stereo_change = defaultdict(ChangeDict,{
            1: ChangeDict({
                Change.BROKEN: Tetrahedral(
                    (1, 0, 2, 3, 7), -1
                )
            }),
            3: ChangeDict({
                Change.BROKEN: Tetrahedral(
                    (3, 1, 4, 5, 6), -1
                )
            }),})
        bond_stereo_change = defaultdict(ChangeDict,{
            Bond({1, 3}): ChangeDict({
                Change.FORMED: PlanarBond(
                    (4, 5, 3, 1, 0, 7), 0
                )})})
        assert scrg.atoms_with_attributes == atoms
        assert scrg.bonds_with_attributes == bonds
        assert scrg.stereo == stereo
        assert scrg._atom_stereo_change == atom_stereo_change
        assert scrg._bond_stereo_change == bond_stereo_change

    @pytest.fixture
    def scrg_stereo_inversion(self, data_path):
        """ A SCRG with just the inversion of a stereocenter, reactant and product are tetrahedral but
        the transition state is square planar. """
        reactant_geo = Geometry.from_xyz_file(data_path / "fluoro_chloro_bromomethane_r.xyz")
        product_geo = Geometry.from_xyz_file(data_path / "fluoro_chloro_bromomethane_s.xyz")
        ts_geo =  Geometry.from_xyz_file(data_path / "fluoro_chloro_bromomethane_ts.xyz")
        
        scrg = self._TestClass.from_reactant_product_and_ts_geometry(
            reactant_geo, product_geo, ts_geo
        )
        return scrg

    def test_creation_from_xyz_atom_stereocenter_inversion(
        self, scrg_stereo_inversion
    ):
        scrg = scrg_stereo_inversion
        assert scrg.get_atom_stereo_change(0) == {
            Change.BROKEN: Tetrahedral(
                (0, 1, 2, 3, 4), 1
            ),
            Change.FORMED: Tetrahedral(
                (0, 1, 2, 3, 4), -1
            ),
        }
        assert scrg._bond_stereo_change == {}

    def test_relabel_reaction_atoms(self, scrg_stereo_inversion):
        scrg = scrg_stereo_inversion
        scrg.relabel_atoms({0: 11, 1: 10, 2: 20, 3: 30, 4: 40}, copy=False)
        assert scrg.get_atom_stereo_change(11) == {
            Change.BROKEN: Tetrahedral(
                (11, 10, 20, 30, 40), 1
            ),
            Change.FORMED: Tetrahedral(
                (11, 10, 20, 30, 40), -1
            ),
        }
        assert scrg._bond_stereo_change == {}

    def test_relabel_reaction_atoms_copy(self, scrg_stereo_inversion):
        scrg = scrg_stereo_inversion
        new_scrg = scrg.relabel_atoms(
            {0: 11, 1: 10, 2: 20, 3: 30, 4: 40}, copy=True
        )
        assert new_scrg.get_atom_stereo_change(11) == {
            Change.BROKEN: Tetrahedral(
                (11, 10, 20, 30, 40), 1
            ),
            Change.FORMED: Tetrahedral(
                (11, 10, 20, 30, 40), -1
            ),
        }
        assert scrg._bond_stereo_change == {}

    @pytest.fixture
    def chiral_reaction_chiral_ts_scrg1(
        self, chiral_reactant_geo, chiral_product_geo1, chiral_ts_geo1
    ):
        return self._TestClass.from_reactant_product_and_ts_geometry(
            chiral_reactant_geo, chiral_product_geo1, chiral_ts_geo1
        )

    @pytest.fixture
    def chiral_reaction_chiral_ts_scrg2(
        self, chiral_reactant_geo, chiral_product_geo2, chiral_ts_geo2
    ):
        return self._TestClass.from_reactant_product_and_ts_geometry(
            chiral_reactant_geo, chiral_product_geo2, chiral_ts_geo2
        )

    def test_isomorphism_same_reactant_and_product_without_ts(
        self, chiral_reaction_scrg1, chiral_reaction_scrg2
    ):
        assert chiral_reaction_scrg1.product().is_isomorphic(
            chiral_reaction_scrg2.product()
        )
        assert chiral_reaction_scrg1.reactant().is_isomorphic(
            chiral_reaction_scrg2.reactant()
        )
        assert not chiral_reaction_scrg1.is_isomorphic(chiral_reaction_scrg2)

    def test_isomorphism_same_reactant_and_product_but_different_ts(
        self, chiral_reaction_chiral_ts_scrg1, chiral_reaction_chiral_ts_scrg2
    ):
        assert chiral_reaction_chiral_ts_scrg1.reactant().is_isomorphic(
            chiral_reaction_chiral_ts_scrg2.reactant()
        )
        assert chiral_reaction_chiral_ts_scrg1.product().is_isomorphic(
            chiral_reaction_chiral_ts_scrg2.product()
        )

        assert not chiral_reaction_chiral_ts_scrg1.is_isomorphic(
            chiral_reaction_chiral_ts_scrg2
        )

    def test_reverse_reaction(self, chiral_reaction_scrg1):
        reversed_reaction = chiral_reaction_scrg1.reverse_reaction()
        assert (
            reversed_reaction.get_broken_bonds()
            == chiral_reaction_scrg1.get_formed_bonds()
        )
        assert (
            reversed_reaction.get_formed_bonds()
            == chiral_reaction_scrg1.get_broken_bonds()
        )
        
        double_reverset_reaction = reversed_reaction.reverse_reaction()
        assert (chiral_reaction_scrg1.atom_stereo == double_reverset_reaction.atom_stereo)
        assert (chiral_reaction_scrg1.bond_stereo == double_reverset_reaction.bond_stereo)
        assert (chiral_reaction_scrg1.bond_stereo_changes == double_reverset_reaction.bond_stereo_changes)
        assert (chiral_reaction_scrg1.atom_stereo_changes == double_reverset_reaction.atom_stereo_changes)

        assert double_reverset_reaction == chiral_reaction_scrg1

    def test_hash_stereo_reaction(
        self, chiral_reaction_scrg1, chiral_reaction_scrg2
    ):
        assert (
            hash(chiral_reaction_scrg1)
            != hash(chiral_reaction_scrg2)
        )

    def test_hash_stereo_reaction_with_ts(
        self, chiral_reaction_chiral_ts_scrg1, chiral_reaction_chiral_ts_scrg2
    ):
        assert (
            hash(chiral_reaction_chiral_ts_scrg1)
            != hash(chiral_reaction_chiral_ts_scrg2)
        )

    def test_hash_enantiomers(
        self, enantiomer_graph1, enantiomer_graph2
    ):
        assert (hash(enantiomer_graph1)
                != hash(enantiomer_graph2))
