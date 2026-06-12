import pytest

from stereomolgraph import (
    Bond,
    CondensedReactionGraph,
)
from stereomolgraph.graphs.crg import Change

from .test_mg import TestMolGraph


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

    def test_bonds(
        self,
        chiral_reactant_graph,
        chiral_product_graph1,
        chiral_reaction_scrg1,
    ):
        unchanged_b = {
            bond
            for bond in chiral_reaction_scrg1.bonds
            if bond not in chiral_reaction_scrg1.get_broken_bonds()
            and bond not in chiral_reaction_scrg1.get_formed_bonds()
        }

        formed_b = chiral_reaction_scrg1.get_formed_bonds()
        broken_b = chiral_reaction_scrg1.get_broken_bonds()

        assert set(chiral_reactant_graph.bonds) == unchanged_b | broken_b
        assert set(chiral_product_graph1.bonds) == unchanged_b | formed_b

    def test_add_bond_error(self, crg):
        with pytest.raises(TypeError):
            crg.add_bond(0, 1, reaction="test")

    def test_add_bond_with_reaction_attr(self, crg):
        crg.add_bond(0, 1, reaction=Change.FORMED, bond_order=1)
        assert crg.get_bond_attribute(0, 1, attr="reaction") == Change.FORMED
        assert crg.get_bond_attribute(0, 1, attr="bond_order") == 1

    def test_set_bond_attribute_reaction_exception(self, crg):
        with pytest.raises(ValueError):
            crg.set_bond_attribute(0, 1, attr="reaction", value="test")

    def test_set_bond_attribute_reaction(self, crg):
        crg.set_bond_attribute(0, 1, attr="reaction", value=Change.FORMED)
        assert crg.get_bond_attribute(0, 1, attr="reaction") == Change.FORMED

    def test_set_bond_attribute(self, crg):
        crg.set_bond_attribute(0, 1, attr="bond_label", value="test")
        assert crg.get_bond_attribute(0, 1, attr="bond_label") == "test"

    def test_add_formed_bond(self, crg):
        crg.add_formed_bond(1, 2)
        assert crg.get_bond_attribute(1, 2, attr="reaction") == Change.FORMED

    def test_add_broken_bond(self, crg):
        crg.add_broken_bond(1, 2)
        assert crg.get_bond_attribute(2, 1, attr="reaction") == Change.BROKEN

    def test_get_formed_bonds(self, crg):
        assert set(crg.get_formed_bonds()) == {
            Bond((0, 1)),
        }

    def test_get_broken_bonds(self, crg):
        assert set(crg.get_broken_bonds()) == {
            Bond((2, 3)),
        }

    def test_active_atoms_crg(self, crg):
        assert set(crg.active_atoms(additional_layer=0)) == {0, 1, 2, 3}
        assert set(crg.active_atoms(additional_layer=1)) == {0, 1, 2, 3}

        crg.add_atom(4, "O")
        crg.add_bond(4, 2)
        crg.add_atom(5, "H")
        crg.add_bond(5, 4)

        assert set(crg.active_atoms(additional_layer=0)) == {0, 1, 2, 3}
        assert set(crg.active_atoms(additional_layer=1)) == {0, 1, 2, 3, 4}

    def test_active_atoms(self):
        crg = self._TestClass()
        crg.add_atom(1, "C")
        crg.add_atom(0, "C")
        crg.add_broken_bond(0, 1)
        for i in range(2, 10):
            crg.add_atom(i, "C")
            crg.add_bond(i, i - 1)
        for j in range(0, 5):
            assert set(crg.active_atoms(additional_layer=j)) == {*range(j + 2)}

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
                atom_attr for atom_attr in attr_dict.keys() if atom_attr != "atom_type"
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
                atom_attr for atom_attr in attr_dict.keys() if atom_attr != "atom_type"
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
        return self._TestClass.from_geometries(chiral_reactant_geo, chiral_product_geo1)

    @pytest.fixture
    def chiral_reaction_scrg2(self, chiral_reactant_geo, chiral_product_geo2):
        return self._TestClass.from_geometries(chiral_reactant_geo, chiral_product_geo2)

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
