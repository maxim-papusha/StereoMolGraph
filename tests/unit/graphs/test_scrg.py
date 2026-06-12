from collections import defaultdict

import pytest
import rdkit.Chem
import rdkit.Chem.rdchem

from stereomolgraph import (
    Bond,
    StereoCondensedReactionGraph,
)
from stereomolgraph.coords import Geometry
from stereomolgraph.graphs.crg import Change
from stereomolgraph.graphs.scrg import ChangeDict
from stereomolgraph.ipython import View2D
from stereomolgraph.periodic_table import PERIODIC_TABLE as PTOE
from stereomolgraph.stereodescriptors import (
    PlanarBond,
    SquarePlanar,
    Tetrahedral,
)

from .test_crg import TestCondensedReactionGraph
from .test_smg import TestStereoMolGraph


class TestStereoCondensedReactionGraph(TestStereoMolGraph, TestCondensedReactionGraph):
    _TestClass: type[StereoCondensedReactionGraph] = StereoCondensedReactionGraph

    @pytest.fixture
    def chiral_ts_geo1(self, data_path):
        return Geometry.from_xyz_file(data_path / "conrot_reaction/ts.xyz")

    @pytest.fixture
    def chiral_ts_geo2(self, data_path):
        return Geometry.from_xyz_file(data_path / "disrot_reaction/ts.xyz")

    def test_product_with_attributes(self, crg):
        super().test_product_with_attributes(crg)
        expected_mol_atom_stereo = {
            key: value
            for key, value in crg._atom_stereo.items()
            if value in (None, Change.FORMED)
        }
        assert (
            crg.product(keep_attributes=True)._atom_stereo == expected_mol_atom_stereo
        )

    def test_product_without_attributes(self, crg):
        super().test_product_without_attributes(crg)
        expected_mol_atom_stereo = {
            key: value
            for key, value in crg._atom_stereo.items()
            if value in (None, Change.FORMED)
        }
        assert (
            crg.product(keep_attributes=True)._atom_stereo == expected_mol_atom_stereo
        )

    def test_reactant_with_attributes(self, crg):
        super().test_reactant_with_attributes(crg)
        expected_mol_atom_stereo = {
            key: value
            for key, value in crg._atom_stereo.items()
            if value in (None, Change.BROKEN)
        }
        assert (
            crg.product(keep_attributes=True)._atom_stereo == expected_mol_atom_stereo
        )

    def test_reactant_without_attributes(self, crg):
        super().test_reactant_without_attributes(crg)
        expected_mol_atom_stereo = {
            key: value
            for key, value in crg._atom_stereo.items()
            if value in (None, Change.BROKEN)
        }
        assert (
            crg.product(keep_attributes=True)._atom_stereo == expected_mol_atom_stereo
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
            == chiral_product_graph1.stereo | chiral_product_graph2.stereo
        )
        assert {
            **combined.atom_stereo_changes,
            **combined.bond_stereo_changes,
        } == {
            **chiral_product_graph1.atom_stereo_changes,
            **chiral_product_graph1.bond_stereo_changes,
            **chiral_product_graph2.atom_stereo_changes,
            **chiral_product_graph2.bond_stereo_changes,
        }

    def test_from_chain_of_states_reaction(self, data_path):
        reactant_geo = Geometry.from_xyz_file(
            data_path / "methylamine_phosgenation_trans_r.xyz"
        )
        product_geo = Geometry.from_xyz_file(
            data_path / "methylamine_phosgenation_trans_p.xyz"
        )
        ts_geo = Geometry.from_xyz_file(
            data_path / "methylamine_phosgenation_trans_ts.xyz"
        )

        scrg = self._TestClass.from_geometries(reactant_geo, product_geo, ts_geo)

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
        stereo = {7: Tetrahedral((7, 1, 8, 9, 10), 1)}
        atom_stereo_change = defaultdict(
            ChangeDict,
            {
                1: ChangeDict({Change.BROKEN: Tetrahedral((1, 0, 2, 3, 7), 1)}),
                3: ChangeDict({Change.BROKEN: Tetrahedral((3, 1, 4, 5, 6), 1)}),
            },
        )
        bond_stereo_change = defaultdict(
            ChangeDict,
            {
                Bond({1, 3}): ChangeDict(
                    {Change.FORMED: PlanarBond((4, 5, 3, 1, 0, 7), 0)}
                )
            },
        )
        assert scrg.atoms_with_attributes == atoms
        assert scrg.bonds_with_attributes == bonds
        assert scrg.stereo == stereo
        assert scrg._atom_stereo_change == atom_stereo_change
        assert scrg._bond_stereo_change == bond_stereo_change

    @pytest.fixture
    def scrg_stereo_inversion(self, data_path):
        """A SCRG with just the inversion of a stereocenter, reactant and product are tetrahedral but
        the transition state is square planar."""
        reactant_geo = Geometry.from_xyz_file(
            data_path / "fluoro_chloro_bromomethane_r.xyz"
        )
        product_geo = Geometry.from_xyz_file(
            data_path / "fluoro_chloro_bromomethane_s.xyz"
        )
        ts_geo = Geometry.from_xyz_file(data_path / "fluoro_chloro_bromomethane_ts.xyz")

        scrg = self._TestClass.from_geometries(reactant_geo, product_geo, ts_geo)
        return scrg

    def test_creation_from_xyz_atom_stereocenter_inversion(self, scrg_stereo_inversion):
        scrg = scrg_stereo_inversion
        assert scrg.get_atom_stereo_change(0) == {
            Change.BROKEN: Tetrahedral((0, 1, 2, 3, 4), -1),
            Change.FLEETING: SquarePlanar((0, 4, 2, 3, 1), 0),
            Change.FORMED: Tetrahedral((0, 1, 2, 3, 4), 1),
        }
        assert scrg._bond_stereo_change == {}

    def test_relabel_reaction_atoms(self, scrg_stereo_inversion):
        scrg = scrg_stereo_inversion
        scrg.relabel_atoms({0: 11, 1: 10, 2: 20, 3: 30, 4: 40}, copy=False)
        assert scrg.get_atom_stereo_change(11) == {
            Change.BROKEN: Tetrahedral((11, 10, 20, 30, 40), -1),
            Change.FLEETING: SquarePlanar((11, 40, 20, 30, 10), 0),
            Change.FORMED: Tetrahedral((11, 10, 20, 30, 40), 1),
        }
        assert scrg._bond_stereo_change == {}

    def test_relabel_reaction_atoms_copy(self, scrg_stereo_inversion):
        scrg = scrg_stereo_inversion
        new_scrg = scrg.relabel_atoms({0: 11, 1: 10, 2: 20, 3: 30, 4: 40}, copy=True)
        assert new_scrg.get_atom_stereo_change(11) == {
            Change.BROKEN: Tetrahedral((11, 10, 20, 30, 40), -1),
            Change.FLEETING: SquarePlanar((11, 40, 20, 30, 10), 0),
            Change.FORMED: Tetrahedral((11, 10, 20, 30, 40), 1),
        }
        assert scrg._bond_stereo_change == {}

    @pytest.fixture
    def chiral_reaction_chiral_ts_scrg1(
        self, chiral_reactant_geo, chiral_product_geo1, chiral_ts_geo1
    ):
        return self._TestClass.from_geometries(
            chiral_reactant_geo, chiral_product_geo1, chiral_ts_geo1
        )

    @pytest.fixture
    def chiral_reaction_chiral_ts_scrg2(
        self, chiral_reactant_geo, chiral_product_geo2, chiral_ts_geo2
    ):
        return self._TestClass.from_geometries(
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
        assert chiral_reaction_scrg1.atom_stereo == double_reverset_reaction.atom_stereo
        assert chiral_reaction_scrg1.bond_stereo == double_reverset_reaction.bond_stereo
        assert (
            chiral_reaction_scrg1.bond_stereo_changes
            == double_reverset_reaction.bond_stereo_changes
        )
        assert (
            chiral_reaction_scrg1.atom_stereo_changes
            == double_reverset_reaction.atom_stereo_changes
        )

        assert double_reverset_reaction == chiral_reaction_scrg1

    def test_hash_stereo_reaction(self, chiral_reaction_scrg1, chiral_reaction_scrg2):
        assert hash(chiral_reaction_scrg1.copy(frozen=True)) == hash(
            chiral_reaction_scrg2.copy(frozen=True)
        )

    def test_ts_uses_single_defined_stereo_change(self):
        scrg = self._TestClass()
        scrg.add_atom(0, "C")
        scrg.add_atom(1, "F")
        scrg.add_atom(2, "Br")
        scrg.add_atom(3, "Cl")
        scrg.add_atom(4, "H")
        scrg.add_bond(0, 1)
        scrg.add_bond(0, 2)
        scrg.add_bond(0, 3)
        scrg.add_broken_bond(0, 4)
        scrg.set_atom_stereo_change(broken=Tetrahedral((0, 1, 2, 3, 4), -1))

        ts = scrg.ts()
        strict_ts = scrg.ts(infer_non_fleeting_stereo=False)

        assert ts.get_atom_stereo(0) == Tetrahedral((0, 1, 2, 3, 4), -1)
        assert strict_ts.get_atom_stereo(0) is None

    def test_view2d_scrg_uses_dotted_reaction_bonds_and_highlights_bond_stereo(self):
        scrg = self._TestClass()
        scrg.add_atom(0, "C")
        scrg.add_atom(1, "H")
        scrg.add_atom(2, "H")
        scrg.add_atom(3, "C")
        scrg.add_atom(4, "H")
        scrg.add_atom(5, "H")
        scrg.add_atom(6, "H")
        scrg.add_atom(7, "C")
        scrg.add_atom(8, "H")
        scrg.add_atom(9, "H")
        scrg.add_atom(10, "C")
        scrg.add_atom(11, "H")
        scrg.add_atom(12, "H")
        scrg.add_atom(13, "H")

        scrg.add_bond(0, 1)
        scrg.add_bond(0, 2)
        scrg.add_bond(0, 3)
        scrg.add_bond(3, 4)
        scrg.add_bond(3, 5)
        scrg.add_formed_bond(3, 6)
        scrg.set_bond_stereo_change(broken=PlanarBond((1, 2, 0, 3, 4, 5), 0))

        scrg.add_bond(7, 8)
        scrg.add_bond(7, 9)
        scrg.add_bond(7, 10)
        scrg.add_bond(10, 11)
        scrg.add_bond(10, 12)
        scrg.add_broken_bond(10, 13)
        scrg.set_bond_stereo_change(formed=PlanarBond((8, 9, 7, 10, 11, 12), 0))

        mol, ht = View2D(generate_bond_orders=False)._to_mol(scrg)

        broken_stereo_idx = mol.GetBondBetweenAtoms(0, 3).GetIdx()
        formed_reaction_idx = mol.GetBondBetweenAtoms(3, 6).GetIdx()
        formed_stereo_idx = mol.GetBondBetweenAtoms(7, 10).GetIdx()
        broken_reaction_idx = mol.GetBondBetweenAtoms(10, 13).GetIdx()

        assert (
            mol.GetBondWithIdx(broken_stereo_idx).GetBondType()
            == rdkit.Chem.rdchem.BondType.AROMATIC
        )
        assert (
            mol.GetBondWithIdx(formed_stereo_idx).GetBondType()
            == rdkit.Chem.rdchem.BondType.AROMATIC
        )
        assert (
            mol.GetBondWithIdx(formed_reaction_idx).GetBondType()
            == rdkit.Chem.rdchem.BondType.HYDROGEN
        )
        assert (
            mol.GetBondWithIdx(broken_reaction_idx).GetBondType()
            == rdkit.Chem.rdchem.BondType.HYDROGEN
        )

        assert ht.highlight_bond_colors[broken_stereo_idx] == (1, 0, 0)
        assert ht.highlight_bond_colors[formed_stereo_idx] == (0, 0, 1)
        assert ht.highlight_bond_colors[formed_reaction_idx] == (0, 0, 1)
        assert ht.highlight_bond_colors[broken_reaction_idx] == (1, 0, 0)

        svg = View2D(generate_bond_orders=False).svg(scrg)
        formed_reaction_lines = [
            line
            for line in svg.splitlines()
            if f"class='bond-{formed_reaction_idx} " in line
        ]
        broken_reaction_lines = [
            line
            for line in svg.splitlines()
            if f"class='bond-{broken_reaction_idx} " in line
        ]

        assert any(
            "stroke:#0000FF" in line and "stroke-dasharray" in line
            for line in formed_reaction_lines
        )
        assert any(
            "stroke:#FF0000" in line and "stroke-dasharray" in line
            for line in broken_reaction_lines
        )

    def test_hash_stereo_reaction_with_ts(
        self, chiral_reaction_chiral_ts_scrg1, chiral_reaction_chiral_ts_scrg2
    ):
        assert hash(chiral_reaction_chiral_ts_scrg1.copy(frozen=True)) != hash(
            chiral_reaction_chiral_ts_scrg2.copy(frozen=True)
        )

    def test_hash_enantiomers(self, enantiomer_graph1, enantiomer_graph2):
        assert hash(enantiomer_graph1.copy(frozen=True)) != hash(
            enantiomer_graph2.copy(frozen=True)
        )
