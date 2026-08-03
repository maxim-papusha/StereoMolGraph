import pytest
import rdkit.Chem
import rdkit.Chem.rdchem

from stereomolgraph import (
    Bond,
    StereoMolGraph,
)
from stereomolgraph.graph2rdmol import stereo_mol_graph_to_rdmol
from stereomolgraph.periodic_table import PERIODIC_TABLE as PTOE
from stereomolgraph.stereodescriptors import (
    AtropBond,
    PlanarBond,
    Tetrahedral,
)

from .test_mg import TestMolGraph


class TestStereoMolGraph(TestMolGraph):
    _TestClass: type[StereoMolGraph] = StereoMolGraph

    def test_from_geometries1(self, chiral_reactant_graph):
        graph = chiral_reactant_graph
        expected_atom_stereo = {
            1: Tetrahedral((1, 0, 2, 3, 9), -1),
            0: Tetrahedral((0, 1, 2, 13, 14), 1),
            9: Tetrahedral((9, 1, 10, 11, 12), 1),
            5: Tetrahedral((5, 2, 6, 7, 8), 1),
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
        expected_bonds = {
            Bond(pair)
            for pair in [
                (0, 1),
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
                (9, 12),
            ]
        }
        assert graph.atom_types == tuple(expected_atom_types)
        assert set(graph.bonds) == expected_bonds
        assert all(
            graph.get_atom_stereo(key) == value
            for key, value in expected_atom_stereo.items()
        )

    def test_atom_stereo(self, chiral_product_graph1):
        expected = {
            1: Tetrahedral((1, 0, 3, 9, 13), -1),
            5: Tetrahedral((5, 2, 6, 7, 8), 1),
            9: Tetrahedral((9, 1, 10, 11, 12), 1),
        }
        expected2 = {Bond((0, 2)): PlanarBond((1, 14, 0, 2, 4, 5), 0)}
        assert all(
            key in chiral_product_graph1._atom_stereo for key in set(expected.keys())
        )
        assert all(
            key in expected for key in set(chiral_product_graph1._atom_stereo.keys())
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
        g1.set_bond_stereo(PlanarBond((0, 1, 2, 3, 4, 5), 0))
        g2.set_bond_stereo(PlanarBond((1, 0, 2, 3, 4, 5), 0))
        g3.set_bond_stereo(PlanarBond((0, 1, 2, 3, 4, 5), None))
        rdmol_g1, idx_atom_map_dict_g1 = stereo_mol_graph_to_rdmol(g1)
        rdmol_g2, idx_atom_map_dict_g2 = stereo_mol_graph_to_rdmol(g2)
        rdmol_g3, idx_atom_map_dict_g3 = stereo_mol_graph_to_rdmol(g3)

        db1 = rdmol_g1.GetBondBetweenAtoms(2, 3)
        stereo_atoms1 = {idx_atom_map_dict_g1[i] for i in db1.GetStereoAtoms()}
        assert stereo_atoms1 == {0, 4} or stereo_atoms1 == {1, 5}
        assert db1.GetStereo() == rdkit.Chem.rdchem.BondStereo.STEREOZ  # type: ignore

        db2 = rdmol_g2.GetBondBetweenAtoms(2, 3)
        stereo_atoms2 = {idx_atom_map_dict_g2[i] for i in db2.GetStereoAtoms()}
        assert stereo_atoms2 == {1, 4} or stereo_atoms2 == {0, 5}
        assert db2.GetStereo() == rdkit.Chem.rdchem.BondStereo.STEREOZ  # type: ignore

        assert g3.get_bond_stereo((2, 3)).parity is None
        db3 = rdmol_g3.GetBondBetweenAtoms(2, 3)

        assert db3.GetStereo() == rdkit.Chem.rdchem.BondStereo.STEREOANY  # type: ignore

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

        mol, _ = stereo_mol_graph_to_rdmol(g)
        chiral_tag = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW  # type: ignore
        assert mol.GetAtomWithIdx(0).GetChiralTag() == chiral_tag

        g.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), -1))
        mol, _ = stereo_mol_graph_to_rdmol(g)
        chiral_tag = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW  # type: ignore
        assert mol.GetAtomWithIdx(0).GetChiralTag() == chiral_tag

    @pytest.mark.parametrize(
        "inchi",
        [
            (
                r"InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6+/m1/s1"
            ),
            (r"InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m0/s1"),
            (r"InChI=1S/C2H2Cl2/c3-1-2-4/h1-2H/b2-1+"),
            (r"InChI=1S/C2H2Cl2/c3-1-2-4/h1-2H/b2-1-"),
            (r"InChI=1S/C4H6/c1-3-4-2/h3-4H,1-2H2"),
        ],
        ids=[
            "alpha-D-gulopyranose",
            "(R)-Bromochlorofluoromethane",
            "Trans-1,2-Dichloroethylene",
            "Cis-1,2-Dichloroethylene",
            "Butadiene",
        ],
    )
    def test_from_rdmol_to_rdmol_stereo(self, inchi):
        pytest.skip("RDKit construction tests moved to test_rdkit_conversion.py")

    @pytest.mark.parametrize(
        "inchi",
        [
            (
                r"InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6+/m1/s1"
            ),
            (r"InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"),
        ],
        ids=[
            "alpha-D-gulopyranose",
            "benzene",
        ],
    )
    def test_from_rdmol_eq_from_geometry(self, inchi):
        pytest.skip("RDKit construction tests moved to test_rdkit_conversion.py")

    @pytest.mark.parametrize(
        "smiles",
        [
            "[H][Pt@SP1](F)(Cl)Br",
            "Cl[Pt@SP1](Cl)([NH3])[NH3]",
            "Cl[Pt@SP2](Cl)([NH3])[NH3]",
        ],
        ids=[
            "Hydridofluorochlorobromoplatinum(II)",
            "(SP-4-2)-diamminedichloroplatinum",
            "(SP-4-1)-diamminedichloroplatinum",
        ],
    )
    def test_from_rdmol_to_rdmol_square_planar(self, smiles):
        pytest.skip("RDKit construction tests moved to test_rdkit_conversion.py")

    def test_from_rdmol_square_planar_different(self):
        pytest.skip("RDKit construction tests moved to test_rdkit_conversion.py")

    def test_from_rdmol_square_planar(self):
        pytest.skip("RDKit construction tests moved to test_rdkit_conversion.py")

    def test_from_rdmol_trigonal_bipyramidal(self):
        pytest.skip("RDKit construction tests moved to test_rdkit_conversion.py")

    def test_from_rdmol_octahedral(self):
        pytest.skip("RDKit construction tests moved to test_rdkit_conversion.py")

    def test_from_rdmol_octahedral_compare(self):
        pytest.skip("RDKit construction tests moved to test_rdkit_conversion.py")

    def test_from_atrop(self):
        self._TestClass()
        smg = StereoMolGraph()
        smg.add_atom(0, "H")
        smg.add_atom(1, "Cl")
        smg.add_atom(2, "C")
        smg.add_atom(3, "C")
        smg.add_atom(4, "F")
        smg.add_atom(5, "I")

        smg.add_bond(0, 2)
        smg.add_bond(1, 2)
        smg.add_bond(2, 3)
        smg.add_bond(3, 4)
        smg.add_bond(3, 5)

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
            == chiral_product_graph1.stereo | chiral_product_graph2.stereo
        )

    def test_atom_stereo_is_isomorphic(self, chiral_product_graph1):
        isomorphic_graph = chiral_product_graph1.copy()
        isomorphic_graph.relabel_atoms({0: 1, 1: 0, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7})
        assert chiral_product_graph1.is_isomorphic(isomorphic_graph)

    def test_isomorphic_enantiomers(self, enantiomer_graph1, enantiomer_graph2):
        assert not enantiomer_graph1.is_isomorphic(enantiomer_graph2)

    def test_enantiomer(self, enantiomer_graph1, enantiomer_graph2):
        assert enantiomer_graph1.enantiomer().is_isomorphic(enantiomer_graph2)
        assert not enantiomer_graph1.is_isomorphic(enantiomer_graph2)

    def test_hash_enantiomers(self, enantiomer_graph1, enantiomer_graph2):
        assert enantiomer_graph1._atom_stereo != enantiomer_graph2._atom_stereo
        assert hash(enantiomer_graph1.copy(frozen=True)) != hash(
            enantiomer_graph2.copy(frozen=True)
        )

    def test_inchi_coords(self):
        pytest.skip("RDKit construction tests moved to test_rdkit_conversion.py")
