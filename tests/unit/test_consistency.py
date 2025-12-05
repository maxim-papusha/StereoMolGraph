import rdkit.Chem  # type: ignore
import rdkit.Chem.rdDistGeom  # type: ignore
import pytest

from stereomolgraph import StereoMolGraph
from stereomolgraph.coords import Geometry
from stereomolgraph.rdmol2graph import RDMol2StereoMolGraph
from stereomolgraph.stereodescriptors import (
    TrigonalBipyramidal,
    Octahedral,
)

@pytest.fixture(scope="module")
def rdmol2graph():
    return RDMol2StereoMolGraph(stereo_complete = True,
        use_atom_map_number = False,
        lone_pair_stereo = False,
        resonance = True)

class TestRDKitConversion:


    def test_from_rdmol(self):
        rdmol = rdkit.Chem.MolFromSmiles(
            "[C:1]([O:2][C:33]([C:4]([O:5][H:13])"
            "([H:111])[H:12])([H:9])[H:10])([H:6])"
            "([H:77])[H:8]",
            sanitize=False,
        )
        rdmol = rdkit.Chem.AddHs(rdmol, explicitOnly=True)
        mol_graph = RDMol2StereoMolGraph(stereo_complete = True,
        use_atom_map_number = True,
        lone_pair_stereo = False,
        resonance = True)(rdmol)

        assert set(mol_graph.atoms) == set(
            (1, 2, 33, 4, 5, 6, 77, 8, 9, 10, 111, 12, 13)
        )

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

    @pytest.mark.parametrize(
        "inchi",
        [
            (r"InChI=1S/C3H8O/c1-3(2)4/h3-4H,1-2H3"),
            (
                r"InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3"
            ),
        ],
        ids=["isopropanol", "caffeine"],
    )
    def test_from_rdmol_to_rdmol_not_chiral(self, inchi, rdmol2graph):
        rdmol = rdkit.Chem.MolFromInchi(inchi, sanitize=False, removeHs=False)
        rdmol = rdkit.Chem.AddHs(rdmol, explicitOnly=True)
        assert inchi == rdkit.Chem.MolToInchi(rdmol, treatWarningAsError=True)  # type: ignore

        smg = rdmol2graph(rdmol)
        rdmol2, _ = smg._to_rdmol(generate_bond_orders=True)
        assert inchi == rdkit.Chem.MolToInchi(rdmol2, treatWarningAsError=True)  # type: ignore

    @pytest.mark.parametrize(
        "inchi",
        [
            (r"InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6+/m1/s1"),
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
    def test_from_rdmol_to_rdmol_stereo(self, inchi, rdmol2graph):
        rdmol = rdkit.Chem.MolFromInchi(inchi, sanitize=False)
        rdmol = rdkit.Chem.AddHs(rdmol, explicitOnly=True)
        smg = rdmol2graph(rdmol)
        rdmol2, _ = smg._to_rdmol(generate_bond_orders=True)
        molblock = rdkit.Chem.MolToMolBlock(rdmol2)
        assert inchi == rdkit.Chem.MolBlockToInchi(molblock)  # type: ignore

    @pytest.mark.parametrize(
        "inchi",
        [
            (r"InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6+/m1/s1"),
            (r"InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"),
            (r"InChI=1/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3"),
            (r"InChI=1/C6H9N3O2/c7-5(6(10)11)1-4-2-8-3-9-4/h2-3,5H,1,7H2,(H,8,9)(H,10,11)/t5-/m0/s1/f/h8,10H"),
            (r"InChI=1/C6H9N3O2/c7-5(6(10)11)1-4-2-8-3-9-4/h2-3,5H,1,7H2,(H,8,9)(H,10,11)/t5-/m0/s1/f/h9-10H"),
        ],
        ids=["alpha-D-gulopyranose", "benzene", "caffeine","L-Histidine N3-H", "L-Histidine N1-H"],
    )
    def test_from_rdmol_eq_from_geometry(self, inchi, rdmol2graph):
        rdmol = rdkit.Chem.MolFromInchi(inchi)
        rdmol = rdkit.Chem.AddHs(rdmol, addCoords=True, explicitOnly=True)
        graph_mol = rdmol2graph(rdmol)

        rdkit.Chem.rdDistGeom.EmbedMolecule(rdmol)
        xyz_str = rdkit.Chem.MolToXYZBlock(rdmol)
        geo = Geometry.from_xyz(xyz_str)
        graph_geo = StereoMolGraph.from_geometry(geo)

        assert graph_mol == graph_geo


    @pytest.mark.parametrize(
        "smiles",
        [
            "[H][Pt@SP1](F)(Cl)Br",
            "Cl[Pt@SP1](Cl)([NH3])[NH3]",
            "Cl[Pt@SP2](Cl)([NH3])[NH3]",
        ],
    )
    def test_from_rdmol_to_rdmol_square_planar(self, smiles, rdmol2graph):
        rdmol = rdkit.Chem.MolFromSmiles(smiles, sanitize=True)
        rdmol = rdkit.Chem.AddHs(rdmol, explicitOnly=True)
        smg = rdmol2graph(rdmol)
        rdmol2, _ = smg._to_rdmol(
            generate_bond_orders=True, allow_charged_fragments=True
        )
        rdkit.Chem.SanitizeMol(
            rdmol2, sanitizeOps=rdkit.Chem.SanitizeFlags.SANITIZE_ALL
        )
        rdkit.Chem.SanitizeMol(
            rdmol, sanitizeOps=rdkit.Chem.SanitizeFlags.SANITIZE_ALL
        )
        for atom in rdmol2.GetAtoms():
            atom.SetAtomMapNum(0)
        assert rdkit.Chem.MolToSmiles(rdmol) == rdkit.Chem.MolToSmiles(rdmol2)

    def test_from_rdmol_square_planar_different(self,rdmol2graph):
        smiles = ["Cl[Pt@SP1](Cl)([NH3])[NH3]", "Cl[Pt@SP2](Cl)([NH3])[NH3]"]
        rdmols = [rdkit.Chem.MolFromSmiles(s, sanitize=False) for s in smiles]
        rdmol = [rdkit.Chem.AddHs(mol, explicitOnly=True) for mol in rdmols]
        molgraphs = [rdmol2graph(mol) for mol in rdmol]
        assert molgraphs[0] != molgraphs[1]

    def test_from_rdmol_square_planar(self,rdmol2graph):
        smiles = (
            "C[Pt@SP1](F)(Cl)[H]",
            "C[Pt@SP2](Cl)(F)[H]",
            "C[Pt@SP3](F)([H])Cl",
        )
        mols = [rdkit.Chem.MolFromSmiles(i, sanitize=False) for i in smiles]
        mols = [rdkit.Chem.AddHs(i, explicitOnly=True) for i in mols]
        molgraphs = [rdmol2graph(i) for i in mols]
        assert all(molgraph == molgraphs[0] for molgraph in molgraphs)

    def test_from_rdmol_trigonal_bipyramidal(self, rdmol2graph):
        smiles = (
            "S[As@TB1](F)(Cl)(Br)N",
            "S[As@TB2](F)(Br)(Cl)N",
            "S[As@TB3](F)(Cl)(N)Br",
            "S[As@TB4](F)(Br)(N)Cl",
            "S[As@TB5](F)(N)(Cl)Br",
            "S[As@TB6](F)(N)(Br)Cl",
            "S[As@TB7](N)(F)(Cl)Br",
            "S[As@TB8](N)(F)(Br)Cl",
            "F[As@TB9](S)(Cl)(Br)N",
            "F[As@TB11](S)(Br)(Cl)N",
            "F[As@TB10](S)(Cl)(N)Br",
            "F[As@TB12](S)(Br)(N)Cl",
            "F[As@TB13](S)(N)(Cl)Br",
            "F[As@TB14](S)(N)(Br)Cl",
            "F[As@TB15](Cl)(S)(Br)N",
            "F[As@TB20](Br)(S)(Cl)N",
            "F[As@TB16](Cl)(S)(N)Br",
            "F[As@TB19](Br)(S)(N)Cl",
            "F[As@TB17](Cl)(Br)(S)N",
            "F[As@TB18](Br)(Cl)(S)N",
        )
        mols = [rdkit.Chem.MolFromSmiles(i, sanitize=False) for i in smiles]
        mols = [rdkit.Chem.AddHs(i, explicitOnly=True) for i in mols]
        molgraphs = [rdmol2graph(i) for i in mols]
        assert all(
            any(
                isinstance(stereo, TrigonalBipyramidal)
                for stereo in molgraph.stereo.values()
            )
            for molgraph in molgraphs
        )
        assert all(molgraph == molgraphs[0] for molgraph in molgraphs)

    def test_from_rdmol_octahedral(self, rdmol2graph):
        smiles = (
            "O[Co@OH1](Cl)(C)(N)(F)P",
            "O[Co@OH2](Cl)(F)(N)(C)P",
            "O[Co@OH3](Cl)(C)(N)(P)F",
            "O[Co@OH16](Cl)(F)(N)(P)C",
            "O[Co@OH6](Cl)(C)(P)(N)F",
            "O[Co@OH18](Cl)(F)(P)(N)C",
            "O[Co@OH19](Cl)(P)(C)(N)F",
            "O[Co@OH24](Cl)(P)(F)(N)C",
            "O[Co@OH25](P)(Cl)(C)(N)F",
            "O[Co@OH30](P)(Cl)(F)(N)C",
            "O[Co@OH4](Cl)(C)(F)(N)P",
            "O[Co@OH14](Cl)(F)(C)(N)P",
            "O[Co@OH5](Cl)(C)(F)(P)N",
            "O[Co@OH15](Cl)(F)(C)(P)N",
            "O[Co@OH7](Cl)(C)(P)(F)N",
            "O[Co@OH17](Cl)(F)(P)(C)N",
            "O[Co@OH20](Cl)(P)(C)(F)N",
            "O[Co@OH23](Cl)(P)(F)(C)N",
            "O[Co@OH26](P)(Cl)(C)(F)N",
            "O[Co@OH29](P)(Cl)(F)(C)N",
            "O[Co@OH10](Cl)(N)(F)(C)P",
            "O[Co@OH8](Cl)(N)(C)(F)P",
            "O[Co@OH11](Cl)(N)(F)(P)C",
            "O[Co@OH9](Cl)(N)(C)(P)F",
            "O[Co@OH13](Cl)(N)(P)(F)C",
            "O[Co@OH12](Cl)(N)(P)(C)F",
            "O[Co@OH22](Cl)(P)(N)(F)C",
            "O[Co@OH21](Cl)(P)(N)(C)F",
            "O[Co@OH28](P)(Cl)(N)(F)C",
            "O[Co@OH27](P)(Cl)(N)(C)F",
        )
        mols = [rdkit.Chem.MolFromSmiles(i, sanitize=False) for i in smiles]
        mols = [rdkit.Chem.AddHs(i, explicitOnly=True) for i in mols]
        molgraphs = [rdmol2graph(i) for i in mols]
        assert all(
            any(
                isinstance(stereo, Octahedral)
                for stereo in molgraph.stereo.values()
            )
            for molgraph in molgraphs
        )
        assert all(molgraph == molgraphs[0] for molgraph in molgraphs)

    def test_from_rdmol_octahedral_compare(self, rdmol2graph):
        identical1 = (
            "Cl[Co@OH1](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH2](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH3](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH4](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH5](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH14](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH15](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH16](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH21](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH22](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH27](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH28](N)(N)(O)(Cl)Cl",
        )
        identical2 = (
            "Cl[Co@OH6](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH7](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH17](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH18](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH19](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH20](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH23](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH24](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH25](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH26](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH29](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH30](N)(N)(O)(Cl)Cl",
        )
        identical3 = (
            "Cl[Co@OH8](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH9](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH10](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH11](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH12](N)(N)(O)(Cl)Cl",
            "Cl[Co@OH13](N)(N)(O)(Cl)Cl",
        )

        mols1 = [rdkit.Chem.MolFromSmiles(i, sanitize=False) for i in identical1]
        mols1 = [rdkit.Chem.AddHs(i, explicitOnly=True) for i in mols1]
        molgraphs1 = [rdmol2graph(i) for i in mols1]
        mols2 = [rdkit.Chem.MolFromSmiles(i, sanitize=False) for i in identical2]
        mols2 = [rdkit.Chem.AddHs(i, explicitOnly=True) for i in mols2]
        molgraphs2 = [rdmol2graph(i) for i in mols2]
        mols3 = [rdkit.Chem.MolFromSmiles(i, sanitize=False) for i in identical3]
        mols3 = [rdkit.Chem.AddHs(i, explicitOnly=True) for i in mols3]
        molgraphs3 = [rdmol2graph(i) for i in mols3]

        assert all(molgraph == molgraphs1[0] for molgraph in molgraphs1)
        assert all(molgraph == molgraphs2[0] for molgraph in molgraphs2)
        assert all(molgraph == molgraphs3[0] for molgraph in molgraphs3)
        assert molgraphs1[0] != molgraphs2[0]
        assert molgraphs1[0] != molgraphs3[0]
        assert molgraphs2[0] != molgraphs3[0]

        set1 = {
            rdkit.Chem.MolToSmiles(
                molgraph._to_rdmol()[0],
                canonical=True,
                ignoreAtomMapNumbers=True,
                isomericSmiles=True,
            )
            for molgraph in molgraphs1
        }
        set2 = {
            rdkit.Chem.MolToSmiles(
                molgraph._to_rdmol()[0],
                canonical=True,
                ignoreAtomMapNumbers=True,
                isomericSmiles=True,
            )
            for molgraph in molgraphs2
        }
        set3 = {
            rdkit.Chem.MolToSmiles(
                molgraph._to_rdmol()[0],
                canonical=True,
                ignoreAtomMapNumbers=True,
                isomericSmiles=True,
            )
            for molgraph in molgraphs3
        }

        assert set() == set1 & set2 == set1 & set3 == set2 & set3

    def test_inchi_coords(self, rdmol2graph):
        g_inchi = (
            "InChI=1S/C9H16O6/c1-9(2)14-7-5(12)6(4(11)3-10)13-8(7)15-9/h4-8,10-12H,3H2,1-2H3/t4-,5+,6-,7-,8-/m1/s1"
        )
        g_mol = rdkit.Chem.AddHs(
            rdkit.Chem.MolFromInchi(g_inchi), explicitOnly=True, addCoords=True
        )
        rdkit.Chem.rdDistGeom.EmbedMolecule(g_mol)
        g_xyz_str = rdkit.Chem.MolToXYZBlock(g_mol)
        g_graph = rdmol2graph(g_mol)

        g_geo = Geometry.from_xyz(g_xyz_str)

        g_graph2 = StereoMolGraph.from_geometry(g_geo)
        assert g_graph.is_isomorphic(g_graph2)
