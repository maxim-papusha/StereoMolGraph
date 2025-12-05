import pytest
import rdkit
from rdkit import Chem

from stereomolgraph.rdmol2graph import RDMol2StereoMolGraph
from stereomolgraph.algorithms.color_refine import (
    color_refine_hash_smg,
)

# This hash is supposed to be consistent between
# different python and library verisions

tcases = [
    {
        "name": "caffeine",
        "smiles": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
        "fixedh_inchi": "InChI=1/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
        "expected_hash": 8267820574251983887,
    },
    {
        "name": "ethanol",
        "smiles": "CCO",
        "fixedh_inchi": "InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3",
        "expected_hash": 8386824054974349231,
    },
    {
        "name": "imidazole",
        "smiles": "c1c[nH]cn1",
        "fixedh_inchi": "InChI=1/C3H4N2/c1-2-5-3-4-1/h1-3H,(H,4,5)/f/h4H",
        "expected_hash": 1256177459527738206,
    },
    {
        "name": "lactose",
        "smiles": "OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@H]1O",
        "fixedh_inchi": "InChI=1/C12H22O11/c13-1-3-5(15)6(16)9(19)12(22-3)23-10-4(2-14)21-11(20)8(18)7(10)17/h3-20H,1-2H2/t3-,4-,5+,6+,7-,8-,9-,10-,11-,12+/m1/s1",
        "expected_hash": 7822475553335959685,
    },
    {
        "name": "hydroxy radical",
        "smiles": "[OH]",
        "fixedh_inchi": "InChI=1/HO/h1H",
        "expected_hash": -3770721040748917487,
    },
    {
        "name": "singlet oxygen",
        "smiles": "[O][O]",
        "fixedh_inchi": "InChI=1/O2/c1-2",
        "expected_hash": 4397632973345371281,
    },
    {
        "name": "triplet oxygen",
        "smiles": "[O][O]",
        "fixedh_inchi": "InChI=1/O2/c1-2",
        "expected_hash": 4397632973345371281,
    },
    {
        "name": "methyl radical",
        "smiles": "[CH3]",
        "fixedh_inchi": "InChI=1/CH3/h1H3",
        "expected_hash": -4784546520906446105,
    },
    {
        "name": "hydroxide",
        "smiles": "[OH-]",
        "fixedh_inchi": "InChI=1/H2O/h1H2/p-1",
        "expected_hash": -3770721040748917487,
    },
    {
        "name": "ammonium",
        "smiles": "[NH4+]",
        "fixedh_inchi": "InChI=1/H3N/h1H3/p+1",
        "expected_hash": 3357968715173358365,
    },
    {
        "name": "ammonia",
        "smiles": "N",
        "fixedh_inchi": "InChI=1/H3N/h1H3",
        "expected_hash": 1182729512510341822,
    },
    {
        "name": "NaCl",
        "smiles": "[Cl-].[Na+]",
        "fixedh_inchi": "InChI=1/ClH.Na/h1H;/q;+1/p-1",
        "expected_hash": 1267298407956022791,
    },
    {
        "name": "OH3+",
        "smiles": "[H][O+]([H])[H]",
        "fixedh_inchi": "InChI=1S/H2O/h1H2/p+1",
        "expected_hash": -8711848525738659963,
    },
    {
        "name": "linear trinitrogen",
        "smiles": "[N]=[N+]=[N-]",
        "fixedh_inchi": "InChI=1/N3/c1-3-2",
        "expected_hash": 1748074141585181045,
    },
    {
        "name": "cyclic trinitrogen",
        "smiles": "[N]1N=N1",
        "fixedh_inchi": "InChI=1/N3/c1-2-3-1",
        "expected_hash": -733687595351923137,
    },
    {
        "name": "L-Histidine N3-H",
        "smiles": "O=C([C@H](CC1=CNC=N1)N)O",
        "fixedh_inchi": "InChI=1/C6H9N3O2/c7-5(6(10)11)1-4-2-8-3-9-4/h2-3,5H,1,7H2,(H,8,9)(H,10,11)/t5-/m0/s1/f/h8,10H",
        "expected_hash": 8220861725675873309,
    },
    {
        "name": "L-Histidine N1-H",
        "smiles": "O=C([C@H](CC1=CN=CN1)N)O",
        "fixedh_inchi": "InChI=1/C6H9N3O2/c7-5(6(10)11)1-4-2-8-3-9-4/h2-3,5H,1,7H2,(H,8,9)(H,10,11)/t5-/m0/s1/f/h9-10H",
        "expected_hash": 2704466096105967559,
    },
]


class TestHashConsistency:
    @pytest.fixture(scope="class")
    def rdmol2graph(self):
        return RDMol2StereoMolGraph(
            stereo_complete=True,
            resonance=True,
            lone_pair_stereo=False,
        )

    @pytest.mark.parametrize(
        "name,smiles,fixedh_inchi,expected_hash",
        [
            (c["name"], c["smiles"], c["fixedh_inchi"], c["expected_hash"])
            for c in tcases
        ],
    )
    def test_all(
        self,
        name,
        smiles,
        fixedh_inchi,
        expected_hash,
        rdmol2graph: RDMol2StereoMolGraph,
    ):
        rdmol1 = Chem.MolFromSmiles(smiles)
        rdmol1 = Chem.AddHs(rdmol1, explicitOnly=True)
        graph1 = rdmol2graph(rdmol1)
        rdmol2 = Chem.MolFromInchi(fixedh_inchi, sanitize=True)
        rdmol2 = Chem.AddHs(rdmol2, explicitOnly=True)
        graph2 = rdmol2graph(rdmol2)

        assert graph1 == graph2
        assert (
            color_refine_hash_smg(graph1)
            == color_refine_hash_smg(graph2)
            == expected_hash
        )
        
