import rdkit.Chem  # type: ignore

from stereomolgraph import StereoMolGraph
from stereomolgraph.experimental._sym_num_old import (
    external_symmetry_number,
    external_symmetry_number_bond,
    external_symmetry_number_per_rotatable_bond,
)


def _graph_from_inchi(inchi: str) -> StereoMolGraph:
    rdmol = rdkit.Chem.AddHs(rdkit.Chem.MolFromInchi(inchi))
    return StereoMolGraph.from_rdmol(rdmol, stereo_complete=True)


def test_external_symmetry_number_per_rotatable_bond() -> None:
    graph = _graph_from_inchi(
        "InChI=1S/C16H34/c1-3-5-7-9-11-13-15-16-14-12-10-8-6-4-2/h3-16H2,1-2H3"
    )

    per_bond_values = external_symmetry_number_per_rotatable_bond(graph)
    per_bond = {tuple(sorted(bond)): value for bond, value in per_bond_values.items()}

    assert external_symmetry_number(graph) == 2
    assert per_bond == {
        (0, 2): 3,
        (1, 3): 3,
    }

    for bond, value in per_bond_values.items():
        assert external_symmetry_number_bond(graph, bond) == value
