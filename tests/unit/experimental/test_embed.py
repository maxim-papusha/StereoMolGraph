"""Tests for the embed module: SMG → 3D coords → SMG round-trip."""

from __future__ import annotations

import pytest
import rdkit.Chem as Chem

from stereomolgraph import StereoMolGraph
from stereomolgraph.experimental._embed import EmbedParameters, SMG2Geo
from stereomolgraph.rdmol2graph import RDMol2StereoMolGraph

# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def rdmol2graph() -> RDMol2StereoMolGraph:
    """Converter from RDKit mol to StereoMolGraph with full stereo."""
    return RDMol2StereoMolGraph(
        stereo_complete=True,
        use_atom_map_number=False,
        lone_pair_stereo=False,
        resonance=True,
    )


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _smiles_to_smg(
    smiles: str,
    converter: RDMol2StereoMolGraph | None = None,
) -> StereoMolGraph:
    """Build a StereoMolGraph from a SMILES string with explicit hydrogens."""
    if converter is None:
        converter = RDMol2StereoMolGraph(
            stereo_complete=True,
            use_atom_map_number=False,
            lone_pair_stereo=False,
            resonance=True,
        )
    rdmol = Chem.MolFromSmiles(smiles)
    rdmol = Chem.AddHs(rdmol)
    return converter(rdmol)


# ---------------------------------------------------------------------------
# Parametrized round-trip tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "smiles",
    [
        # Single tetrahedral center
        "F[C@H](Cl)Br",
        "[C@H](F)(Cl)(Br)I",
        "N[C@@H](C)C(=O)O",
        # Two tetrahedral centers
        "C[C@H](O)[C@@H](C)O",
        # E/Z double bond stereo
        "F/C=C/Cl",
        "F/C=C\\Cl",
        "C/C=C\\C",
        # Ring with stereo
        "C1CC[C@@H](O)C1",
    ],
    ids=[
        "(S)-bromochlorofluoromethane",
        "(R)-bromochlorofluoroiodomethane",
        "D-alanine",
        "(2R,3R)-2,3-butanediol",
        "(E)-1,2-difluoroethylene",
        "(Z)-1,2-difluoroethylene",
        "(Z)-2-butene",
        "(R)-cyclopentanol",
    ],
)
def test_embed_roundtrip(
    smiles: str,
    rdmol2graph: RDMol2StereoMolGraph,
) -> None:
    """SMILES → SMG → embed → Geometry → SMG — both SMGs must be equal."""
    smg = _smiles_to_smg(smiles, rdmol2graph)

    # Embed to 3D coordinates (skip UFF optimization to avoid bond-order
    # issues with exotic atom types; ETKDG alone is sufficient for
    # stereo-determining bond lengths).
    embed = SMG2Geo(params=EmbedParameters(optimize=False))
    geos = list(embed(smg))
    assert len(geos) == 1, "Expected exactly one conformer"
    geo = geos[0]

    # Verify the geometry has the right number of atoms
    assert len(geo.atom_types) == smg.n_atoms

    # Rebuild SMG from the embedded geometry
    smg_from_geo = StereoMolGraph.from_geometry(geo)

    # Both SMGs should describe the same molecule
    assert smg == smg_from_geo, (
        f"Round-trip mismatch for {smiles}\n"
        f"Original  atom_stereo: {smg._atom_stereo}\n"
        f"Rebuilt   atom_stereo: {smg_from_geo._atom_stereo}\n"
        f"Original  bond_stereo: {smg._bond_stereo}\n"
        f"Rebuilt   bond_stereo: {smg_from_geo._bond_stereo}"
    )


# ---------------------------------------------------------------------------
# Multiple conformer test
# ---------------------------------------------------------------------------


def test_embed_multiple_conformers(rdmol2graph: RDMol2StereoMolGraph) -> None:
    """Each conformer from a multi-conformer embed must yield an equal SMG."""
    smiles = "CCO"  # ethanol — simple, flexible molecule
    smg = _smiles_to_smg(smiles, rdmol2graph)

    n_conformers = 3
    embed = SMG2Geo(params=EmbedParameters(optimize=False, n_conformers=n_conformers))
    geos = list(embed(smg))
    assert len(geos) == n_conformers

    for i, geo in enumerate(geos):
        smg_from_geo = StereoMolGraph.from_geometry(geo)
        assert smg == smg_from_geo, f"Conformer {i} mismatch for {smiles}"


# ---------------------------------------------------------------------------
# Fixed coordinate constraint test
# ---------------------------------------------------------------------------


def test_embed_with_fixed_coords(rdmol2graph: RDMol2StereoMolGraph) -> None:
    """Embedding with fixed atom coordinates must still produce an equal SMG."""
    smiles = "CCO"
    smg = _smiles_to_smg(smiles, rdmol2graph)

    # Fix the oxygen atom at the origin
    # SMG atom IDs go 0,1,2,...; we need the oxygen ID
    o_id = None
    for a in smg.atoms:
        if smg.get_atom_type(a) == 8:  # oxygen
            o_id = a
            break
    assert o_id is not None, "Ethanol must have an oxygen atom"

    embed = SMG2Geo(params=EmbedParameters(optimize=False))
    geos = list(embed(smg, fixed_coords={o_id: (0.0, 0.0, 0.0)}))
    assert len(geos) == 1
    geo = geos[0]

    # Check that the oxygen is near the origin
    import numpy as np

    o_coords = geo.coords[o_id]
    assert np.allclose(o_coords, [0.0, 0.0, 0.0], atol=1e-3), (
        f"Fixed atom {o_id} not at origin: {o_coords}"
    )

    smg_from_geo = StereoMolGraph.from_geometry(geo)
    assert smg == smg_from_geo


# ---------------------------------------------------------------------------
# Consistency: atom types and bond topology
# ---------------------------------------------------------------------------


def test_embed_preserves_atom_types(rdmol2graph: RDMol2StereoMolGraph) -> None:
    """The embedded geometry must have the same atom types as the SMG."""
    smiles = "C[C@H](O)[C@@H](C)O"
    smg = _smiles_to_smg(smiles, rdmol2graph)

    embed = SMG2Geo(params=EmbedParameters(optimize=False))
    geos = list(embed(smg))
    geo = geos[0]

    assert geo.atom_types == smg.atom_types, (
        f"Atom type mismatch: {geo.atom_types} vs {smg.atom_types}"
    )
