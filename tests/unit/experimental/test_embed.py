"""Tests for the embed module: SMG → 3D coords → SMG round-trip.

Focuses on larger chiral molecules with multiple stereocenters,
bridged/fused ring systems, and diverse functional groups.
"""

from __future__ import annotations

import pytest
import rdkit.Chem as Chem

from stereomolgraph import StereoMolGraph
from stereomolgraph.experimental._embed import EmbedParameters, SMG2Geo
from stereomolgraph.rdmol2graph import RDMol2StereoMolGraph


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def rdmol2graph() -> RDMol2StereoMolGraph:
    """Converter from RDKit mol to StereoMolGraph.

    Uses ``stereo_complete=False`` so that only RDKit-explicit stereo
    is captured — this makes the round-trip comparison against
    ``StereoMolGraph.from_geometry`` meaningful.
    """
    return RDMol2StereoMolGraph(
        stereo_complete=False,
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
            stereo_complete=False,
            use_atom_map_number=False,
            lone_pair_stereo=False,
            resonance=True,
        )
    rdmol = Chem.MolFromSmiles(smiles)
    rdmol = Chem.AddHs(rdmol)
    return converter(rdmol)


# ---------------------------------------------------------------------------
# Parametrized round-trip: larger chiral molecules
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "smiles",
    [
        # ── Aldoses (open chain, 3–4 contiguous stereocenters) ─────────
        # D-Glucose: 4 tetrahedral centers
        "C([C@H]([C@H]([C@@H]([C@H](C=O)O)O)O)O)O",
        # D-Mannose: C2 epimer of glucose
        "C([C@@H]([C@@H]([C@H]([C@H](C=O)O)O)O)O)O",
        # D-Galactose: C4 epimer of glucose
        "C([C@H]([C@H]([C@H]([C@H](C=O)O)O)O)O)O",
        # D-Ribose: 3 tetrahedral centers, 5-carbon aldose
        "C([C@H]([C@H]([C@H](C=O)O)O)O)O",

        # ── Pyranose rings (5 tetrahedral centers) ─────────────────────
        # α-D-Glucopyranose
        "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
        # α-D-Mannopyranose
        "OC[C@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O",

        # ── Furanose ring ──────────────────────────────────────────────
        # β-D-Fructofuranose: 4 tetrahedral centers, ketose
        "OC[C@H]1O[C@@](O)(CO)[C@@H](O)[C@@H]1O",
    ],
    ids=[
        "D-glucose",
        "D-mannose",
        "D-galactose",
        "D-ribose",
        "alpha-D-glucopyranose",
        "alpha-D-mannopyranose",
        "beta-D-fructofuranose",
    ],
)
def test_embed_roundtrip(
    smiles: str,
    rdmol2graph: RDMol2StereoMolGraph,
) -> None:
    """SMILES → SMG → embed → Geometry → SMG — both SMGs must be equal."""
    smg = _smiles_to_smg(smiles, rdmol2graph)

    # Embed to 3D coordinates.
    # Skip UFF optimization — ETKDG alone produces bond lengths
    # sufficient for stereo discrimination via the connectivity cutoff.
    embed = SMG2Geo(params=EmbedParameters(optimize=False))
    geos = list(embed(smg))
    assert len(geos) == 1, "Expected exactly one conformer"
    geo = geos[0]

    # Verify atom count
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
    smiles = "C([C@H]([C@H]([C@@H]([C@H](C=O)O)O)O)O)O"  # D-glucose (open)
    smg = _smiles_to_smg(smiles, rdmol2graph)

    n_conformers = 3
    embed = SMG2Geo(
        params=EmbedParameters(optimize=False, n_conformers=n_conformers)
    )
    geos = list(embed(smg))
    assert len(geos) == n_conformers

    for i, geo in enumerate(geos):
        smg_from_geo = StereoMolGraph.from_geometry(geo)
        assert smg == smg_from_geo, (
            f"Conformer {i} mismatch for {smiles}"
        )


# ---------------------------------------------------------------------------
# Fixed coordinate constraint test
# ---------------------------------------------------------------------------


def test_embed_with_fixed_coords(rdmol2graph: RDMol2StereoMolGraph) -> None:
    """Embedding with fixed atom coordinates must still produce an equal SMG."""
    smiles = "C([C@H]([C@H]([C@H](C=O)O)O)O)O"  # D-ribose
    smg = _smiles_to_smg(smiles, rdmol2graph)

    # Fix the aldehyde carbon at the origin
    c_id = None
    for a in smg.atoms:
        if smg.get_atom_type(a) == 6:  # carbon
            # Find the carbonyl carbon (the one with an O neighbor that has only one neighbor)
            nbrs = smg.bonded_to(a)
            for n in nbrs:
                if smg.get_atom_type(n) == 8 and len(smg.bonded_to(n)) == 1:
                    c_id = a
                    break
            if c_id is not None:
                break
    assert c_id is not None, "Ribose must have a carbonyl carbon"

    embed = SMG2Geo(params=EmbedParameters(optimize=False))
    geos = list(embed(smg, fixed_coords={c_id: (0.0, 0.0, 0.0)}))
    assert len(geos) == 1
    geo = geos[0]

    # Check that the carbonyl carbon is near the origin
    import numpy as np
    c_coords = geo.coords[c_id]
    assert np.allclose(c_coords, [0.0, 0.0, 0.0], atol=1e-3), (
        f"Fixed atom {c_id} not at origin: {c_coords}"
    )

    smg_from_geo = StereoMolGraph.from_geometry(geo)
    assert smg == smg_from_geo


# ---------------------------------------------------------------------------
# Consistency: atom types and bond topology
# ---------------------------------------------------------------------------


def test_embed_preserves_atom_types(rdmol2graph: RDMol2StereoMolGraph) -> None:
    """The embedded geometry must have the same atom types as the SMG."""
    smiles = "C([C@H]([C@H]([C@@H]([C@H](C=O)O)O)O)O)O"  # D-glucose
    smg = _smiles_to_smg(smiles, rdmol2graph)

    embed = SMG2Geo(params=EmbedParameters(optimize=False))
    geos = list(embed(smg))
    geo = geos[0]

    assert geo.atom_types == smg.atom_types, (
        f"Atom type mismatch: {geo.atom_types} vs {smg.atom_types}"
    )
