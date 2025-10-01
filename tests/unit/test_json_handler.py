import json
from pathlib import Path

import pytest

from stereomolgraph import (
    CondensedReactionGraph,
    MolGraph,
    StereoCondensedReactionGraph,
    StereoMolGraph,
)
from stereomolgraph.coords import Geometry
from stereomolgraph.experimental import JSONHandler
from stereomolgraph.graphs.crg import Change
from stereomolgraph.stereodescriptors import Tetrahedral


def _build_sample_molgraph() -> MolGraph:
    mol = MolGraph()
    for idx, atom_type in enumerate("CHHHCHCHH"):
        mol.add_atom(idx, atom_type)
    for a1, a2 in (
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        (4, 5),
        (4, 6),
        (5, 7),
        (5, 8),
    ):
        mol.add_bond(a1, a2)
    return mol


def _build_sample_crg() -> CondensedReactionGraph:
    crg = CondensedReactionGraph()
    for idx, atom_type in enumerate("HNOH"):
        crg.add_atom(idx, atom_type)
    crg.add_bond(0, 1, reaction=Change.FORMED)
    crg.add_bond(1, 2)
    crg.add_bond(2, 3, reaction=Change.BROKEN)
    crg.add_bond(1, 3)
    return crg


@pytest.fixture
def sample_molgraph() -> MolGraph:
    return _build_sample_molgraph()


@pytest.fixture
def sample_stereo_molgraph(sample_molgraph: MolGraph) -> StereoMolGraph:
    smg = StereoMolGraph(sample_molgraph)
    smg.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), -1))
    return smg


@pytest.fixture
def sample_crg() -> CondensedReactionGraph:
    return _build_sample_crg()


@pytest.fixture
def sample_scrg(data_path: Path) -> StereoCondensedReactionGraph:
    reactant = Geometry.from_xyz_file(
        data_path / "methylamine_phosgenation_trans_r.xyz"
    )
    product = Geometry.from_xyz_file(
        data_path / "methylamine_phosgenation_trans_p.xyz"
    )
    transition_state = Geometry.from_xyz_file(
        data_path / "methylamine_phosgenation_trans_ts.xyz"
    )
    return StereoCondensedReactionGraph.from_geometries(
        reactant,
        product,
        transition_state,
    )


def test_json_roundtrip_molgraph(sample_molgraph: MolGraph) -> None:
    serialized = JSONHandler.json_serialize(sample_molgraph)
    deserialized = JSONHandler.json_deserialize(serialized)

    assert isinstance(deserialized, MolGraph)
    assert deserialized == sample_molgraph


def test_json_roundtrip_stereo_molgraph(
    sample_stereo_molgraph: StereoMolGraph,
) -> None:
    serialized = JSONHandler.json_serialize(sample_stereo_molgraph)
    deserialized = JSONHandler.json_deserialize(serialized)

    assert isinstance(deserialized, StereoMolGraph)
    assert deserialized == sample_stereo_molgraph


def test_json_roundtrip_condensed_reaction_graph(
    sample_crg: CondensedReactionGraph,
) -> None:
    serialized = JSONHandler.json_serialize(sample_crg)
    deserialized = JSONHandler.json_deserialize(serialized)

    assert isinstance(deserialized, CondensedReactionGraph)
    assert deserialized == sample_crg
    assert deserialized.get_formed_bonds() == sample_crg.get_formed_bonds()
    assert deserialized.get_broken_bonds() == sample_crg.get_broken_bonds()


def test_json_roundtrip_stereo_condensed_reaction_graph(
    sample_scrg: StereoCondensedReactionGraph,
) -> None:
    serialized = JSONHandler.json_serialize(sample_scrg)
    deserialized = JSONHandler.json_deserialize(serialized)

    assert isinstance(deserialized, StereoCondensedReactionGraph)
    assert deserialized == sample_scrg

    assert deserialized.reactant() == sample_scrg.reactant()
    assert deserialized.product() == sample_scrg.product()

    for atom in sample_scrg.atom_stereo_changes:
        assert (
            deserialized.atom_stereo_changes[atom].keys()
            == sample_scrg.atom_stereo_changes[atom].keys()
        )

    for bond in sample_scrg.bond_stereo_changes:
        assert (
            deserialized.bond_stereo_changes[bond].keys()
            == sample_scrg.bond_stereo_changes[bond].keys()
        )


def test_json_roundtrip_from_dict(sample_molgraph: MolGraph) -> None:
    as_dict = JSONHandler.as_dict(sample_molgraph)
    serialized = json.dumps(as_dict)
    deserialized = JSONHandler.json_deserialize(serialized)

    assert isinstance(deserialized, MolGraph)
    assert deserialized == sample_molgraph
