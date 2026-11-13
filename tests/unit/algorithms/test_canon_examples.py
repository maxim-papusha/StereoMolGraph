"""Canonization of the examples already used by the four graph test suites."""

from random import Random

import pytest

from stereomolgraph import (
    CondensedReactionGraph,
    MolGraph,
    StereoCondensedReactionGraph,
    StereoMolGraph,
)
from stereomolgraph.coords import Geometry
from stereomolgraph.algorithms.canon import canon_atom_num
from stereomolgraph.graphs.crg import Change
from tests.unit.graphs import test_crg as crg_examples
from tests.unit.graphs import test_mg as mg_examples
from tests.unit.graphs import test_scrg as scrg_examples

from test_canon import _ORDER_CHANGES, _permuted


def _canonical_result(graph, mapping):
    """Serialize the supported graph data independently of search certificates.

    Comparing graph equality would hide unspecified parity and descriptor-class
    differences. Symmetric atoms may exchange numbers, so compare numbered data
    rather than atom-to-number dictionaries.
    """

    def site_key(site):
        if isinstance(site, frozenset):
            return "bond", tuple(sorted(mapping[atom] for atom in site))
        return "atom", (mapping[site],)

    def descriptor_key(stereo):
        descriptor = stereo.__class__(
            tuple(0 if atom is None else mapping[atom] for atom in stereo.atoms),
            stereo.parity,
        )
        return (
            stereo.__class__.__module__,
            stereo.__class__.__qualname__,
            descriptor.canonical_form(),
        )

    atoms = tuple(sorted(zip((mapping[a] for a in graph.atoms), graph.atom_types)))
    bonds = tuple(
        sorted(
            (
                tuple(sorted(mapping[a] for a in bond)),
                attributes["reaction"].value if attributes.get("reaction") else "",
            )
            for bond, attributes in graph.bonds_with_attributes.items()
        )
    )
    stereo = (
        tuple(
            sorted(
                (site_key(site), descriptor_key(s)) for site, s in graph.stereo.items()
            )
        )
        if isinstance(graph, StereoMolGraph)
        else ()
    )
    changes = (
        tuple(
            sorted(
                (site_key(site), change.value, descriptor_key(descriptor))
                for records in (graph.atom_stereo_changes, graph.bond_stereo_changes)
                for site, phases in records.items()
                for change, descriptor in phases.items()
                if descriptor is not None
            )
        )
        if isinstance(graph, StereoCondensedReactionGraph)
        else ()
    )
    return atoms, bonds, stereo, changes


def _assert_enumeration_requirements(graph, mapping):
    """Check public ordering rules without using the search's candidate labels."""
    assert set(mapping) == set(graph.atoms)
    assert set(mapping.values()) == set(range(1, len(graph) + 1))
    types = dict(zip(graph.atoms, graph.atom_types))
    changed_bonds = {
        bond
        for bond, attrs in graph.bonds_with_attributes.items()
        if attrs.get("reaction")
    }
    active = {
        atom
        for bond in changed_bonds
        if graph.bonds_with_attributes[bond]["reaction"]
        in (Change.FORMED, Change.BROKEN)
        for atom in bond
    }
    stereo_atoms = (
        {
            atom
            for records in (graph.atom_stereo_changes, graph.bond_stereo_changes)
            for phases in records.values()
            for stereo in phases.values()
            if stereo is not None
            for atom in stereo.atoms
            if atom is not None
        }
        if isinstance(graph, StereoCondensedReactionGraph)
        else set()
    ) - active
    assert {mapping[a] for a in active} == set(range(1, len(active) + 1))
    assert {mapping[a] for a in stereo_atoms} == set(
        range(len(active) + 1, len(active) + len(stereo_atoms) + 1)
    )

    assigned = set()
    remaining = set(graph.atoms)
    for atom in sorted(mapping, key=mapping.get):
        frontier = {
            nbr for previous in assigned for nbr in graph.neighbors[previous]
        } & remaining
        candidates = frontier or remaining
        for group in (active, stereo_atoms):
            group_remaining = group & remaining
            if group_remaining:
                candidates = candidates & group_remaining or group_remaining
                break
        assert atom in candidates  # Finish the frontier before starting elsewhere.
        heavy = {a for a in candidates if types[a] != 1}
        if heavy:
            assert atom in heavy
            candidates = heavy
        changed_connections = {
            a
            for a in candidates
            if any(frozenset((a, b)) in changed_bonds for b in assigned)
        }
        if changed_connections:
            assert atom in changed_connections
        assigned.add(atom)
        remaining.remove(atom)


def _assert_invariant(graph, changes):
    reference_mapping = canon_atom_num(graph)
    _assert_enumeration_requirements(graph, reference_mapping)
    expected = _canonical_result(graph, reference_mapping)
    rng = Random(2026)
    for _ in range(3):
        permuted, relabeling = _permuted(graph, rng, changes)
        # First verify the test transformation preserves every supported feature.
        corresponding_mapping = {
            relabeling[atom]: number for atom, number in reference_mapping.items()
        }
        assert _canonical_result(permuted, corresponding_mapping) == expected
        for pruning in (False, True):
            mapping = canon_atom_num(permuted, prune_automorphisms=pruning)
            _assert_enumeration_requirements(permuted, mapping)
            assert _canonical_result(permuted, mapping) == expected


_PERMUTATIONS = [pytest.param((change,), id=change) for change in _ORDER_CHANGES] + [
    pytest.param(_ORDER_CHANGES, id="all-orders")
]


class _GraphExamples:
    # Reuse just the fixtures, without inheriting and collecting graph API tests.
    empty_mol_graph = mg_examples.TestMolGraph.empty_mol_graph
    mol_graph = mg_examples.TestMolGraph.mol_graph
    water_graph = mg_examples.TestMolGraph.water_graph
    enantiomer_graph1 = mg_examples.TestMolGraph.enantiomer_graph1
    enantiomer_graph2 = mg_examples.TestMolGraph.enantiomer_graph2
    chiral_reactant_geo = mg_examples.TestMolGraph.chiral_reactant_geo
    chiral_reactant_graph = mg_examples.TestMolGraph.chiral_reactant_graph
    chiral_product_geo1 = mg_examples.TestMolGraph.chiral_product_geo1
    chiral_product_geo2 = mg_examples.TestMolGraph.chiral_product_geo2
    chiral_product_graph1 = mg_examples.TestMolGraph.chiral_product_graph1
    chiral_product_graph2 = mg_examples.TestMolGraph.chiral_product_graph2


class TestMolecularExamples(_GraphExamples):
    @pytest.fixture(
        params=[
            MolGraph,
            StereoMolGraph,
            CondensedReactionGraph,
            StereoCondensedReactionGraph,
        ],
        ids=lambda graph_type: graph_type.__name__,
        autouse=True,
    )
    def graph_class(self, request):
        self._TestClass = request.param

    @pytest.mark.parametrize(
        "example",
        [
            "empty_mol_graph",
            "mol_graph",
            "water_graph",
            "enantiomer_graph1",
            "enantiomer_graph2",
            "chiral_reactant_graph",
            "chiral_product_graph1",
            "chiral_product_graph2",
        ],
    )
    @pytest.mark.parametrize("changes", _PERMUTATIONS)
    def test_enumeration_and_order_invariance(self, request, example, changes):
        _assert_invariant(request.getfixturevalue(example), changes)

    def test_enantiomers_and_products(
        self,
        enantiomer_graph1,
        enantiomer_graph2,
        chiral_product_graph1,
        chiral_product_graph2,
    ):
        def result(graph):
            return _canonical_result(graph, canon_atom_num(graph))

        assert result(chiral_product_graph1) == result(chiral_product_graph2)
        same_enantiomer_result = result(enantiomer_graph1) == result(enantiomer_graph2)
        assert same_enantiomer_result == (
            not isinstance(enantiomer_graph1, StereoMolGraph)
        )


class TestReactionExamples(_GraphExamples):
    crg = crg_examples.TestCondensedReactionGraph.crg
    chiral_reaction_scrg1 = (
        crg_examples.TestCondensedReactionGraph.chiral_reaction_scrg1
    )
    chiral_reaction_scrg2 = (
        crg_examples.TestCondensedReactionGraph.chiral_reaction_scrg2
    )

    @pytest.fixture(
        params=[CondensedReactionGraph, StereoCondensedReactionGraph],
        ids=lambda graph_type: graph_type.__name__,
        autouse=True,
    )
    def graph_class(self, request):
        self._TestClass = request.param

    @pytest.mark.parametrize(
        "example", ["crg", "chiral_reaction_scrg1", "chiral_reaction_scrg2"]
    )
    @pytest.mark.parametrize("changes", _PERMUTATIONS)
    def test_enumeration_and_order_invariance(self, request, example, changes):
        _assert_invariant(request.getfixturevalue(example), changes)

    def test_reaction_stereo_distinction(
        self, chiral_reaction_scrg1, chiral_reaction_scrg2
    ):
        first = _canonical_result(
            chiral_reaction_scrg1, canon_atom_num(chiral_reaction_scrg1)
        )
        second = _canonical_result(
            chiral_reaction_scrg2, canon_atom_num(chiral_reaction_scrg2)
        )
        assert (first == second) == (self._TestClass is CondensedReactionGraph)


class TestStereoReactionExamples(_GraphExamples):
    _TestClass = StereoCondensedReactionGraph
    chiral_ts_geo1 = scrg_examples.TestStereoCondensedReactionGraph.chiral_ts_geo1
    chiral_ts_geo2 = scrg_examples.TestStereoCondensedReactionGraph.chiral_ts_geo2
    scrg_stereo_inversion = (
        scrg_examples.TestStereoCondensedReactionGraph.scrg_stereo_inversion
    )
    chiral_reaction_chiral_ts_scrg1 = (
        scrg_examples.TestStereoCondensedReactionGraph.chiral_reaction_chiral_ts_scrg1
    )
    chiral_reaction_chiral_ts_scrg2 = (
        scrg_examples.TestStereoCondensedReactionGraph.chiral_reaction_chiral_ts_scrg2
    )

    @pytest.fixture
    def methylamine_reaction(self, data_path):
        # The chain-of-states example in test_scrg.py includes ordinary stereo,
        # multiple atom stereo changes, a bond stereo change and changed bonds.
        states = [
            Geometry.from_xyz_file(
                data_path / f"methylamine_phosgenation_trans_{state}.xyz"
            )
            for state in ("r", "p", "ts")
        ]
        return self._TestClass.from_geometries(*states)

    @pytest.mark.parametrize(
        "example",
        [
            "scrg_stereo_inversion",
            "methylamine_reaction",
            "chiral_reaction_chiral_ts_scrg1",
            "chiral_reaction_chiral_ts_scrg2",
        ],
    )
    @pytest.mark.parametrize("changes", _PERMUTATIONS)
    def test_enumeration_and_order_invariance(self, request, example, changes):
        _assert_invariant(request.getfixturevalue(example), changes)

    def test_transition_states_remain_distinct(
        self, chiral_reaction_chiral_ts_scrg1, chiral_reaction_chiral_ts_scrg2
    ):
        first = _canonical_result(
            chiral_reaction_chiral_ts_scrg1,
            canon_atom_num(chiral_reaction_chiral_ts_scrg1),
        )
        second = _canonical_result(
            chiral_reaction_chiral_ts_scrg2,
            canon_atom_num(chiral_reaction_chiral_ts_scrg2),
        )
        assert first != second
