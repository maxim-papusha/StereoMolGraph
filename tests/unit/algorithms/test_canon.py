"""Correctness of certificate and automorphism pruning, independent of atom IDs."""

from itertools import combinations
from random import Random

import pytest
from rdkit import Chem

from stereomolgraph import (
    CondensedReactionGraph,
    MolGraph,
    StereoCondensedReactionGraph,
    StereoMolGraph,
)
from stereomolgraph.algorithms import canon
from stereomolgraph.graphs.crg import Change
from stereomolgraph.algorithms.canon import (
    _candidate_label,
    _CanonState,
    _graph_certificate,
    _prune_orbits,
    _search,
    canon_atom_num,
    _initialize,
)
from stereomolgraph.stereodescriptors import PlanarBond, SquarePlanar, Tetrahedral


def _graph(n, edges):
    graph = MolGraph()
    for atom in range(n):
        graph.add_atom(atom, "C")
    for a, b in edges:
        graph.add_bond(a, b)
    return graph


def _score(params, mapping):
    partial = {}
    labels = []
    for atom in sorted(mapping, key=mapping.get):
        labels.append(_candidate_label(atom, params, partial))
        partial[atom] = mapping[atom]
    return labels + [()], _graph_certificate(params, mapping)


def _exhaustive_score(params, mapping=None, labels=()):
    """Small recursive oracle: no prefix pruning, symmetry, or frontier undo."""
    mapping = {} if mapping is None else mapping
    remaining = params.nbrhds.keys() - mapping.keys()
    if not remaining:
        return list(labels) + [()], _graph_certificate(params, mapping)
    frontier = {
        nbr for atom in mapping for nbr in params.nbrhds[atom] if nbr in remaining
    } or remaining
    for group in (params.active_atoms, params.stereo_change_atoms):
        group_remaining = remaining & group
        if group_remaining:
            frontier = frontier & group_remaining or group_remaining
            break
    candidates = {atom: _candidate_label(atom, params, mapping) for atom in frontier}
    best_label = min(candidates.values())
    return min(
        _exhaustive_score(
            params,
            {**mapping, atom: len(mapping) + 1},
            (*labels, best_label),
        )
        for atom, label in candidates.items()
        if label == best_label
    )


_ORDER_CHANGES = (
    "ids",
    "atoms",
    "bonds",
    "stereo",
    "stereo_changes",
    "descriptor_atoms",
)


def _reordered(items, rng):
    """Randomize order, forcing a nonidentity permutation when possible."""
    indices = list(range(len(items)))
    rng.shuffle(indices)
    if len(indices) > 1 and indices == list(range(len(items))):
        indices = indices[1:] + indices[:1]
    return [items[index] for index in indices]


def _permuted(graph, rng, changes=_ORDER_CHANGES):
    """Rebuild a graph with independently chosen ID and insertion-order changes."""
    atoms = list(graph.atoms)
    numbers = (
        rng.sample(range(-10000, 10000), len(atoms)) if "ids" in changes else atoms
    )
    mapping = dict(zip(atoms, numbers))
    if "atoms" in changes:
        atoms = _reordered(atoms, rng)
    result = graph.__class__()
    for atom in atoms:
        result.add_atom(mapping[atom], **graph.atoms_with_attributes[atom])
    bonds = list(graph.bonds)
    if "bonds" in changes:
        bonds = _reordered(bonds, rng)
    for bond in bonds:
        a, b = bond
        if "bonds" in changes and rng.randrange(2):
            a, b = b, a
        result.add_bond(mapping[a], mapping[b], **graph.bonds_with_attributes[bond])

    def remap_stereo(stereo):
        atoms = stereo.atoms
        if "descriptor_atoms" in changes:
            equivalents = [
                tuple(atoms[index] for index in permutation)
                for permutation in stereo.PERMUTATION_GROUP
                if tuple(atoms[index] for index in permutation) != atoms
            ]
            if equivalents:
                atoms = rng.choice(equivalents)
        return stereo.__class__(
            tuple(None if atom is None else mapping[atom] for atom in atoms),
            stereo.parity,
        )

    if isinstance(graph, StereoMolGraph):
        for records, setter in (
            (graph.atom_stereo, result.set_atom_stereo),
            (graph.bond_stereo, result.set_bond_stereo),
        ):
            descriptors = list(records.values())
            if "stereo" in changes:
                descriptors = _reordered(descriptors, rng)
            for stereo in descriptors:
                setter(remap_stereo(stereo))
    if isinstance(graph, StereoCondensedReactionGraph):
        for records, setter, target in (
            (
                graph.atom_stereo_changes,
                result.set_atom_stereo_change,
                result.atom_stereo_changes,
            ),
            (
                graph.bond_stereo_changes,
                result.set_bond_stereo_change,
                result.bond_stereo_changes,
            ),
        ):
            sites = list(records.items())
            if "stereo_changes" in changes:
                sites = _reordered(sites, rng)
            for site, change_records in sites:
                setter(
                    **{
                        change.value: remap_stereo(stereo)
                        for change, stereo in change_records.items()
                        if stereo is not None
                    }
                )
                if "stereo_changes" in changes:
                    mapped_site = (
                        frozenset(mapping[atom] for atom in site)
                        if isinstance(site, frozenset)
                        else mapping[site]
                    )
                    # Setters normalize keyword order; reorder the stored change
                    # dictionary too, so phase-order invariance is actually tested.
                    stored = target[mapped_site]
                    items = _reordered(list(stored.items()), rng)
                    stored.clear()
                    stored.update(items)
    return result, mapping


def _shuffled(graph, rng):
    return _permuted(graph, rng)[0]


def test_all_small_graphs_match_exhaustive_search():
    for n in range(5):
        edges = list(combinations(range(n), 2))
        for mask in range(1 << len(edges)):
            graph = _graph(n, [e for i, e in enumerate(edges) if mask & (1 << i)])
            params = _initialize(graph)
            expected = _exhaustive_score(params)
            for pruning in (False, True):
                mapping = canon_atom_num(graph, prune_automorphisms=pruning)
                assert _score(params, mapping) == expected, (n, mask, pruning)


def test_random_typed_graphs_match_exhaustive_search_and_relabeling():
    rng = Random(2026)
    for _ in range(40):
        n = rng.randrange(1, 8)
        graph = MolGraph()
        for atom in range(n):
            graph.add_atom(atom, rng.choice(("C", "N", "O", "H")))
        for a, b in combinations(range(n), 2):
            if rng.random() < 0.35:
                graph.add_bond(a, b)
        params = _initialize(graph)
        assert _score(params, canon_atom_num(graph)) == _exhaustive_score(params)
        shuffled = _shuffled(graph, rng)
        assert _score(_initialize(shuffled), canon_atom_num(shuffled)) == _score(
            params, canon_atom_num(graph)
        )


@pytest.mark.parametrize("pruning", [False, True])
def test_ascending_numbers_and_nonhydrogen_frontier_priority(pruning):
    graph = MolGraph()
    for atom, element in ((10, "C"), (20, "O"), (30, "H"), (40, "H")):
        graph.add_atom(atom, element)
    for atom in (20, 30, 40):
        graph.add_bond(10, atom)
    mapping = canon_atom_num(graph, prune_automorphisms=pruning)
    # The unique high-degree carbon starts; oxygen outranks frontier hydrogens.
    assert mapping[10] == 1
    assert mapping[20] == 2
    assert {mapping[30], mapping[40]} == {3, 4}
    assert list(mapping.values()) == [1, 2, 3, 4]


@pytest.mark.parametrize("pruning", [False, True])
def test_component_hydrogens_precede_next_component(pruning):
    graph = MolGraph()
    for atom, element in ((0, "C"), (1, "H"), (2, "O"), (3, "H")):
        graph.add_atom(atom, element)
    graph.add_bond(0, 1)
    graph.add_bond(2, 3)
    rng = Random(123)
    for _ in range(5):
        mapping = canon_atom_num(graph, prune_automorphisms=pruning)
        types = dict(zip(graph.atoms, graph.atom_types))
        order = sorted(graph.atoms, key=mapping.get)
        assert [types[atom] == 1 for atom in order] == [False, True, False, True]
        assert order[1] in graph.neighbors[order[0]]
        assert order[3] in graph.neighbors[order[2]]
        graph = _shuffled(graph, rng)


def test_nonhydrogen_priority_precedes_color_rarity():
    graph = MolGraph()
    for atom, element in ((0, "C"), (1, "C"), (2, "H")):
        graph.add_atom(atom, element)
    params = _initialize(graph)
    params.atom_color_frequency = {0: 2, 1: 2, 2: 1}
    assert _candidate_label(0, params, {}) < _candidate_label(2, params, {})


@pytest.mark.parametrize("kind", ["star", "cycle", "clique", "isolated", "components"])
def test_symmetric_graphs_prune_equivalent_siblings(kind, monkeypatch):
    n = 6
    edges = {
        "star": [(0, a) for a in range(1, n)],
        "cycle": [(a, (a + 1) % n) for a in range(n)],
        "clique": list(combinations(range(n), 2)),
        "isolated": [],
        "components": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)],
    }[kind]
    graph = _graph(n, edges)
    params = _initialize(graph)
    leaves = []
    record_leaf = canon._record_leaf

    def observe_leaf(params, state):
        leaves.append(state.current_mapping.copy())
        record_leaf(params, state)

    monkeypatch.setattr(canon, "_record_leaf", observe_leaf)
    pruned = _search(params)
    pruned_leaves = len(leaves)
    leaves.clear()
    baseline = _search(params, prune_automorphisms=False)
    assert _score(params, pruned.best_mapping) == _score(params, baseline.best_mapping)
    assert pruned_leaves < len(leaves)
    assert pruned.automorphisms


def test_orbit_pruning_requires_pointwise_prefix_stabilizer():
    state = _CanonState(
        current_mapping={0: 4, 1: 3},
        queue=[{3}],
        explored=[{2}],
        automorphisms=[{0: 1, 1: 0, 2: 3, 3: 2}],
    )
    _prune_orbits(state)
    assert state.queue[-1] == {3}  # Setwise fixation is insufficient.
    state.automorphisms.append({0: 0, 1: 1, 2: 3, 3: 2})
    _prune_orbits(state)
    assert state.queue[-1] == set()


def test_orbit_pruning_uses_generator_compositions():
    state = _CanonState(
        queue=[{1, 2}],
        explored=[{0}],
        automorphisms=[{0: 1, 1: 2, 2: 0}],
    )
    _prune_orbits(state)
    assert not state.queue[-1]


@pytest.mark.parametrize(
    "smiles",
    ["C", "CC", "C1CCCCC1", "F[C@](Cl)(Br)I", "F/C=C/F", "F/C=C\\F", "N[C@@](C)(O)F"],
)
def test_stereo_canonicalization_survives_relabeling(smiles):
    graph = StereoMolGraph.from_rdmol(Chem.AddHs(Chem.MolFromSmiles(smiles)))
    expected = _score(_initialize(graph), canon_atom_num(graph))
    rng = Random(42)
    for _ in range(5):
        shuffled = _shuffled(graph, rng)
        params = _initialize(shuffled)
        for pruning in (True, False):
            mapping = canon_atom_num(shuffled, prune_automorphisms=pruning)
            assert _score(params, mapping) == expected


def test_certificate_preserves_exact_stereo_and_placeholder():
    graph = StereoMolGraph(_graph(4, [(0, 1), (0, 2), (0, 3)]))
    graph.set_atom_stereo(Tetrahedral((0, 1, 2, 3, None), 1))
    params = _initialize(graph)
    mapping = {a: a + 1 for a in graph.atoms}
    original = _graph_certificate(params, mapping)
    assert _graph_certificate(params, {0: 1, 1: 3, 2: 2, 3: 4}) != original
    # Inversion plus an odd ligand exchange represents the same orientation.
    graph.set_atom_stereo(Tetrahedral((0, 2, 1, 3, None), -1))
    assert _graph_certificate(_initialize(graph), mapping) == original
    graph.set_atom_stereo(Tetrahedral((0, 1, 2, 3, None), None))
    assert _graph_certificate(_initialize(graph), mapping) != original


def test_certificate_distinguishes_descriptor_classes_with_undefined_parity():
    graph = StereoMolGraph(_graph(5, [(0, a) for a in range(1, 5)]))
    mapping = {a: a + 1 for a in graph.atoms}
    graph.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), None))
    original = _graph_certificate(_initialize(graph), mapping)
    graph.set_atom_stereo(SquarePlanar((0, 1, 2, 3, 4), None))
    assert _graph_certificate(_initialize(graph), mapping) != original


def test_stereo_with_lone_pair_and_negative_ids_survives_relabeling():
    graph = StereoMolGraph(_graph(4, [(0, 1), (0, 2), (0, 3)]))
    graph.set_atom_stereo(Tetrahedral((0, 1, 2, 3, None), 1))
    expected = _score(_initialize(graph), canon_atom_num(graph))
    rng = Random(123)
    for _ in range(10):
        shuffled = _shuffled(graph, rng)
        params = _initialize(shuffled)
        assert _score(params, canon_atom_num(shuffled)) == expected
        assert _exhaustive_score(params) == expected


@pytest.mark.parametrize("pruning", [False, True])
def test_disconnected_graph_uses_one_search_without_extraction(monkeypatch, pruning):
    graph = _graph(3, [(0, 1)])

    def reject_extraction(*args, **kwargs):
        raise AssertionError("Canonization must not split or extract components")

    monkeypatch.setattr(MolGraph, "connected_components", reject_extraction)
    monkeypatch.setattr(MolGraph, "subgraph", reject_extraction)
    mapping = canon_atom_num(graph, prune_automorphisms=pruning)
    params = _initialize(graph)
    assert _score(params, mapping) == _exhaustive_score(params)
    assert mapping[2] == 1  # The rare isolated atom precedes the larger component.


def test_disconnected_lone_pair_stereo_is_preserved():
    graph = StereoMolGraph()
    for atom, element in enumerate(("C", "F", "Cl", "Br", "O")):
        graph.add_atom(atom, element)
    for atom in (1, 2, 3):
        graph.add_bond(0, atom)
    descriptor = Tetrahedral((0, 1, 2, 3, None), 1)
    graph.set_atom_stereo(descriptor)
    expected = _graph_certificate(_initialize(graph), canon_atom_num(graph))
    rng = Random(17)
    for _ in range(10):
        shuffled = _shuffled(graph, rng)
        for pruning in (False, True):
            assert (
                _graph_certificate(
                    _initialize(shuffled),
                    canon_atom_num(shuffled, prune_automorphisms=pruning),
                )
                == expected
            )
    assert graph.get_atom_stereo(0) is descriptor


def test_disconnected_enantiomers_and_repeated_components_survive_shuffle():
    smiles = "F[C@](Cl)(Br)I.F[C@@](Cl)(Br)I.CC.CC.[He]"
    graph = StereoMolGraph.from_rdmol(Chem.AddHs(Chem.MolFromSmiles(smiles)))
    expected = _graph_certificate(_initialize(graph), canon_atom_num(graph))
    rng = Random(19)
    for _ in range(10):
        shuffled = _shuffled(graph, rng)
        assert (
            _graph_certificate(_initialize(shuffled), canon_atom_num(shuffled))
            == expected
        )


def test_stereo_linking_components_uses_whole_graph_search():
    graph = StereoMolGraph(_graph(5, []))
    graph.set_atom_stereo(Tetrahedral((0, 1, 2, 3, 4), 1))
    params = _initialize(graph)
    assert canon_atom_num(graph) == _search(params).best_mapping


@pytest.fixture(params=[CondensedReactionGraph, StereoCondensedReactionGraph])
def graph_type(request):
    return request.param


def _assert_active_first(graph, mapping):
    active = {
        atom
        for bond, attributes in graph.bonds_with_attributes.items()
        if attributes.get("reaction") in (Change.FORMED, Change.BROKEN)
        for atom in bond
    }
    assert {mapping[atom] for atom in active} == set(range(1, len(active) + 1))
    if isinstance(graph, StereoCondensedReactionGraph):
        stereo_atoms = {
            atom
            for records in (graph.atom_stereo_changes, graph.bond_stereo_changes)
            for changes in records.values()
            for stereo in changes.values()
            if stereo is not None
            for atom in stereo.atoms
            if atom is not None
        } - active
        assert {mapping[atom] for atom in stereo_atoms} == set(
            range(len(active) + 1, len(active) + len(stereo_atoms) + 1)
        )


def _bond_signature(graph, mapping):
    """Independent serialization of atom types, bond endpoints and change tags."""
    return (
        tuple(sorted((mapping[a], t) for a, t in zip(graph.atoms, graph.atom_types))),
        tuple(
            sorted(
                (
                    tuple(sorted(mapping[a] for a in bond)),
                    attributes["reaction"].value if attributes.get("reaction") else "",
                )
                for bond, attributes in graph.bonds_with_attributes.items()
            )
        ),
    )


@pytest.mark.parametrize("change", [Change.FORMED, Change.BROKEN])
@pytest.mark.parametrize("pruning", [False, True])
def test_active_hydrogen_precedes_inactive_heavy_atoms(graph_type, change, pruning):
    graph = graph_type()
    for atom, element in enumerate(("C", "H", "N", "O")):
        graph.add_atom(atom, element)
    graph.add_bond(0, 1, reaction=change)
    graph.add_bond(0, 2)
    graph.add_bond(2, 3)
    mapping = canon_atom_num(graph, prune_automorphisms=pruning)
    assert mapping[0] == 1
    assert mapping[1] == 2
    _assert_active_first(graph, mapping)


@pytest.mark.parametrize("connected", [False, True])
@pytest.mark.parametrize("pruning", [False, True])
def test_separated_active_regions_precede_inactive_frontier(
    graph_type, connected, pruning
):
    graph = graph_type()
    for atom in range(7):
        graph.add_atom(atom, "C")
    graph.add_formed_bond(0, 1)
    graph.add_broken_bond(3, 4)
    graph.add_bond(1, 2)
    graph.add_bond(4, 5)
    if connected:
        graph.add_bond(2, 3)  # Inactive bridge between two active regions.
    mapping = canon_atom_num(graph, prune_automorphisms=pruning)
    _assert_active_first(graph, mapping)


def test_fleeting_bonds_do_not_make_atoms_active(graph_type):
    graph = graph_type()
    for atom in range(4):
        graph.add_atom(atom, "C")
    graph.add_fleeting_bond(0, 1)
    graph.add_formed_bond(2, 3)
    assert _initialize(graph).active_atoms == frozenset((2, 3))
    _assert_active_first(graph, canon_atom_num(graph))


@pytest.mark.parametrize("change", [Change.FORMED, Change.BROKEN])
@pytest.mark.parametrize("element", ["C", "H"])
def test_changed_bond_connection_precedes_regular_connection(
    graph_type, change, element, monkeypatch
):
    monkeypatch.setattr(canon, "color_refine_crg", lambda graph: [0] * len(graph))
    monkeypatch.setattr(canon, "color_refine_scrg", lambda graph: [0] * len(graph))
    graph = graph_type()
    for atom in range(6):
        graph.add_atom(atom, element if atom == 1 else "C")
    graph.add_bond(0, 1, reaction=change)
    graph.add_formed_bond(2, 3)
    for atom in (2, 4, 5):
        graph.add_bond(0, atom)
    # All four active atoms compete in the first group. Atom 2 has greater
    # degree than atom 1 but is connected to the starting atom by a regular bond.
    for pruning in (False, True):
        mapping = canon_atom_num(graph, prune_automorphisms=pruning)
        assert mapping[0] == 1
        assert mapping[1 if element == "C" else 2] == 2
        assert _score(_initialize(graph), mapping) == _exhaustive_score(
            _initialize(graph)
        )
    expected = _bond_signature(graph, mapping)
    rng = Random(32)
    for _ in range(5):
        shuffled = _shuffled(graph, rng)
        for pruning in (False, True):
            mapping = canon_atom_num(shuffled, prune_automorphisms=pruning)
            assert _bond_signature(shuffled, mapping) == expected


@pytest.mark.parametrize("kind", ["atom", "bond"])
def test_stereo_change_type_directly_orders_equal_color_candidates(kind, monkeypatch):
    monkeypatch.setattr(canon, "color_refine_scrg", lambda graph: [0] * len(graph))
    graph = StereoCondensedReactionGraph()
    for start, phase in ((0, "broken"), (5, "formed")):
        for atom in range(start, start + 4):
            graph.add_atom(atom, "C")
        if kind == "atom":
            for atom in range(start + 1, start + 4):
                graph.add_bond(start, atom)
            graph.set_atom_stereo_change(
                **{
                    phase: Tetrahedral(
                        (start, start + 1, start + 2, start + 3, None), 1
                    )
                }
            )
        else:
            for atom in range(start, start + 3):
                graph.add_bond(atom, atom + 1)
            graph.set_bond_stereo_change(
                **{
                    phase: PlanarBond(
                        (start, None, start + 1, start + 2, start + 3, None), 0
                    )
                }
            )
    params = _initialize(graph)
    assert _candidate_label(1, params, {}) < _candidate_label(6, params, {})
    expected = _exhaustive_score(params)
    rng = Random(13)
    for _ in range(5):
        for pruning in (False, True):
            mapping = canon_atom_num(graph, prune_automorphisms=pruning)
            assert _score(_initialize(graph), mapping) == expected
        if _ == 0:
            assert min(mapping, key=mapping.get) in (0, 1, 2, 3)
        graph = _shuffled(graph, rng)


@pytest.mark.parametrize(
    "original_descriptor, descriptor",
    [
        (
            Tetrahedral((0, 1, 2, 3, None), 1),
            Tetrahedral((0, 1, 2, 3, None), -1),
        ),
        (
            Tetrahedral((0, 1, 2, 3, None), 1),
            Tetrahedral((0, 1, 2, 3, None), None),
        ),
        (
            Tetrahedral((0, 1, 2, 3, None), None),
            SquarePlanar((0, 1, 2, 3, None), None),
        ),
    ],
)
def test_stereo_change_descriptor_directly_affects_priority(
    original_descriptor, descriptor, monkeypatch
):
    monkeypatch.setattr(canon, "color_refine_scrg", lambda graph: [0] * len(graph))
    graph = StereoCondensedReactionGraph()
    for atom in range(4):
        graph.add_atom(atom, "C")
    for atom in (1, 2, 3):
        graph.add_bond(0, atom)
    graph.set_atom_stereo_change(formed=original_descriptor)
    mapping = {1: 1, 2: 2, 3: 3}
    original = _candidate_label(0, _initialize(graph), mapping)
    graph.set_atom_stereo_change(formed=descriptor)
    assert _candidate_label(0, _initialize(graph), mapping) != original


def test_bond_stereo_change_orientation_directly_affects_priority(monkeypatch):
    monkeypatch.setattr(canon, "color_refine_scrg", lambda graph: [0] * len(graph))
    graph = StereoCondensedReactionGraph()
    for atom in range(4):
        graph.add_atom(atom, "C")
    for atom in range(3):
        graph.add_bond(atom, atom + 1)
    graph.set_bond_stereo_change(formed=PlanarBond((0, None, 1, 2, 3, None), 0))
    mapping = {0: 1, 2: 2, 3: 3}
    original = _candidate_label(1, _initialize(graph), mapping)
    graph.set_bond_stereo_change(formed=PlanarBond((None, 0, 1, 2, 3, None), 0))
    assert _candidate_label(1, _initialize(graph), mapping) != original


@pytest.mark.parametrize("n", [0, 1, 3])
def test_reaction_without_bond_changes_has_no_active_atoms(graph_type, n):
    graph = graph_type()
    for atom in range(n):
        graph.add_atom(atom, "C")
    for atom in range(n - 1):
        graph.add_bond(atom, atom + 1)
    assert not _initialize(graph).active_atoms
    assert set(canon_atom_num(graph).values()) == set(range(1, n + 1))


def test_random_reactions_match_exhaustive_search_and_shuffle(graph_type):
    rng = Random(2026)
    changes = (None, Change.FORMED, Change.BROKEN, Change.FLEETING)
    for _ in range(40):
        graph = graph_type()
        n = rng.randrange(1, 7)
        for atom in range(n):
            graph.add_atom(atom, rng.choice(("C", "N", "O", "H")))
        for a, b in combinations(range(n), 2):
            if rng.random() < 0.4:
                change = rng.choice(changes)
                graph.add_bond(a, b, **({} if change is None else {"reaction": change}))
        params = _initialize(graph)
        expected = _exhaustive_score(params)
        reference = None
        for pruning in (False, True):
            mapping = canon_atom_num(graph, prune_automorphisms=pruning)
            assert _score(params, mapping) == expected
            _assert_active_first(graph, mapping)
            signature = _bond_signature(graph, mapping)
            reference = signature if reference is None else reference
            assert signature == reference
        for _ in range(3):
            shuffled = _shuffled(graph, rng)
            mapping = canon_atom_num(shuffled)
            _assert_active_first(shuffled, mapping)
            assert _bond_signature(shuffled, mapping) == reference


def test_exact_certificate_rejects_exchanging_formed_and_broken_bonds(
    graph_type, monkeypatch
):
    # Make all colors collide: reaction identity must come from exact bond tags.
    monkeypatch.setattr(canon, "color_refine_crg", lambda graph: [0] * len(graph))
    monkeypatch.setattr(canon, "color_refine_scrg", lambda graph: [0] * len(graph))
    graph = graph_type()
    for atom in range(4):
        graph.add_atom(atom, "C")
    graph.add_formed_bond(0, 1)
    graph.add_bond(1, 2)
    graph.add_broken_bond(2, 3)
    graph.add_bond(3, 0)
    params = _initialize(graph)
    identity = {0: 1, 1: 2, 2: 3, 3: 4}
    rotated = {0: 3, 1: 4, 2: 1, 3: 2}
    assert _graph_certificate(params, identity) != _graph_certificate(params, rotated)
    assert _bond_signature(graph, canon_atom_num(graph)) == _bond_signature(
        graph, canon_atom_num(graph, prune_automorphisms=False)
    )


@pytest.mark.parametrize("kind", ["atom", "bond"])
def test_reaction_stereo_changes_and_lone_pairs_survive_shuffle(kind):
    graph = StereoCondensedReactionGraph()
    if kind == "atom":
        for atom, element in enumerate(("N", "C", "F", "Cl", "He")):
            graph.add_atom(atom, element)
        for atom in (1, 2, 3):
            graph.add_bond(0, atom)
        graph.add_formed_bond(1, 4)
        graph.set_atom_stereo_change(
            broken=Tetrahedral((0, 1, 2, 3, None), -1),
            formed=Tetrahedral((0, 1, 2, 3, None), 1),
        )
        assert _initialize(graph).active_atoms == frozenset((1, 4))
    else:
        for atom, element in enumerate(("C", "N", "N", "C", "F")):
            graph.add_atom(atom, element)
        for a, b in ((0, 1), (1, 2), (2, 3)):
            graph.add_bond(a, b)
        graph.add_formed_bond(0, 4)
        graph.set_bond_stereo_change(
            broken=PlanarBond((0, None, 1, 2, 3, None), 0),
            formed=PlanarBond((None, 0, 1, 2, 3, None), 0),
        )
        assert _initialize(graph).active_atoms == frozenset((0, 4))
    expected = _graph_certificate(_initialize(graph), canon_atom_num(graph))
    rng = Random(45)
    for _ in range(10):
        shuffled = _shuffled(graph, rng)
        for pruning in (False, True):
            mapping = canon_atom_num(shuffled, prune_automorphisms=pruning)
            _assert_active_first(shuffled, mapping)
            assert _graph_certificate(_initialize(shuffled), mapping) == expected


def test_exact_certificate_distinguishes_stereo_changes(monkeypatch):
    monkeypatch.setattr(canon, "color_refine_scrg", lambda graph: [0] * len(graph))
    graph = StereoCondensedReactionGraph()
    for atom in range(5):
        graph.add_atom(atom, "C")
    for atom in range(1, 5):
        graph.add_bond(0, atom)
    mapping = {atom: atom + 1 for atom in graph.atoms}
    graph.set_atom_stereo_change(formed=Tetrahedral((0, 1, 2, 3, 4), 1))
    assert not _initialize(graph).active_atoms
    original = _graph_certificate(_initialize(graph), mapping)
    graph.set_atom_stereo_change(formed=Tetrahedral((0, 1, 2, 3, 4), -1))
    assert _graph_certificate(_initialize(graph), mapping) != original
    graph.set_atom_stereo_change(broken=Tetrahedral((0, 1, 2, 3, 4), 1))
    assert _graph_certificate(_initialize(graph), mapping) != original


@pytest.mark.parametrize("kind", ["atom", "bond"])
@pytest.mark.parametrize("connected", [False, True])
@pytest.mark.parametrize("bond_change", [False, True])
@pytest.mark.parametrize("phase", ["formed", "broken", "fleeting"])
def test_bond_changes_then_stereo_changes_then_rest(
    kind, connected, bond_change, phase
):
    graph = StereoCondensedReactionGraph()
    for atom, element in enumerate(
        ("C", "H", "N", "C", "F", "H", "O", "N", "C", "F", "H", "C")
    ):
        graph.add_atom(atom, element)
    graph.add_bond(0, 11)  # Lower-priority frontier after visiting atom 0.
    graph.add_bond(0, 1)
    if bond_change:
        graph.add_formed_bond(0, 1)
        graph.add_broken_bond(1, 3)  # Atom 3 belongs to both change groups.
    for start in (2, 7):
        if kind == "atom":
            for ligand in range(start + 1, start + 4):
                graph.add_bond(start, ligand)
            descriptor = Tetrahedral((start, start + 1, start + 2, start + 3, None), 1)
            graph.set_atom_stereo_change(**{phase: descriptor})
        else:
            for atom in range(start, start + 3):
                graph.add_bond(atom, atom + 1)
            descriptor = PlanarBond(
                (start, None, start + 1, start + 2, start + 3, None), 0
            )
            graph.set_bond_stereo_change(**{phase: descriptor})
    if connected:
        for a, b in ((0, 2), (3, 6), (6, 8)):
            graph.add_bond(a, b)

    params = _initialize(graph)
    assert params.stereo_change_atoms == frozenset((2, 3, 4, 5, 7, 8, 9, 10))
    expected_score = _exhaustive_score(params)
    expected_certificate = expected_score[1]
    rng = Random(17)
    for pruning in (False, True):
        mapping = canon_atom_num(graph, prune_automorphisms=pruning)
        _assert_active_first(graph, mapping)
        assert _score(params, mapping) == expected_score
        assert _graph_certificate(params, mapping) == expected_certificate
        # Stereo-change hydrogens must precede every atom in the final group.
        assert max(mapping[5], mapping[10]) < min(mapping[6], mapping[11])
    for _ in range(3):
        graph = _shuffled(graph, rng)
        for pruning in (False, True):
            mapping = canon_atom_num(graph, prune_automorphisms=pruning)
            _assert_active_first(graph, mapping)
            assert _score(_initialize(graph), mapping) == expected_score
