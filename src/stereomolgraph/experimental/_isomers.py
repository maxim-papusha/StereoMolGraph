from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np

from stereomolgraph import AtomId, Bond, StereoCondensedReactionGraph, StereoMolGraph
from stereomolgraph.algorithms.color_refine import color_refine_smg
from stereomolgraph.graphs.crg import Change

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator


def unique_generator(input_generator: Iterator) -> Iterator:
    """
    A generator that yields unique objects from another generator.

    Args:
        input_generator: A generator yielding hashable objects

    Yields:
        Only the first occurrence of each unique object from the
        input generator
    """
    seen_hash = set()
    for item in input_generator:
        item_hash = hash(item)
        if item_hash not in seen_hash:
            seen_hash.add(item_hash)
            yield item


@dataclass(frozen=True)
class _StereoSite:
    kind: Literal["atom", "bond"]
    key: AtomId | Bond


@dataclass(frozen=True)
class _StereoSiteInfo:
    site: _StereoSite
    signature: tuple[object, ...]
    atom_scope: frozenset[AtomId]
    neighborhood: frozenset[AtomId]
    ambiguous: bool
    isomers: tuple[object, ...]


def _unique_ordered(values: Iterable[AtomId]) -> tuple[AtomId, ...]:
    return tuple(dict.fromkeys(values))


def _unique_bonds(values: Iterable[Bond]) -> tuple[Bond, ...]:
    return tuple(dict.fromkeys(Bond(bond) for bond in values))


def _color_map(graph: StereoMolGraph) -> dict[AtomId, int]:
    initial_color_array = np.array(graph.atom_types, dtype=np.int64)
    color_array = color_refine_smg(graph, atom_labels=initial_color_array)
    return {atom: int(color) for atom, color in zip(graph.atoms, color_array)}


def _color_token(
    atom: None | AtomId,
    color_map: dict[AtomId, int],
) -> tuple[int, int]:
    if atom is None:
        return (1, 0)
    return (0, color_map[atom])


def _site_stereo(graph: StereoMolGraph, site: _StereoSite):
    if site.kind == "atom":
        stereo = graph.get_atom_stereo(site.key)
    else:
        stereo = graph.get_bond_stereo(site.key)

    if stereo is None:
        raise ValueError(f"Stereo site {site} is not present in the graph")
    return stereo


def _set_site_stereo(graph: StereoMolGraph, site: _StereoSite, stereo: object):
    if site.kind == "atom":
        graph.set_atom_stereo(stereo)
    else:
        graph.set_bond_stereo(stereo)


def _anchor_atoms(stereo: object) -> tuple[AtomId, ...]:
    if hasattr(stereo, "central_atom"):
        return (stereo.central_atom,)
    if hasattr(stereo, "bond"):
        return tuple(sorted(stereo.bond))
    raise TypeError(f"Unsupported stereo type: {type(stereo)!r}")


def _variable_atoms(stereo: object) -> tuple[None | AtomId, ...]:
    if hasattr(stereo, "central_atom"):
        return tuple(stereo.atoms[1:])

    atoms = tuple(stereo.atoms)
    if len(atoms) >= 4 and frozenset(atoms[2:4]) == stereo.bond:
        return atoms[:2] + atoms[4:]

    anchor_set = set(stereo.bond)
    return tuple(atom for atom in atoms if atom not in anchor_set)


def _stereo_sort_key(stereo: object) -> tuple[object, ...]:
    if hasattr(stereo, "canonical_form"):
        canonical = stereo.canonical_form()
    else:
        canonical = (getattr(stereo, "parity", None), tuple(stereo.atoms))
    return (stereo.__class__.__name__, canonical)


def _site_signature(
    stereo: object,
    color_map: dict[AtomId, int],
) -> tuple[tuple[object, ...], bool, frozenset[AtomId]]:
    anchors = _anchor_atoms(stereo)
    variable_atoms = _variable_atoms(stereo)
    anchor_colors = tuple(sorted(_color_token(atom, color_map) for atom in anchors))
    variable_colors = tuple(
        sorted(_color_token(atom, color_map) for atom in variable_atoms)
    )
    ambiguous = len(variable_colors) != len(set(variable_colors))
    atom_scope = frozenset(
        atom for atom in (*anchors, *variable_atoms) if atom is not None
    )
    signature = (stereo.__class__.__name__, anchor_colors, variable_colors)
    return signature, ambiguous, atom_scope


def _site_neighborhood(
    graph: StereoMolGraph,
    atom_scope: frozenset[AtomId],
) -> frozenset[AtomId]:
    neighborhood = set(atom_scope)
    for atom in atom_scope:
        neighborhood.update(graph.bonded_to(atom))
    return frozenset(neighborhood)


def _collect_stereo_sites(
    graph: StereoMolGraph,
    atoms: None | Iterable[AtomId] = None,
    bonds: None | Iterable[Bond] = None,
) -> tuple[_StereoSiteInfo, ...]:
    color_map = _color_map(graph)
    infos: list[_StereoSiteInfo] = []

    atom_ids = _unique_ordered(graph.atoms if atoms is None else atoms)
    for atom in atom_ids:
        stereo = graph.get_atom_stereo(atom)
        if stereo is None or stereo.parity is not None:
            continue

        signature, ambiguous, atom_scope = _site_signature(stereo, color_map)
        site = _StereoSite("atom", atom)
        infos.append(
            _StereoSiteInfo(
                site=site,
                signature=signature,
                atom_scope=atom_scope,
                neighborhood=_site_neighborhood(graph, atom_scope),
                ambiguous=ambiguous,
                isomers=tuple(sorted(stereo.get_isomers(), key=_stereo_sort_key)),
            )
        )

    bond_ids = _unique_bonds(graph.bonds if bonds is None else bonds)
    for bond in bond_ids:
        stereo = graph.get_bond_stereo(bond)
        if stereo is None or stereo.parity is not None:
            continue

        signature, ambiguous, atom_scope = _site_signature(stereo, color_map)
        site = _StereoSite("bond", bond)
        infos.append(
            _StereoSiteInfo(
                site=site,
                signature=signature,
                atom_scope=atom_scope,
                neighborhood=_site_neighborhood(graph, atom_scope),
                ambiguous=ambiguous,
                isomers=tuple(sorted(stereo.get_isomers(), key=_stereo_sort_key)),
            )
        )

    return tuple(infos)


def _add_edge(
    adjacency: dict[_StereoSite, set[_StereoSite]],
    source: _StereoSite,
    target: _StereoSite,
):
    if source != target:
        adjacency[source].add(target)


def _dependency_graph(
    infos: tuple[_StereoSiteInfo, ...],
) -> dict[_StereoSite, set[_StereoSite]]:
    adjacency = {info.site: set() for info in infos}

    signature_groups: dict[tuple[object, ...], list[_StereoSiteInfo]] = {}
    for info in infos:
        signature_groups.setdefault(info.signature, []).append(info)

    for group in signature_groups.values():
        if len(group) < 2:
            continue
        for info_a, info_b in itertools.permutations(group, 2):
            _add_edge(adjacency, info_a.site, info_b.site)

    for info_a, info_b in itertools.combinations(infos, 2):
        if (
            info_a.neighborhood.isdisjoint(info_b.atom_scope)
            and info_b.neighborhood.isdisjoint(info_a.atom_scope)
        ):
            continue

        if info_a.ambiguous:
            _add_edge(adjacency, info_a.site, info_b.site)
        if info_b.ambiguous:
            _add_edge(adjacency, info_b.site, info_a.site)

    return adjacency


def _strongly_connected_components(
    adjacency: dict[_StereoSite, set[_StereoSite]],
) -> tuple[tuple[_StereoSite, ...], ...]:
    visited: set[_StereoSite] = set()
    order: list[_StereoSite] = []

    def dfs(node: _StereoSite):
        visited.add(node)
        for neighbor in adjacency[node]:
            if neighbor not in visited:
                dfs(neighbor)
        order.append(node)

    for node in adjacency:
        if node not in visited:
            dfs(node)

    reverse = {node: set() for node in adjacency}
    for node, neighbors in adjacency.items():
        for neighbor in neighbors:
            reverse[neighbor].add(node)

    components: list[tuple[_StereoSite, ...]] = []
    visited.clear()

    def reverse_dfs(node: _StereoSite, component: list[_StereoSite]):
        visited.add(node)
        component.append(node)
        for neighbor in reverse[node]:
            if neighbor not in visited:
                reverse_dfs(neighbor, component)

    for node in reversed(order):
        if node in visited:
            continue
        component: list[_StereoSite] = []
        reverse_dfs(node, component)
        components.append(tuple(component))

    return tuple(components)


def _component_order(
    components: tuple[tuple[_StereoSite, ...], ...],
    adjacency: dict[_StereoSite, set[_StereoSite]],
    info_by_site: dict[_StereoSite, _StereoSiteInfo],
) -> tuple[tuple[_StereoSite, ...], ...]:
    if not components:
        return ()

    site_order = {
        info.site: idx for idx, info in enumerate(info_by_site.values())
    }
    component_indices = {
        site: idx
        for idx, component in enumerate(components)
        for site in component
    }
    component_edges = {idx: set() for idx in range(len(components))}
    component_reverse = {idx: set() for idx in range(len(components))}
    for site, neighbors in adjacency.items():
        source_idx = component_indices[site]
        for neighbor in neighbors:
            target_idx = component_indices[neighbor]
            if source_idx == target_idx:
                continue
            component_edges[source_idx].add(target_idx)
            component_reverse[target_idx].add(source_idx)

    def component_rank(component_index: int) -> int:
        return min(site_order[site] for site in components[component_index])

    deferred: set[int] = set()
    for idx, component in enumerate(components):
        if (
            len(component) == 1
            and not info_by_site[component[0]].ambiguous
            and not component_edges[idx]
            and not component_reverse[idx]
        ):
            deferred.add(idx)

    indegree = {idx: 0 for idx in range(len(components)) if idx not in deferred}
    for source_idx, targets in component_edges.items():
        if source_idx in deferred:
            continue
        for target_idx in targets:
            if target_idx in deferred:
                continue
            indegree[target_idx] += 1

    ready = sorted(
        (idx for idx, degree in indegree.items() if degree == 0),
        key=component_rank,
    )
    ordered: list[tuple[_StereoSite, ...]] = []
    while ready:
        component_idx = ready.pop(0)
        ordered.append(
            tuple(sorted(components[component_idx], key=lambda site: site_order[site]))
        )
        for target_idx in sorted(component_edges[component_idx], key=component_rank):
            if target_idx in deferred:
                continue
            indegree[target_idx] -= 1
            if indegree[target_idx] == 0:
                ready.append(target_idx)
                ready.sort(key=component_rank)

    ordered.extend(
        tuple(sorted(components[idx], key=lambda site: site_order[site]))
        for idx in sorted(deferred, key=component_rank)
    )
    return tuple(ordered)


def _dependency_components(
    graph: StereoMolGraph,
    atoms: None | Iterable[AtomId] = None,
    bonds: None | Iterable[Bond] = None,
) -> tuple[_StereoSiteInfo, ...] | tuple[tuple[_StereoSite, ...], ...]:
    infos = _collect_stereo_sites(graph, atoms=atoms, bonds=bonds)
    info_by_site = {info.site: info for info in infos}
    adjacency = _dependency_graph(infos)
    components = _strongly_connected_components(adjacency)
    return _component_order(components, adjacency, info_by_site)


def _refinement_signature(graph: StereoMolGraph) -> tuple[int, ...]:
    initial_color_array = np.array(graph.atom_types, dtype=np.int64)
    colors = color_refine_smg(graph, atom_labels=initial_color_array)
    return tuple(int(color) for color in colors)


def _append_unique_graph(
    buckets: dict[tuple[int, ...], list[StereoMolGraph]],
    unique_graphs: list[StereoMolGraph],
    candidate: StereoMolGraph,
):
    signature = _refinement_signature(candidate)
    bucket = buckets.setdefault(signature, [])
    if any(existing == candidate for existing in bucket):
        return
    bucket.append(candidate)
    unique_graphs.append(candidate)


def _collapse_enantiomers(
    states: Iterable[StereoMolGraph],
) -> Iterator[StereoMolGraph]:
    unique: list[StereoMolGraph] = []
    for state in states:
        mirrored = state.enantiomer()
        if any(existing == state or existing == mirrored for existing in unique):
            continue
        unique.append(state)
        yield state


def _generate_stereoisomers_bruteforce(
    graph: StereoMolGraph,
    enantiomers: bool = True,
    atoms: None | Iterable[AtomId] = None,
    bonds: None | Iterable[Bond] = None,
) -> Iterator[StereoMolGraph]:
    if atoms is None:
        atom_stereos = (
            stereo.get_isomers()
            for a in graph.atoms
            if ((stereo := graph.get_atom_stereo(a)) and stereo.parity is None)
        )
    else:
        atom_stereos = (
            stereo.get_isomers()
            for a in atoms
            if (stereo := graph.get_atom_stereo(a)) is not None
        )

    if bonds is None:
        bond_stereos = (
            stereo.get_isomers()
            for b in graph.bonds
            if ((stereo := graph.get_bond_stereo(b)) and stereo.parity is None)
        )
    else:
        bond_stereos = (
            stereo.get_isomers()
            for b in bonds
            if (stereo := graph.get_bond_stereo(b)) is not None
        )

    seen: set[StereoMolGraph] = set()
    enantiomers_seen: set[StereoMolGraph] = set()

    for a_stereos, b_stereos in itertools.product(
        itertools.product(*atom_stereos), itertools.product(*bond_stereos)
    ):
        stereoisomer = graph.copy()
        for a_stereo in a_stereos:
            stereoisomer.set_atom_stereo(a_stereo)
        for b_stereo in b_stereos:
            stereoisomer.set_bond_stereo(b_stereo)

        if stereoisomer not in seen:
            seen.add(stereoisomer)
            yield stereoisomer

            if enantiomers:
                enantiomer = stereoisomer.enantiomer()

                if enantiomer != stereoisomer and enantiomer not in enantiomers_seen:
                    enantiomers_seen.add(enantiomer)
                    yield enantiomer


def generate_stereoisomers_staged(
    graph: StereoMolGraph,
    enantiomers: bool = True,
    atoms: None | Iterable[AtomId] = None,
    bonds: None | Iterable[Bond] = None,
) -> Iterator[StereoMolGraph]:
    """Generates unique stereoisomers by staged activation of dependent blocks.

    This is a first SCC-based implementation of the dependency-driven outline:
    unresolved stereodescriptors are grouped using Stereo-WL colors as a proxy
    for dependence, then activated block-wise with pruning after each block.
    The schedule is intentionally conservative, so larger components still
    preserve correctness even if dependence is over-approximated.
    """
    base_graph = graph.copy()
    infos = _collect_stereo_sites(base_graph, atoms=atoms, bonds=bonds)
    if not infos:
        if enantiomers:
            yield base_graph
        else:
            yield from _collapse_enantiomers([base_graph])
        return

    info_by_site = {info.site: info for info in infos}
    adjacency = _dependency_graph(infos)
    components = _component_order(
        _strongly_connected_components(adjacency),
        adjacency,
        info_by_site,
    )

    states = [base_graph]
    for component in components:
        next_states: list[StereoMolGraph] = []
        buckets: dict[tuple[int, ...], list[StereoMolGraph]] = {}
        component_isomers = [info_by_site[site].isomers for site in component]
        for state in states:
            for assignment in itertools.product(*component_isomers):
                candidate = state.copy()
                for site, stereo in zip(component, assignment):
                    _set_site_stereo(candidate, site, stereo)
                _append_unique_graph(buckets, next_states, candidate)
        states = next_states

    if enantiomers:
        yield from states
    else:
        yield from _collapse_enantiomers(states)


def generate_stereoisomers(
    graph: StereoMolGraph,
    enantiomers: bool = True,
    atoms: None | Iterable[AtomId] = None,
    bonds: None | Iterable[Bond] = None,
    strategy: Literal["bruteforce", "staged"] = "bruteforce",
) -> Iterator[StereoMolGraph]:
    """Generates unique stereoisomers of a StereoMolGraph.

    Args:
        graph: Input graph.
        enantiomers: If True, enumerate both members of enantiomeric pairs.
        atoms: Optional subset of atoms to consider for stereoisomerism.
        bonds: Optional subset of bonds to consider for stereoisomerism.
        strategy: Enumeration strategy. ``"bruteforce"`` preserves the
            existing exhaustive generator, while ``"staged"`` activates
            stereodescriptors SCC-by-SCC with pruning between stages.

    Yields:
        StereoMolGraph: Each unique stereoisomer.
    """
    if strategy == "bruteforce":
        yield from _generate_stereoisomers_bruteforce(
            graph,
            enantiomers=enantiomers,
            atoms=atoms,
            bonds=bonds,
        )
        return

    if strategy == "staged":
        yield from generate_stereoisomers_staged(
            graph,
            enantiomers=enantiomers,
            atoms=atoms,
            bonds=bonds,
        )
        return

    raise ValueError(
        f"Unknown stereoisomer generation strategy: {strategy!r}"
    )


def generate_fleeting_stereoisomers(
    graph: StereoCondensedReactionGraph,
    enantiomers: bool = True,
    atoms: None | Iterable[AtomId] = None,
    bonds: None | Iterable[Bond] = None,
) -> Iterator[StereoCondensedReactionGraph]:
    """Generates fleeting stereoisomers of a reaction graph.

    Only includes stereocenters which have a parity of None for the
    fleeting change.
    If a parity is set, it is not changed.

    Args:
        graph: The reaction graph to generate isomers from
        enantiomers: If True, both enantiomers are included (default: True)
        atoms: Optional subset of atoms to consider for stereoisomerism
        bonds: Optional subset of bonds to consider for stereoisomerism

    Yields:
        StereoCondensedReactionGraph: Each unique fleeting stereoisomer
    """
    if atoms is None:
        atom_stereos = (
            stereo.get_isomers()
            for a in graph.atoms
            if (
                (stereo_change_dict := graph.get_atom_stereo_change(a))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        )
    else:
        atom_stereos = (
            stereo.get_isomers()
            for a in atoms
            if (
                (stereo_change_dict := graph.get_atom_stereo_change(a))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        )

    if bonds is None:
        bond_stereos = (
            stereo.get_isomers()
            for b in graph.bonds
            if (
                (stereo_change_dict := graph.get_bond_stereo_change(b))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        )
    else:
        bond_stereos = (
            stereo.get_isomers()
            for b in bonds
            if (
                (stereo_change_dict := graph.get_bond_stereo_change(b))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        )

    seen_isomers: set[StereoCondensedReactionGraph] = set()
    seen_enantiomers: set[StereoCondensedReactionGraph] = set()

    for a_stereos, b_stereos in itertools.product(
        itertools.product(*atom_stereos), itertools.product(*bond_stereos)
    ):
        stereoisomer = graph.copy()
        for a_stereo in a_stereos:
            stereoisomer.set_atom_stereo_change(fleeting=a_stereo)
        for b_stereo in b_stereos:
            stereoisomer.set_bond_stereo_change(fleeting=b_stereo)

        if stereoisomer not in seen_isomers:
            seen_isomers.add(stereoisomer)
            yield stereoisomer

            if enantiomers:
                enantiomer = stereoisomer.enantiomer()
                if enantiomer != stereoisomer and enantiomer not in seen_enantiomers:
                    seen_enantiomers.add(enantiomer)
                    yield enantiomer
