from __future__ import annotations

import itertools
from collections import defaultdict, deque
from collections.abc import Hashable, Iterable, Iterator, Mapping
from heapq import heappop, heappush
from typing import TYPE_CHECKING, TypeVar

from stereomolgraph import StereoCondensedReactionGraph, StereoMolGraph
from stereomolgraph.algorithms.color_refine import color_refine_smg
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.graphs.scrg import Change

T = TypeVar("T")


if TYPE_CHECKING:
    from typing import Optional

    from stereomolgraph.graphs.smg import (
        AtomId,
        AtomStereo,
        Bond,
        BondStereo,
        StereoMolGraph,
    )

    Atom_Or_Bond = AtomId | Bond


def generate_stereo_permutations(
    smg: StereoMolGraph, atoms: Iterable[int], bonds: Iterable[frozenset[int]]
) -> Iterator[StereoMolGraph]: ...


def unique_generator(
    source: Iterator[T],
) -> Iterator[T]:
    """Yield only the first occurrence of each element in *source*.

    A fast *hash_func* groups candidates into buckets, and *eq_func* resolves
    hash collisions by pairwise comparison. The function keeps only the minimal
    state required to guard against duplicates and releases each unique element
    lazily as soon as it is identified.
    """

    buckets: dict[Hashable, list[T]] = {}

    for candidate in source:
        key = hash(candidate)
        bucket = buckets.get(key)

        if bucket is None:
            buckets[key] = [candidate]
            yield candidate
            continue

        for seen in bucket:
            # TODO: use refined colors!
            if candidate == seen:
                break
        else:
            bucket.append(candidate)
            yield candidate


def strongly_connected_components(
    directed_multigraph: Mapping[int, Iterable[int]],
) -> list[frozenset[int]]:
    """Return SCCs in edge-respecting order, preferring larger components.

    The graph is first decomposed with Tarjan's algorithm. The resulting
    condensation DAG is then traversed with a size-aware Kahn topological
    ordering so that whenever multiple SCCs are ready, the largest is emitted
    first. This keeps traversal aligned with edge directions while favouring
    richer components.
    """

    adjacency: dict[int, set[int]] = {}
    all_nodes: set[int] = set()

    for node, neighbours in directed_multigraph.items():
        node_neighbors = adjacency.setdefault(node, set())
        for neighbour in neighbours:
            node_neighbors.add(neighbour)
            all_nodes.add(neighbour)
        all_nodes.add(node)

    for node in all_nodes:
        adjacency.setdefault(node, set())

    index_counter = 0
    indices: dict[int, int] = {}
    lowlinks: dict[int, int] = {}
    stack: list[int] = []
    on_stack: set[int] = set()
    components: list[frozenset[int]] = []

    def strongconnect(node: int) -> None:
        nonlocal index_counter
        indices[node] = index_counter
        lowlinks[node] = index_counter
        index_counter += 1

        stack.append(node)
        on_stack.add(node)

        for neighbour in adjacency[node]:
            if neighbour not in indices:
                strongconnect(neighbour)
                lowlinks[node] = min(lowlinks[node], lowlinks[neighbour])
            elif neighbour in on_stack:
                lowlinks[node] = min(lowlinks[node], indices[neighbour])

        if lowlinks[node] == indices[node]:
            component: set[int] = set()
            while True:
                w = stack.pop()
                on_stack.remove(w)
                component.add(w)
                if w == node:
                    break
            components.append(frozenset(component))

    for node in adjacency:
        if node not in indices:
            strongconnect(node)

    if not components:
        return []

    component_index: dict[int, int] = {}
    for idx, component in enumerate(components):
        for member in component:
            component_index[member] = idx

    dag_successors: list[set[int]] = [set() for _ in components]
    indegree: list[int] = [0] * len(components)

    for node, neighbours in adjacency.items():
        src_idx = component_index[node]
        for neighbour in neighbours:
            dst_idx = component_index[neighbour]
            if src_idx == dst_idx:
                continue
            if dst_idx not in dag_successors[src_idx]:
                dag_successors[src_idx].add(dst_idx)
                indegree[dst_idx] += 1

    heap: list[tuple[int, tuple[int, ...], int]] = []
    for idx, component in enumerate(components):
        if indegree[idx] == 0:
            heappush(
                heap,
                (-len(component), tuple(sorted(component)), idx),
            )

    ordered_components: list[frozenset[int]] = []
    while heap:
        _, _, idx = heappop(heap)
        ordered_components.append(components[idx])
        for successor in dag_successors[idx]:
            indegree[successor] -= 1
            if indegree[successor] == 0:
                successor_component = components[successor]
                heappush(
                    heap,
                    (
                        -len(successor_component),
                        tuple(sorted(successor_component)),
                        successor,
                    ),
                )

    return ordered_components


Color = int


def fast_stereoisomer(
    graph: StereoMolGraph,
    atoms: Optional[Iterable[AtomId]] = None,
    bonds: Optional[Iterable[Bond]] = None,
) -> Iterator[StereoMolGraph]:
    atoms = graph.atoms if atoms is None else tuple(atoms)
    bonds = graph.bonds if bonds is None else tuple(bonds)

    id_index_dict: dict[AtomId, int] = {
        atom: idx for idx, atom in enumerate(graph.atoms)
    }
    init_colors: list[int] = [int(c) for c in color_refine_smg(graph)]

    atom_stereodescriptor_colors: dict[Color, set[AtomId]] = defaultdict(set)
    bond_stereodescriptor_colors: dict[Color, set[Bond]] = defaultdict(set)

    for atom, color in zip(graph.atoms, init_colors):
        if (
            a_stereo := graph.get_atom_stereo(atom)
        ) is not None and a_stereo.parity is None:
            atom_stereodescriptor_colors[color].add(atom)

    for bond, b_stereo in graph.bond_stereo.items():
        if b_stereo.parity is None:
            bond_color = hash(
                frozenset(init_colors[id_index_dict[atom]] for atom in bond)
            )
            bond_stereodescriptor_colors[bond_color].add(bond)

    unique_stereodescriptors: set[Atom_Or_Bond] = set()
    for color, atoms in atom_stereodescriptor_colors.items():
        if len(atoms) == 1:
            unique_stereodescriptors.add(min(atoms))
    for color, bonds in bond_stereodescriptor_colors.items():
        if len(bonds) == 1:
            unique_stereodescriptors.add(min(bonds))

    dcg: dict[Color, set[Color]] = defaultdict(set)
    for atom in graph.atoms:
        if atom in unique_stereodescriptors:
            continue
        a_stereo = graph.get_atom_stereo(atom)
        nbrs: Iterable[int] = (
            a_stereo.atoms if a_stereo else graph.bonded_to(atom)
        )
        for nbr in nbrs:
            if nbr not in unique_stereodescriptors:
                dcg[init_colors[id_index_dict[atom]]].add(
                    init_colors[id_index_dict[nbr]]
                )

    for bond, b_stereo in graph.bond_stereo.items():
        if bond in unique_stereodescriptors:
            continue
        for nbr in b_stereo.atoms:
            bond_color: int = hash(
                frozenset(init_colors[id_index_dict[atom]] for atom in bond)
            )
            if nbr not in bond:
                dcg[bond_color].add(init_colors[id_index_dict[nbr]])

    sccs: list[frozenset[int]] = strongly_connected_components(dcg)

    def set_atom_stereos(
        smg: StereoMolGraph,
        atom_stereos: Iterable[AtomStereo],
    ) -> StereoMolGraph:
        smg = smg.copy()
        for atom_stereo in atom_stereos:
            smg.set_atom_stereo(atom_stereo)
        return smg

    def set_bond_stereos(
        smg: StereoMolGraph,
        bond_stereos: Iterable[BondStereo],
    ) -> StereoMolGraph:
        smg = smg.copy()
        for bond_stereo in bond_stereos:
            smg.set_bond_stereo(bond_stereo)
        return smg

    def perm_gen(
        smg: StereoMolGraph,
        last_colors: Iterable[int],
        sccs: list[frozenset[int]],
        unique_stereodescriptors: set[int | frozenset[int]],
    ) -> Iterator[StereoMolGraph]:
        sccs = sccs.copy()
        scc: frozenset[int | frozenset[int]] = sccs.pop()

        if next(iter(scc)) in atom_stereodescriptor_colors:
            a_stereos = [
                a_stereo.get_isomers()  # type: ignore
                for color in scc
                for atom in atom_stereodescriptor_colors[color]
                if (a_stereo := smg.get_atom_stereo(atom)).parity is None  # type: ignore
            ]
            perm_generator = itertools.product(*a_stereos)
            smg_gen: Iterator[StereoMolGraph] = unique_generator(
                set_atom_stereos(smg, perm_stereos)
                for perm_stereos in perm_generator
            )

            if not sccs:
                yield from smg_gen
            else:
                yield from perm_gen(
                    smg,
                    last_colors,
                    sccs,
                    unique_stereodescriptors,
                )
        elif next(iter(scc)) in bond_stereodescriptor_colors:
            b_stereos = [
                b_stereo.get_isomers()  # type: ignore
                for color in scc
                for bond in bond_stereodescriptor_colors[color]
                if (b_stereo := smg.get_bond_stereo(bond)).parity is None  # type: ignore
            ]
            perm_generator = itertools.product(*b_stereos)
            smg_gen: Iterator[StereoMolGraph] = unique_generator(
                set_bond_stereos(smg, perm_stereos)
                for perm_stereos in perm_generator
            )

            if not sccs:
                yield from smg_gen
            else:
                yield from perm_gen(
                    smg,
                    last_colors,
                    sccs,
                    unique_stereodescriptors,
                )

    yield from perm_gen(graph, init_colors, sccs, unique_stereodescriptors)


def generate_stereoisomers(
    graph: StereoMolGraph,
    enantiomers: bool = True,
    atoms: Optional[Iterable[AtomId]] = None,
    bonds: Optional[Iterable[Bond]] = None,
) -> Iterator[StereoMolGraph]:
    """Generates all unique stereoisomers of a StereoMolGraph by generation of
    all combinations of parities. Only includes stereocenters which have a
    parity of None. If a parity is set, it is not changed.

    If include_enantiomers is True, both enantiomers of a stereoisomer are
    included, if it is False, only one enantiomer is included.

    Args:
        enantiomers: If True, both enantiomers are included
        atoms: Optional subset of atoms to consider for stereoisomerism
        bonds: Optional subset of bonds to consider for stereoisomerism

    Yields:
        StereoMolGraph: Each unique stereoisomer (and enantiomer if requested)
    """
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

    seen_hash: set[int] = set()
    enantiomers_seen_hash: set[int] = set()

    for a_stereos, b_stereos in itertools.product(
        itertools.product(*atom_stereos), itertools.product(*bond_stereos)
    ):
        stereoisomer = graph.copy()
        for a_stereo in a_stereos:
            stereoisomer.set_atom_stereo(a_stereo)
        for b_stereo in b_stereos:
            stereoisomer.set_bond_stereo(b_stereo)

        stereoisomer_hash = hash(stereoisomer)
        if stereoisomer_hash not in seen_hash:
            seen_hash.add(stereoisomer_hash)
            yield stereoisomer

            if not enantiomers:
                enantiomer = stereoisomer.enantiomer()
                enantiomer_hash = hash(enantiomer)
                if enantiomer_hash not in enantiomers_seen_hash:
                    enantiomers_seen_hash.add(enantiomer_hash)
                    yield enantiomer


def generate_fleeting_stereoisomers(
    graph: StereoCondensedReactionGraph,
    enantiomers: bool = True,
    atoms: Optional[Iterable[AtomId]] = None,
    bonds: Optional[Iterable[Bond]] = None,
) -> Iterator[StereoCondensedReactionGraph]:
    """Generate all unique fleeting stereoisomers of a reaction graph.

    Only includes stereocenters which have a parity of None for the fleeting
    change. If a parity is set, it is not changed.

    Args:
        graph: The reaction graph to generate isomers from
        enantiomers: If True, both enantiomers are included (default: True)
        atoms: Optional subset of atoms to consider for stereoisomerism
        bonds: Optional subset of bonds to consider for stereoisomerism

    Yields:
        StereoCondensedReactionGraph: Each unique fleeting stereoisomer
    """
    # Get atom stereoisomers
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

    # Get bond stereoisomers
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

    seen_isomers_hash = set()
    seen_enantiomers_hash = set()

    for a_stereos, b_stereos in itertools.product(
        itertools.product(*atom_stereos), itertools.product(*bond_stereos)
    ):
        stereoisomer = graph.copy()
        for a_stereo in a_stereos:
            stereoisomer.set_atom_stereo_change(fleeting=a_stereo)
        for b_stereo in b_stereos:
            stereoisomer.set_bond_stereo_change(fleeting=b_stereo)
        
        stereoisomer_hash = hash(stereoisomer)
        if stereoisomer_hash not in seen_isomers_hash:
            seen_isomers_hash.add(stereoisomer_hash)
            yield stereoisomer

            if not enantiomers:
                enantiomer = stereoisomer.enantiomer()
                enantiomer_hash = hash(enantiomer)
                if enantiomer_hash not in seen_enantiomers_hash:
                    seen_enantiomers_hash.add(enantiomer_hash)
                    yield enantiomer


def topological_symmetry_number(graph: StereoMolGraph) -> int:
    """
    Calculated from the number of graph isomorphisms which conserve the
    stereo information.
    symmetry_number = internal_symmetry_number * rotational_symmetry_number
    TODO: add paper reference
    """

    if any(stereo.parity is None for stereo in graph.stereo.values()):
        raise NotImplementedError(
            "all stereocenters have to be defined"
            " to calculate the symmetry number"
        )
    colorings = color_refine_smg(graph)
    mappings = vf2pp_all_isomorphisms(
        graph,
        graph,
        atom_labels=(colorings, colorings),
        stereo=True,
    )
    return deque(enumerate(mappings, 1), maxlen=1)[0][0]
