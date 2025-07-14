"""Implementation of the isomorphism and subgraph isomorphism algorithms.
Based of VF2++ from https://doi.org/10.1016/j.dam.2018.02.018"""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import TYPE_CHECKING, NamedTuple

from stereomolgraph.algorithms.color_refine import color_refine_mg, label_hash

if TYPE_CHECKING:
    from collections.abc import Callable, Hashable, Iterable, Iterator, Mapping
    from typing import Optional, TypeVar

    from stereomolgraph.graphs import (
        AtomId,
        CondensedReactionGraph,
        MolGraph,
        StereoCondensedReactionGraph,
        StereoMolGraph,
    )
    from stereomolgraph.graphs.scrg import StereoChange
    from stereomolgraph.stereodescriptors import Stereo

    KT = TypeVar("KT", bound=Hashable, covariant=True)
    VT = TypeVar("VT", bound=Hashable, covariant=True)

class _Parameters(NamedTuple):
    """
    Parameters of the algorithm.
    :ivar g1_nbrhd: Neighborhood of the first graph
    :ivar g2_nbrhd: Neighborhood of the second graph
    :ivar g1_labels: Labels of the first graph
    :ivar g2_labels: Labels of the second graph
    :ivar nodes_of_g1Labels: Nodes of the first graph grouped by labels
    :ivar nodes_of_g2Labels: Nodes of the second graph grouped by labels
    :ivar g1_degree: Degree of the first graph
    :ivar g2_nodes_of_degree: Nodes of the second graph grouped by degree
    :ivar g1_stereo: Stereochemistry of the first graph
    :ivar g2_stereo: Stereochemistry of the second graph
    :ivar g1_stereo_changes: Stereochemistry changes of the first graph
    :ivar g2_stereo_changes: Stereochemistry changes of the second graph
    """

    # Neighborhood
    g1_nbrhd: Mapping[AtomId, set[AtomId]]
    g2_nbrhd: Mapping[AtomId, set[AtomId]]
    # atomid: label
    g1_labels: Mapping[AtomId, int]
    g2_labels: Mapping[AtomId, int]
    # label: set of atomids
    nodes_of_g1Labels: Mapping[int, set[AtomId]]
    nodes_of_g2Labels: Mapping[int, set[AtomId]]
    # degree: set of atomids
    g1_degree: Mapping[AtomId, int]
    g2_nodes_of_degree: Mapping[int, set[AtomId]]
    # atomid: list of stereos containing the atom
    g1_stereo: Mapping[AtomId, list[Stereo]]
    g2_stereo: Mapping[AtomId, list[Stereo]]

    g1_stereo_changes: Mapping[AtomId, Mapping[StereoChange, list[Stereo]]]
    g2_stereo_changes: Mapping[AtomId, Mapping[StereoChange, list[Stereo]]]


class _State(NamedTuple):
    """
    State of the algorithm.
    :ivar mapping:
    :ivar inverted_mapping:
    :ivar frontier1: neighbors of mapped atoms in g1
    :ivar external1: atoms not in mapping and not in frontier1
    :ivar frontier2: neighbors of mapped atoms in g2
    :ivar external2: atoms not in mapping and not in frontier2
    """

    mapping: dict[int, int]
    inverted_mapping: dict[int, int]

    frontier1: set[int]  # neighbors of mapped atoms in g1
    external1: set[int]  # atoms not in mapping and not in frontier1

    frontier2: set[int]  # neighbors of mapped atoms in g2
    external2: set[int]  # atoms not in mapping and not in frontier2


def group_keys_by_value(many_to_one: dict[KT, VT]) -> dict[VT, set[KT]]:
    """Inverts a many-to-one mapping to create a one-to-many mapping.

    Converts a dictionary where multiple keys may point to the same value
    into a dictionary where each original value maps to a set of all original
    keys.

    >>> group_keys_by_value({"a": 1, "b": 1, "c": 2, "d": 3, "e": 3})
    {1: {'a', 'b'}, 2: {'c'}, 3: {'e', 'd'}}

    :param many_to_one: Dictionary to invert
    :return: Inverted dictionary where each value maps to a set of keys
    """
    inverted: defaultdict[VT, set[KT]] = defaultdict(set)
    for key, value in many_to_one.items():
        inverted[value].add(key)
    return dict(inverted)


def bfs_layers(
    neighbor_dict: Mapping[AtomId, Iterable[AtomId]],
    sources: Iterable[AtomId] | AtomId,
) -> Iterator[list[AtomId]]:
    """
    Generates layers of the graph starting from the source atoms.
    Each layer contains all nodes that are at the same distance from the
    sources.
    The first layer contains the sources.

    :param neighbor_dict: Dictionary of neighbors for each atom
    :param sources: Sources to start from
    """
    if isinstance(sources, int):
        sources = [sources]
    
    current_layer = list(sources)
    visited = set(sources)

    if any(source not in neighbor_dict for source in current_layer):
        raise ValueError("Source atom not in molecule")

    # this is basically BFS, except that the current layer only stores the
    # nodes at same distance from sources at each iteration
    while current_layer:
        yield current_layer
        next_layer: list[AtomId] = []
        for node in current_layer:
            for child in neighbor_dict[node]:
                if child not in visited:
                    visited.add(child)
                    next_layer.append(child)
        current_layer = next_layer


def _sanity_check_and_init(
    g1: StereoMolGraph | MolGraph,
    g2: StereoMolGraph | MolGraph,
    labels: Optional[Iterable[str]] = None,
    color_refine: bool | int = False,
    stereo: bool = False,
    stereo_change: bool = False,
    subgraph: bool = False,
) -> Optional[tuple[_Parameters, _State]]:
    if stereo_change and not stereo:
        raise ValueError("Stereo change is only available for stereo graphs.")
    if subgraph and color_refine:
        raise ValueError(
            "Subgraph isomorphism is not compatible with color refinement."
        )

    g1_nbrhd = g1.neighbors
    g2_nbrhd = g2.neighbors

    if len(g1_nbrhd) == 0 or len(g2_nbrhd) == 0:
        return None

    elif not subgraph:
        if len(g1_nbrhd) != len(g2_nbrhd):
            return None

        if sorted(len(n) for n in g1_nbrhd.values()) != sorted(
            len(n) for n in g2_nbrhd.values()
        ):
            return None

    elif subgraph:
        if len(g1_nbrhd) > len(g2_nbrhd):
            ValueError("The second graph must be larger than the first one.")
        elif len(g1_nbrhd) == len(g2_nbrhd):
            raise ValueError(
                "Both graphs have the same number of atoms. "
                "Do not use subgraph isomorphism in this case."
            )
        elif Counter(len(nbr) for nbr in g1_nbrhd.values()) > Counter(
            len(nbr) for nbr in g2_nbrhd.values()
        ):
            return None

    atom_labels = ("atom_type", *labels) if labels else ("atom_type",)

    max_iter = color_refine
    
    if color_refine: # and not stereo and not stereo_change:
        g1_labels = color_refine_mg(
            g1, atom_labels=atom_labels, max_iter=max_iter
        )
        g2_labels = color_refine_mg(
            g2, atom_labels=atom_labels, max_iter=max_iter
        )

    #elif color_refine and stereo and not stereo_change:
    #
        # TODO: Implement color refinement for stereo graphs
        #g1_labels = color_refine_smg(
        #    g1, atom_labels=atom_labels, max_iter=max_iter
        #)
        #g2_labels = color_refine_smg(
        #    g2, atom_labels=atom_labels, max_iter=max_iter
        #)

    #elif color_refine and stereo and stereo_change:
    #
        #g1_labels = color_refine_scrg(
        #    g1, atom_labels=atom_labels, max_iter=max_iter
        #)
        #g2_labels = color_refine_scrg(
        ##    g2, atom_labels=atom_labels, max_iter=max_iter
        #)

    elif not color_refine:
        g1_labels = label_hash(g1, atom_labels=atom_labels)
        g2_labels = label_hash(g2, atom_labels=atom_labels)

    else:
        raise ValueError("Invalid combination of parameters.")
    
    
    g1_labels_counter = Counter(g1_labels.values())
    g2_labels_counter = Counter(g2_labels.values())

    if not subgraph and g1_labels_counter != g2_labels_counter:
        return None

    elif subgraph and g1_labels_counter > g2_labels_counter:
        return None

    g1_stereo: defaultdict[AtomId, list[Stereo]] = defaultdict(list)
    g2_stereo: defaultdict[AtomId, list[Stereo]] = defaultdict(list)

    if stereo:
        if TYPE_CHECKING:
            assert isinstance(g1, StereoMolGraph)
            assert isinstance(g2, StereoMolGraph)

        for s in g1.stereo.values():
            for atom in s.atoms:
                g1_stereo[atom].append(s)

        for s in g2.stereo.values():
            for atom in s.atoms:
                g2_stereo[atom].append(s)

    g1_stereo_changes: defaultdict[AtomId, defaultdict[StereoChange, list[Stereo]]] = defaultdict(lambda: defaultdict(list)) # type: ignore

    g2_stereo_changes: defaultdict[AtomId, defaultdict[StereoChange, list[Stereo]]] = defaultdict(lambda: defaultdict(list)) # type: ignore

    if stereo_change:
        if TYPE_CHECKING:
            assert isinstance(g1, StereoCondensedReactionGraph)
            assert isinstance(g2, StereoCondensedReactionGraph)
        
        for _, stereo_change_dict in g1.atom_stereo_changes.items():
            for stereo_change_enum, atom_stereo in stereo_change_dict.items():
                if atom_stereo is not None:
                    for atom in atom_stereo.atoms:
                        g1_stereo_changes[atom][stereo_change_enum].append(
                            atom_stereo
                        )

        for _, stereo_change_dict in g2.atom_stereo_changes.items():
            for stereo_change_enum, atom_stereo in stereo_change_dict.items():
                if atom_stereo is not None:
                    for atom in atom_stereo.atoms:
                        g2_stereo_changes[atom][stereo_change_enum].append(
                            atom_stereo
                        )

        for _, stereo_change_dict in g1.bond_stereo_changes.items():
            for stereo_change_enum, bond_stereo in stereo_change_dict.items():
                if bond_stereo is not None:
                    for atom in bond_stereo.atoms:
                        g1_stereo_changes[atom][stereo_change_enum].append(
                            bond_stereo
                        )

        for _, stereo_change_dict in g2.bond_stereo_changes.items():
            for stereo_change_enum, bond_stereo in stereo_change_dict.items():
                if bond_stereo is not None:
                    for atom in bond_stereo.atoms:
                        g2_stereo_changes[atom][stereo_change_enum].append(
                            bond_stereo
                        )

    g1_degree = {a: len(n_set) for a, n_set in g1_nbrhd.items()}
    g2_degree = {a: len(n_set) for a, n_set in g2_nbrhd.items()}

    params = _Parameters(
        g1_nbrhd,
        g2_nbrhd,
        g1_labels,
        g2_labels,
        group_keys_by_value(g1_labels),
        group_keys_by_value(g2_labels),
        g1_degree,
        group_keys_by_value(g2_degree),
        g1_stereo,
        g2_stereo,
        g1_stereo_changes,
        g2_stereo_changes,
    )

    state = _State({}, {}, set(), set(g1_nbrhd), set(), set(g2_nbrhd))

    return params, state


def _wrap_all(
    *funcs: Callable[[AtomId, AtomId, _State, _Parameters], bool],
) -> Callable[[AtomId, AtomId, _State, _Parameters], bool]:
    def wrapper(a: AtomId, b: AtomId, state: _State, params: _Parameters) -> bool:
        return all(f(a, b, state, params) for f in funcs)
    
    return wrapper


def vf2pp_all_isomorphisms(
    g1: MolGraph
    | StereoMolGraph
    | CondensedReactionGraph
    | StereoCondensedReactionGraph,
    g2: MolGraph
    | StereoMolGraph
    | CondensedReactionGraph
    | StereoCondensedReactionGraph,
    labels: Optional[Iterable[str]] = None,
    bond_labels: Optional[Iterable[str]] = None,
    color_refine: bool | int = True,
    stereo: bool = False,
    stereo_change: bool = False,
    subgraph: bool = False,
) -> Iterator[dict[AtomId, AtomId]]:
    if params_state := _sanity_check_and_init(
        g1, g2, labels, color_refine, stereo, stereo_change, subgraph
    ):
        params, state = params_state
    else:
        return  # if no isomorphisms return like an empty generator

    # setup helper function based on input parameters

    feasibility_funcs: list[Callable[[AtomId, AtomId, _State, _Parameters], bool]] = []
    if subgraph:
        feasibility_funcs.append(_subgraph_feasibility)
        if stereo:
            feasibility_funcs.append(_subgraph_stereo_feasibility)
        if stereo_change:
            feasibility_funcs.append(_subgraph_stereo_change_feasibility)
            
    elif not subgraph:
        feasibility_funcs.append(_graph_feasibility)
        if stereo:
            feasibility_funcs.append(_stereo_feasibility)
        if stereo_change:
            feasibility_funcs.append(_stereo_change_feasibility)
    else:
        raise ValueError("Invalid combination of parameters.")

    feasibility = _wrap_all(*feasibility_funcs) # type: ignore
    revert_state = _revert_state
    update_state = _update_state
    find_candidates = (_find_subgraph_candidates if subgraph
                       else _find_candidates)

    # to avoid overhead
    mapping = state.mapping
    inverted_mapping = state.inverted_mapping
    termination_length = len(g1)

    # Initialize the stack
    node_order: list[AtomId] = _matching_order(params)
    candidates: set[AtomId] = find_candidates(node_order[0], state, params)

    stack: list[tuple[AtomId, set[AtomId]]] = []
    stack.append((node_order[0], candidates))

    # Index of the node from the order, currently being examined
    matching_atom_index = 1

    while stack:
        matching_atom, candidates = stack[-1]
        candidate = candidates.pop() if candidates else None
        if candidate is None:
            # If no remaining candidates, return to a previous state,
            # and follow another branch
            stack.pop()
            matching_atom_index -= 1
            if stack:
                # Pop the previously added u-v pair, and look for
                # a different candidate _v for u
                last_atom1, _ = stack[-1]
                last_atom2 = mapping.pop(last_atom1)
                inverted_mapping.pop(last_atom2)
                revert_state(last_atom1, last_atom2, state, params)
            continue

        mapping[matching_atom] = candidate
        inverted_mapping[candidate] = matching_atom

        if feasibility(matching_atom, candidate, state, params):
            if len(mapping) == termination_length:
                yield mapping.copy()
                mapping.pop(matching_atom)
                inverted_mapping.pop(candidate)
                continue

            update_state(matching_atom, candidate, state, params)
            # Append the next node and its candidates to the stack
            matching_atom = node_order[matching_atom_index]
            candidates = find_candidates(matching_atom, state, params)
            stack.append((matching_atom, candidates))
            matching_atom_index += 1

        else:  # if not feaseble
            mapping.pop(matching_atom)
            inverted_mapping.pop(candidate)


def _graph_feasibility(
    u: AtomId, v: AtomId, state: _State, params: _Parameters
):
    g1_nbrhd, g2_nbrhd, g1_labels, g2_labels, *_ = params # type: ignore
    _, _, frontier1, external1, frontier2, external2 = state

    t1_labels: list[int] = []
    t2_labels: list[int] = []
    t1_tilde_labels: list[int] = []
    t2_tilde_labels: list[int] = []

    for n in g1_nbrhd[u]:
        if n in external1:
            t1_tilde_labels.append(g1_labels[n])
        elif n in frontier1:
            t1_labels.append(g1_labels[n])

    for n in g2_nbrhd[v]:
        if n in external2:
            t2_tilde_labels.append(g2_labels[n])
        elif n in frontier2:
            t2_labels.append(g2_labels[n])

    t1_labels.sort()
    t2_labels.sort()

    if t1_labels != t2_labels:
        return False

    t1_tilde_labels.sort()
    t2_tilde_labels.sort()

    if t1_tilde_labels != t2_tilde_labels:
        return False

    return True

def _subgraph_feasibility(
    u: AtomId, v: AtomId, state: _State, params: _Parameters
) -> bool:
    g1_nbrhd, g2_nbrhd, g1_labels, g2_labels, *_ = params # type: ignore
    _, _, frontier1, external1, frontier2, external2 = state

    counter1 = Counter(g1_labels[n] for n in g1_nbrhd[u] if n in frontier1)
    counter2 = Counter(g2_labels[n] for n in g2_nbrhd[v] if n in frontier2)
    if counter1 > counter2:
        return False

    counter1 = Counter(g1_labels[n] for n in g1_nbrhd[u] if n in external1)
    counter2 = Counter(g2_labels[n] for n in g2_nbrhd[v] if n in external2)
    return counter1 <= counter2

def _stereo_feasibility(
    u: AtomId, v: AtomId, state: _State, params: _Parameters
) -> bool:
    s1 = [
        stereo.__class__(
            atoms=tuple([state.mapping[a] for a in stereo.atoms]),
            parity=stereo.parity,
        )
        for stereo in params.g1_stereo[u]
        if all([a in state.mapping for a in stereo.atoms])
    ]

    s2 = [
        stereo
        for stereo in params.g2_stereo[v]
        if all([a in state.inverted_mapping for a in stereo.atoms])
    ]

    if len(s2) != len(s1):
        return False

    if all(s in s2 for s in s1):
        return True
    return False

def _subgraph_stereo_feasibility(
    u: AtomId, v: AtomId, state: _State, params: _Parameters
) -> bool: ...

def _stereo_change_feasibility(
    u: AtomId, v: AtomId, state: _State, params: _Parameters
) -> bool:
    s1 = {
        (
            stereo_change,
            stereo.__class__(
                atoms=tuple([state.mapping[a] for a in stereo.atoms]),
                parity=stereo.parity,
            ),
        )
        for stereo_change, stereo_list in params.g1_stereo_changes[u].items()
        for stereo in stereo_list
        if stereo is not None # type: ignore
        and all([a in state.mapping for a in stereo.atoms])
    }

    s2 = {
        (stereo_change, stereo)
        for stereo_change, stereo_list in params.g2_stereo_changes[u].items()
        for stereo in stereo_list
        if stereo is not None # type: ignore
        and all([a in state.inverted_mapping for a in stereo.atoms])
    }

    if s1 == s2:
        return True
    return False

def _subgraph_stereo_change_feasibility(
    u: AtomId, v: AtomId, state: _State, params: _Parameters
) -> bool: ...
    # TODO: Implement subgraph stereo change feasibility check

def _matching_order(params: _Parameters) -> list[AtomId]:
    g1_nbrhd, _, g1_labels, _, _, nodes_of_g2Labels, *_ = params

    V1_unordered = set(g1_labels.keys())
    label_rarity = {
        label: len(nodes) for label, nodes in nodes_of_g2Labels.items()
    }
    used_degrees = {node: 0 for node in g1_nbrhd}
    node_order: list[int] = []

    while V1_unordered:
        max_rarity = min(label_rarity[g1_labels[x]] for x in V1_unordered)
        rarest_nodes = [
            n for n in V1_unordered if label_rarity[g1_labels[n]] == max_rarity
        ]
        max_node: AtomId = max(rarest_nodes, key={a: len(n_set) for a, n_set in g1_nbrhd.items()}.get,) # type: ignore
        assert isinstance(max_node, int)
        for dlevel_nodes in bfs_layers(g1_nbrhd, max_node):
            nodes_to_add = dlevel_nodes.copy()
            while nodes_to_add:
                max_used_degree = max(used_degrees[n] for n in nodes_to_add)
                max_used_degree_nodes = [
                    n
                    for n in nodes_to_add
                    if used_degrees[n] == max_used_degree
                ]
                max_degree = max(
                    len(g1_nbrhd[n]) for n in max_used_degree_nodes
                )
                max_degree_nodes = [
                    n
                    for n in max_used_degree_nodes
                    if len(g1_nbrhd[n]) == max_degree
                ]
                next_node = min(
                    max_degree_nodes, key=lambda x: label_rarity[g1_labels[x]]
                )

                node_order.append(next_node)
                for node in g1_nbrhd[next_node]:
                    used_degrees[node] += 1

                nodes_to_add.remove(next_node)
                label_rarity[g1_labels[next_node]] -= 1
                V1_unordered.discard(next_node)

    return node_order


def _find_candidates(
    u: AtomId, state: _State, params: _Parameters
) -> set[AtomId]:
    (
        g1_nbrhd,
        g2_nbrhd,
        g1_labels,
        _,
        _,
        nodes_of_g2Labels,
        g1_deg,
        g2_a_of_deg,
        *_,
    ) = params
    mapping, inverted_mapping, *_, external2 = state

    covered_nbrs = [nbr for nbr in g1_nbrhd[u] if nbr in mapping]

    if not covered_nbrs:
        candidates = set(nodes_of_g2Labels[g1_labels[u]])
        candidates.intersection_update(external2, g2_a_of_deg[g1_deg[u]])
        candidates.difference_update(inverted_mapping)

        return candidates

    nbr1 = covered_nbrs[0]
    candidates = set(g2_nbrhd[mapping[nbr1]])

    for nbr1 in covered_nbrs[1:]:
        candidates.intersection_update(g2_nbrhd[mapping[nbr1]])

    candidates.difference_update(inverted_mapping)
    candidates.intersection_update(
        nodes_of_g2Labels[g1_labels[u]], g2_a_of_deg[g1_deg[u]]
    )
    return candidates


def _find_subgraph_candidates(
    u: AtomId, state: _State, params: _Parameters
) -> set[AtomId]:
    (
        g1_nbrhd,
        g2_nbrhd,
        g1_labels,
        _,
        _,
        nodes_of_g2Labels,
        *_,
    ) = params
    mapping, inverted_mapping, *_, external2 = state

    covered_nbrs = [nbr for nbr in g1_nbrhd[u] if nbr in mapping]

    if not covered_nbrs:
        candidates = set(nodes_of_g2Labels[g1_labels[u]])
        candidates.intersection_update(external2)
        candidates.difference_update(inverted_mapping)
        # candidates.intersection_update(g2_a_of_deg[degree] for degree
        #                               in g2_a_of_deg if degree < g1_deg[u])
        return candidates

    nbr1 = covered_nbrs[0]
    common_nodes = set(g2_nbrhd[mapping[nbr1]])

    for nbr1 in covered_nbrs[1:]:
        common_nodes.intersection_update(g2_nbrhd[mapping[nbr1]])

    common_nodes.difference_update(inverted_mapping)
    common_nodes.intersection_update(nodes_of_g2Labels[g1_labels[u]])

    # common_nodes.intersection_update(g2_a_of_deg[degree] for degree
    #                                 in g2_a_of_deg if degree < g1_deg[u])

    return common_nodes


def _update_state(
    new_atom1: AtomId, new_atom2: AtomId, state: _State, params: _Parameters
) -> None:
    g1_nbrhd, g2_nbrhd, *_ = params
    mapping, inverted_mapping, frontier1, external1, frontier2, external2 = (
        state
    )

    unmapped_neighbors1 = {n for n in g1_nbrhd[new_atom1] if n not in mapping}
    frontier1 |= unmapped_neighbors1
    external1 -= unmapped_neighbors1
    frontier1.discard(new_atom1)
    external1.discard(new_atom1)

    unmapped_neighbors2 = {
        n for n in g2_nbrhd[new_atom2] if n not in inverted_mapping
    }
    frontier2 |= unmapped_neighbors2
    external2 -= unmapped_neighbors2
    frontier2.discard(new_atom2)
    external2.discard(new_atom2)

    return


def _revert_state(
    last_atom1: AtomId, last_atom2: AtomId, state: _State, params: _Parameters
) -> None:
    # If the node we want to remove from the mapping, has at least one
    # covered neighbor, add it to frontier1.
    g1_nbrhd, g2_nbrhd, *_ = params
    mapping, inverted_mapping, frontier1, external1, frontier2, external2 = (
        state
    )

    has_covered_neighbor = False
    for neighbor in g1_nbrhd[last_atom1]:
        if neighbor in mapping:
            # if a neighbor of the excluded node1 is in the mapping,
            # keep node1 in frontier1
            has_covered_neighbor = True
            frontier1.add(last_atom1)
        else:
            # check if its neighbor has another connection with a covered node.
            # If not, only then exclude it from frontier1
            if any(nbr in mapping for nbr in g1_nbrhd[neighbor]):
                continue
            frontier1.discard(neighbor)
            external1.add(neighbor)

    # Case where the node is not present in neither the mapping nor frontier1.
    # By definition, it should belong to external1
    if not has_covered_neighbor:
        external1.add(last_atom1)

    has_covered_neighbor = False
    for neighbor in g2_nbrhd[last_atom2]:
        if neighbor in inverted_mapping:
            has_covered_neighbor = True
            frontier2.add(last_atom2)
        else:
            if any(nbr in inverted_mapping for nbr in g2_nbrhd[neighbor]):
                continue
            frontier2.discard(neighbor)
            external2.add(neighbor)

    if not has_covered_neighbor:
        external2.add(last_atom2)

