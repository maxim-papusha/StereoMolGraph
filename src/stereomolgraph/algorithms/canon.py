"""Frontier canonization with lexicographic and automorphism pruning.

The label sequence plays the role of McKay's path indicator (Section 2.26),
and equal numbered-graph certificates reveal automorphisms (Section 2.21):
https://users.cecs.anu.edu.au/~bdm/papers/pgi.pdf

This is a specialized frontier search: colors are refined once initially,
rather than individualizing and refining an ordered partition at every node.
"""

from collections import Counter, defaultdict
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from stereomolgraph import (
    AtomId,
    Bond,
    CondensedReactionGraph,
    MolGraph,
    StereoCondensedReactionGraph,
    StereoMolGraph,
)
from stereomolgraph.algorithms.circular import (
    color_refine_crg,
    color_refine_mg,
    color_refine_scrg,
    color_refine_smg,
)
from stereomolgraph.graphs.crg import Change
from stereomolgraph.stereodescriptors import Stereo


__all__ = ["canon_atom_num"]

_Label = tuple[Any, ...]
_Certificate = tuple[Any, ...]


@dataclass
class _Parameters:
    """Fixed graph data and precomputed labels shared by all search branches."""

    nbrhds: Mapping[AtomId, frozenset[AtomId]]
    atom_color_frequency: Mapping[AtomId, int]
    degrees: Mapping[AtomId, int]
    colors: Mapping[AtomId, int]
    atom_types: Mapping[AtomId, int]
    active_atoms: frozenset[AtomId]
    stereo_change_atoms: frozenset[AtomId]
    bond_changes: Mapping[Bond, str]

    stereo_of_atoms: Mapping[AtomId, list[Stereo]]
    stereo: Mapping[AtomId | Bond, Stereo]

    stereo_changes: list[tuple[AtomId | Bond, Change, Stereo]]
    stereo_changes_of_atoms: Mapping[AtomId, list[tuple[Change, Stereo]]]


@dataclass
class _CanonState:
    """Mutable search state, with one queue and label entry per search level.

    Dict insertion order lets popitem() undo the last atom assignment.
    """

    best_labels: list[_Label] = field(default_factory=list)
    best_mapping: dict[AtomId, int] = field(default_factory=dict)

    current_labels: list[_Label] = field(default_factory=list)
    current_mapping: dict[AtomId, int] = field(default_factory=dict)

    frontier: set[AtomId] = field(default_factory=set)

    queue: list[set[AtomId]] = field(default_factory=list)
    explored: list[set[AtomId]] = field(default_factory=list)

    best_certificate: _Certificate = ()
    automorphisms: list[dict[AtomId, AtomId]] = field(default_factory=list)
    prune_automorphisms: bool = True


def _stereo_record(
    site: AtomId | Bond, stereo: Stereo, mapping: Mapping[AtomId, int]
) -> tuple[Any, ...]:
    site_label = (
        ("atom", (mapping[site],))
        if not isinstance(site, frozenset)
        else ("bond", tuple(sorted(mapping[atom] for atom in site)))
    )
    descriptor = stereo.__class__(
        atoms=tuple(0 if atom is None else mapping[atom] for atom in stereo.atoms),
        parity=stereo.parity,
    )
    return (
        site_label,
        stereo.__class__.__module__,
        stereo.__class__.__qualname__,
        descriptor.canonical_form(),
    )


def _graph_certificate(
    params: _Parameters, mapping: Mapping[AtomId, int]
) -> _Certificate:
    """Encode a complete numbered graph without relying on color hashes alone.

    Stereo records include their attachment site and descriptor class. None is
    encoded as 0, outside the assigned numbers 1..n. Stereo.__eq__ is unsuitable
    here because unspecified parity matches defined parity there.
    """
    atoms = tuple(
        sorted(
            (
                mapping[atom],
                params.atom_types[atom],
                params.colors[atom],
                tuple(
                    sorted(
                        (
                            mapping[nbr],
                            params.bond_changes.get(frozenset((atom, nbr)), ""),
                        )
                        for nbr in params.nbrhds[atom]
                    )
                ),
            )
            for atom in params.nbrhds
        )
    )
    stereo_records = tuple(
        sorted(
            _stereo_record(site, stereo, mapping)
            for site, stereo in params.stereo.items()
        )
    )
    stereo_changes = tuple(
        sorted(
            (change.value, _stereo_record(site, stereo, mapping))
            for site, change, stereo in params.stereo_changes
        )
    )
    return atoms, stereo_records, stereo_changes


def _record_leaf(params: _Parameters, state: _CanonState) -> None:
    certificate = _graph_certificate(params, state.current_mapping)
    if state.prune_automorphisms and certificate == state.best_certificate:
        inverse = {number: atom for atom, number in state.current_mapping.items()}
        # Equal exact certificates prove that the two numberings differ by an
        # automorphism: current_mapping^-1 composed with best_mapping.
        permutation = {
            atom: inverse[number] for atom, number in state.best_mapping.items()
        }
        if any(a != b for a, b in permutation.items()) and (
            permutation not in state.automorphisms
        ):
            state.automorphisms.append(permutation)

    # Exact certificates break any ties left by the hashed traversal labels.
    if not state.best_labels or (state.current_labels, certificate) < (
        state.best_labels,
        state.best_certificate,
    ):
        state.best_labels = state.current_labels.copy()
        state.best_certificate = certificate
        state.best_mapping = state.current_mapping.copy()


def _prune_orbits(state: _CanonState) -> None:
    """Discard siblings equivalent to an explored sibling under proven symmetry.

    Only generators fixing the numbered prefix pointwise may be used. They
    generate a subgroup of its stabilizer, so missed symmetries cost time only.
    Recompute after each child because that subtree may reveal new generators.
    """
    if not state.automorphisms or not state.explored[-1]:
        return
    generators = [
        permutation
        for permutation in state.automorphisms
        if all(permutation[atom] == atom for atom in state.current_mapping)
    ]
    orbit = state.explored[-1].copy()
    pending = list(orbit)
    while pending:
        atom = pending.pop()
        for permutation in generators:
            image = permutation[atom]
            if image not in orbit:
                orbit.add(image)
                pending.append(image)
    pruned = state.queue[-1] & orbit
    state.queue[-1] -= pruned


def _canonical_stereo_label(
    atom_id: AtomId,
    stereo: Stereo,
    current_mapping: Mapping[AtomId, int],
    colors: Mapping[AtomId, int],
) -> tuple[Any, ...]:
    """Order markers: self=0, assigned=1..n, unassigned=n+1, lone pair=n+2."""
    unmapped_atom_num = len(colors) + 1
    atoms = []
    for atom in stereo.atoms:
        if atom == atom_id:
            number = 0
        elif atom is None:
            number = unmapped_atom_num + 1
        else:
            number = current_mapping.get(atom, unmapped_atom_num)
        atoms.append((number, colors.get(atom, -1)))
    return stereo.__class__(
        atoms=tuple(atoms),
        parity=stereo.parity,
    ).canonical_form()


def _candidate_label(
    atom_id: AtomId,
    params: _Parameters,
    current_mapping: Mapping[AtomId, int],
) -> _Label:
    """Prefer reaction groups, non-H atoms, then changed-bond connections.

    Rare colors and high degree resolve further ties. Neighbor records put
    reaction-tagged bonds first, then earlier assigned numbers. Changed stereo
    descriptors are compared before ordinary stereo and the final color.
    """
    unmapped_atom_num = len(params.nbrhds) + 1
    if atom_id in params.active_atoms:
        reaction_group = 0
    elif atom_id in params.stereo_change_atoms:
        reaction_group = 1
    else:
        reaction_group = 2

    neighbor_label = []
    changed_connection = False
    for nbr in params.nbrhds[atom_id]:
        change = params.bond_changes.get(frozenset((atom_id, nbr)), "")
        neighbor_label.append(
            (
                not change,
                current_mapping.get(nbr, unmapped_atom_num),
                params.colors[nbr],
                change,
            )
        )
        if change and nbr in current_mapping:
            changed_connection = True
    neighbor_label.sort()
    stereo_change_label = sorted(
        (
            change.value,
            stereo.__class__.__module__,
            stereo.__class__.__qualname__,
            _canonical_stereo_label(atom_id, stereo, current_mapping, params.colors),
        )
        for change, stereo in params.stereo_changes_of_atoms.get(atom_id, [])
    )
    stereo_label = sorted(
        _canonical_stereo_label(
            atom_id,
            stereo,
            current_mapping=current_mapping,
            colors=params.colors,
        )
        for stereo in params.stereo_of_atoms.get(atom_id, [])
    )

    return (
        reaction_group,
        params.atom_types[atom_id] == 1,
        not changed_connection,
        params.atom_color_frequency[atom_id],
        -params.degrees[atom_id],
        neighbor_label,
        stereo_change_label,
        stereo_label,
        params.colors[atom_id],
    )


def _update_queue(params: _Parameters, state: _CanonState) -> None:
    """Push a search level containing atoms tied for the best next-step label.

    An empty candidate set marks a complete or pruned branch.
    """

    frontier = state.frontier or (params.nbrhds.keys() - state.current_mapping.keys())
    for group in (params.active_atoms, params.stereo_change_atoms):
        remaining = group - state.current_mapping.keys()
        if remaining:
            # Restart within this group rather than crossing a lower-priority atom.
            frontier = frontier & remaining or remaining
            break

    candidate_labels = {
        atom_id: _candidate_label(
            atom_id,
            params=params,
            current_mapping=state.current_mapping,
        )
        for atom_id in frontier
    }

    atoms = set()
    best_step_label = ()

    if candidate_labels:
        best_step_label = min(candidate_labels.values())
        potential_new_labels = state.current_labels + [best_step_label]

        if (
            not state.best_labels
            or potential_new_labels <= state.best_labels[: len(potential_new_labels)]
        ):
            atoms = {
                atom_id
                for atom_id, candidate_label in candidate_labels.items()
                if candidate_label == best_step_label
            }

    state.queue.append(atoms)
    state.explored.append(set())
    state.current_labels.append(best_step_label)

    if not atoms and len(params.nbrhds) == len(state.current_mapping):
        _record_leaf(params, state)


def _initialize(g: MolGraph) -> _Parameters:
    """Prepare fixed colors, neighborhoods, and stereo records for the search."""

    if isinstance(g, StereoCondensedReactionGraph):
        refine = color_refine_scrg
    elif isinstance(g, CondensedReactionGraph):
        refine = color_refine_crg
    elif isinstance(g, StereoMolGraph):
        refine = color_refine_smg
    else:
        refine = color_refine_mg

    colors = {}
    if g.atoms:
        colors = {atom: int(label) for atom, label in zip(g.atoms, refine(g))}

    degrees = {atom: len(nbrs) for atom, nbrs in g.neighbors.items()}

    label_class_sizes = Counter(colors.values())
    atom_color_frequency = {atom: label_class_sizes[colors[atom]] for atom in g.atoms}

    stereo_of_atoms: dict[AtomId, list[Stereo]] = defaultdict(list)
    if isinstance(g, StereoMolGraph):
        for stereo in g.stereo.values():
            for atom in stereo.atoms:
                if atom is not None:
                    stereo_of_atoms[atom].append(stereo)

    stereo_changes = []
    stereo_changes_of_atoms: dict[AtomId, list[tuple[Change, Stereo]]] = defaultdict(
        list
    )
    if isinstance(g, StereoCondensedReactionGraph):
        for changes in (g.atom_stereo_changes, g.bond_stereo_changes):
            for site, records in changes.items():
                for change, stereo in records.items():
                    if stereo is not None:
                        stereo_changes.append((site, change, stereo))
                        for atom in {atom for atom in stereo.atoms if atom is not None}:
                            stereo_changes_of_atoms[atom].append((change, stereo))

    bond_changes = {}
    active_atoms: set[AtomId] = set()
    if isinstance(g, CondensedReactionGraph):
        for bond, attributes in g.bonds_with_attributes.items():
            change = attributes.get("reaction")
            if change is not None:
                bond_changes[bond] = change.value
            if change in (Change.FORMED, Change.BROKEN):
                active_atoms.update(bond)

    return _Parameters(
        nbrhds=g.neighbors,
        colors=colors,
        degrees=degrees,
        atom_color_frequency=atom_color_frequency,
        stereo_of_atoms=stereo_of_atoms,
        atom_types=dict(zip(g.atoms, g.atom_types)),
        active_atoms=frozenset(active_atoms),
        stereo_change_atoms=frozenset(stereo_changes_of_atoms),
        bond_changes=bond_changes,
        stereo=g.stereo if isinstance(g, StereoMolGraph) else {},
        stereo_changes=stereo_changes,
        stereo_changes_of_atoms=stereo_changes_of_atoms,
    )


def _search(
    params: _Parameters,
    *,
    prune_automorphisms: bool = True,
) -> _CanonState:
    """Search for the smallest numbering under fixed graph parameters."""

    state = _CanonState(
        prune_automorphisms=prune_automorphisms,
    )
    _update_queue(params, state)

    while state.queue:
        if state.prune_automorphisms:
            _prune_orbits(state)
        # Extend the current numbering.
        if state.queue[-1]:
            next_atom = state.queue[-1].pop()
            state.explored[-1].add(next_atom)

            state.current_mapping[next_atom] = len(state.current_mapping) + 1

            state.frontier.discard(next_atom)
            state.frontier |= params.nbrhds[next_atom] - state.current_mapping.keys()

            _update_queue(params, state)

        # Backtrack once this level is exhausted.
        else:
            state.queue.pop()
            state.explored.pop()
            state.current_labels.pop()
            if state.queue:
                atom_id, _ = state.current_mapping.popitem()

                not_frontier_anymore = {
                    nbr
                    for nbr in params.nbrhds[atom_id]
                    if state.current_mapping.keys().isdisjoint(params.nbrhds[nbr])
                }
                state.frontier -= not_frontier_anymore
                if not state.current_mapping.keys().isdisjoint(params.nbrhds[atom_id]):
                    state.frontier.add(atom_id)

    return state


def canon_atom_num(
    g: MolGraph,
    *,
    prune_automorphisms: bool = True,
) -> Mapping[AtomId, int]:
    """Assign canonical numbers from 1 in a single search of the whole graph.

    For reaction graphs, number all atoms incident to formed or broken bonds
    first, then remaining atoms included in stereo-change descriptors, then the
    rest. These priorities hold across disconnected or separated regions;
    lone-pair placeholders are excluded.
    Within each group, prefer non-hydrogens in the frontier, then atoms linked
    to numbered atoms by reaction-tagged bonds. Stereo-change descriptors and
    their change types contribute directly to traversal labels. Ordinary molecular
    graphs finish the current component before choosing another start atom.
    Worse label prefixes and, by default, automorphism-equivalent siblings are
    pruned.

    Uses atom types, connectivity, stored stereo, and reaction bond/stereo changes.
    Other atom/bond attributes are not encoded.
    """
    return _search(_initialize(g), prune_automorphisms=prune_automorphisms).best_mapping
