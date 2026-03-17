from collections import Counter, defaultdict
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from stereomolgraph import AtomId, MolGraph, StereoMolGraph
from stereomolgraph.algorithms.color_refine import (
    color_refine_mg,
    color_refine_smg,
)
from stereomolgraph.stereodescriptors import Stereo

CanonNum = int

LabelClassSize = int
Degree = int
Color = int
Label = tuple[Any, ...]


@dataclass
class _Parameters:
    """Parameters that are invariant across the search.

    :ivar nbrhds: The neighborhood of each atom in the graph.
    """

    nbrhds: Mapping[AtomId, set[AtomId]]
    atom_color_frequency: Mapping[AtomId, LabelClassSize]
    degrees: Mapping[AtomId, Degree]
    colors: Mapping[AtomId, Color]

    stereo_of_atoms: Mapping[AtomId, list[Stereo]]

    # TODO: stereo, bond_change, stereo_change, ...


@dataclass
class _CanonState:
    """State of the search that is updated at each step.
    Everywhere where the order of elements should not matter sets are used instead of
    lists. Python dicts keep the insertion order and are FILO"""

    best_labels: list[Label] = field(default_factory=list)
    best_mapping: dict[AtomId, CanonNum] = field(default_factory=dict)

    current_labels: list[Label] = field(default_factory=list)
    current_mapping: dict[AtomId, CanonNum] = field(default_factory=dict)

    frontier: set[AtomId] = field(default_factory=set)

    queue: list[set[AtomId]] = field(default_factory=list)


def better_than(label_like1, label_like2):
    """Tuples and lists are compared element wise lexicographically.
    If all shared elements are equal, the shorter one is considered smaller."""
    return label_like1 > label_like2


def _canonical_stereo_label(
    atom_id: AtomId,
    stereo: Stereo,
    current_mapping: Mapping[AtomId, CanonNum],
    colors: Mapping[AtomId, Color],
) -> tuple[Any, ...]:
    unmapped_atom_num = -1
    self_marker = len(colors) + 1  # n_atoms + 1, always > any assigned canon num
    return stereo.__class__(
        atoms=tuple(
            (
                self_marker
                if a == atom_id
                else current_mapping.get(a, unmapped_atom_num),
                colors.get(a, unmapped_atom_num),
            )
            for a in stereo.atoms
        ),
        parity=stereo.parity,
    ).canonical_form()


def _candidate_label(
    atom_id: AtomId,
    params: _Parameters,
    current_mapping: Mapping[AtomId, CanonNum],
) -> Label:
    unmapped_atom_num: CanonNum = -1

    neighbor_label = sorted(
        (
            (current_mapping.get(nbr, unmapped_atom_num), params.colors[nbr])
            for nbr in params.nbrhds[atom_id]
        ),
        reverse=True,
    )
    stereo_label = sorted(
        (
            _canonical_stereo_label(
                atom_id,
                stereo,
                current_mapping=current_mapping,
                colors=params.colors,
            )
            for stereo in params.stereo_of_atoms.get(atom_id, [])
        ),
        reverse=True,
    )

    return (
        -1 * params.atom_color_frequency[atom_id],
        params.degrees[atom_id],
        neighbor_label,
        stereo_label,
        params.colors[atom_id],
    )


def update_queue(params: _Parameters, state: _CanonState) -> None:
    """Return frontier atoms tied for the best next-step label.

    An empty result means the current branch is complete or can be pruned.
    """

    frontier = (
        state.frontier
        if state.frontier
        else set(params.nbrhds) - state.current_mapping.keys()
    )

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
        best_step_label = max(candidate_labels.values())
        potential_new_labels = state.current_labels + [best_step_label]

        if not better_than(
            state.best_labels[: len(potential_new_labels)],
            potential_new_labels,
        ):
            atoms = {
                atom_id
                for atom_id, candidate_label in candidate_labels.items()
                if candidate_label == best_step_label
            }

    state.queue.append(atoms)
    state.current_labels.append(best_step_label)

    if (
        not atoms
        and len(params.nbrhds) == len(state.current_mapping)
        and better_than(state.current_labels, state.best_labels)
    ):
        state.best_labels = state.current_labels.copy()
        state.best_mapping = state.current_mapping.copy()


def initialize(g: MolGraph) -> _Parameters:
    """Prepare immutable inputs and initial state for non-stereo canonical
    enumeration."""
    nbrhd = g.neighbors

    if isinstance(g, StereoMolGraph):
        colors = {atom: int(label) for atom, label in zip(g.atoms, color_refine_smg(g))}
    elif isinstance(g, MolGraph):
        colors = {atom: int(label) for atom, label in zip(g.atoms, color_refine_mg(g))}

    degrees = {atom: len(nbrs) for atom, nbrs in nbrhd.items()}

    label_class_sizes = Counter(colors.values())
    atom_color_frequency = {atom: label_class_sizes[colors[atom]] for atom in g.atoms}

    stereo_of_atoms: dict[AtomId, list[Stereo]] = defaultdict(list)
    if isinstance(g, StereoMolGraph):
        for s in g.stereo.values():
            for atom in s.atoms:
                if atom is not None:
                    stereo_of_atoms[atom].append(s)

    param = _Parameters(
        nbrhds=nbrhd,
        colors=colors,
        degrees=degrees,
        atom_color_frequency=atom_color_frequency,
        stereo_of_atoms=stereo_of_atoms,
    )

    return param


def canon_atom_num(g: MolGraph) -> Mapping[AtomId, CanonNum]:
    """Return canonical atom order using best-first branching and tie backtracking."""

    canon_num: CanonNum = len(g.atoms)

    params = initialize(g)
    state = _CanonState()
    update_queue(params, state)

    while state.queue:
        # execute step forward
        if state.queue[-1]:
            # Enumeration start with "n_atoms" and counts down to 1,
            # so that the first mapped atom gets the highest canon number.
            next_atom = state.queue[-1].pop()

            state.current_mapping[next_atom] = canon_num
            canon_num -= 1

            # incremental frontier update
            state.frontier.discard(next_atom)
            state.frontier |= params.nbrhds[next_atom] - state.current_mapping.keys()

            # plans next step
            update_queue(params, state)

        # execute step back
        else:
            _empty_set = state.queue.pop()

            _last_label = state.current_labels.pop()
            if state.queue:
                atom_id, _ = state.current_mapping.popitem()
                canon_num += 1

                not_frontier_anymore = {
                    nbr
                    for nbr in params.nbrhds[atom_id]
                    if not any(
                        nbr2 in state.current_mapping for nbr2 in params.nbrhds[nbr]
                    )
                }
                state.frontier -= not_frontier_anymore
                if params.nbrhds[atom_id].intersection(state.current_mapping.keys()):
                    state.frontier.add(atom_id)

    return state.best_mapping
