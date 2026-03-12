from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from stereomolgraph import AtomId, MolGraph
from stereomolgraph.algorithms.color_refine import color_refine_mg


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

    # TODO: stereo, bond_change, stereo_change, ...


@dataclass
class _CanonState:
    best_labels: list[Label] = field(default_factory=list)
    best_mapping: dict[AtomId, CanonNum] = field(default_factory=dict)

    current_labels: list[Label] = field(default_factory=list)
    current_mapping: dict[AtomId, CanonNum] = field(default_factory=dict)

    frontier: set[AtomId] = field(default_factory=set)

    queue: list[set[AtomId]] = field(default_factory=list)

def better_than(label_like1, label_like2):
    return label_like1 < label_like2

def get_candidates(
    params: _Parameters, state: _CanonState
) -> tuple[set[AtomId], Label]:
    """Return the a set of candidate atoms to investigate.
    The candidates are determined based on the current state of the search:
    """
    n_atoms = len(params.nbrhds)
    frontier = (
        state.frontier
        if state.frontier
        else set(params.nbrhds) - state.current_mapping.keys()
    )

    step_labels = [
        (
            sorted(
                state.current_mapping.get(nbr, n_atoms + 1)
                for nbr in params.nbrhds[atom_id]
            ),
            params.atom_color_frequency[atom_id],
            params.degrees[atom_id],
            params.colors[atom_id],
            atom_id,
        )
        for atom_id in frontier
    ]
    if not step_labels:
        assert len(state.current_mapping) == len(params.nbrhds), (
            len(state.current_mapping),
            len(params.nbrhds),
        )
        return set(), ()

    # Compare using label parts except atom_id (last element),
    # so all atoms with the same minimal label are kept.
    best_step_label = min(lbl[:-1] for lbl in step_labels)

    if (
        state.best_labels
        and state.current_mapping
        and len(state.best_labels) < len(state.current_mapping)
        and better_than(
            state.best_labels[len(state.current_mapping) - 1], best_step_label
        )
    ):
        return set(), ()

    atoms = {lbl[-1] for lbl in step_labels if lbl[:-1] == best_step_label}

    return atoms, best_step_label

def initialize(g: MolGraph) -> tuple[_Parameters, _CanonState]:
    """Prepare immutable inputs and initial state for non-stereo canonical
    enumeration."""
    nbrhd = g.neighbors
    colors = {atom: int(label) for atom, label in zip(g.atoms, color_refine_mg(g))}

    degrees = {atom: len(nbrs) for atom, nbrs in nbrhd.items()}

    label_class_sizes = Counter(colors.values())
    atom_color_frequency = {atom: label_class_sizes[colors[atom]] for atom in g.atoms}

    param = _Parameters(
        nbrhds=nbrhd,
        colors=colors,
        degrees=degrees,
        atom_color_frequency=atom_color_frequency,
    )

    initial_state = _CanonState(frontier=set(g.atoms))

    initial_candidates, initial_label = get_candidates(param, initial_state)

    # Takes first step, because initially frontier contains all atoms
    # and later only the neighbors of the already enumerated ones
    atom = initial_candidates.pop()

    state = _CanonState(
        best_labels=[initial_label],
        best_mapping={atom: 1},
        current_labels=[initial_label],
        current_mapping={atom: 1},
        frontier=nbrhd[atom].copy(),
        queue=[initial_candidates],
    )

    return param, state

def canon_atom_num(g: MolGraph) -> Mapping[AtomId, CanonNum]:
    """Return canonical atom order using best-first branching and tie backtracking."""

    params, state = initialize(g)
    # includes first step.

    while state.queue:
        if len(state.queue[-1]) > len(params.nbrhds) - len(state.current_mapping):
            raise RuntimeError("Too many candidates ?")
        if len(state.queue) > len(params.nbrhds):
            raise RuntimeError("Too many backtracking levels ?")

        # plan step
        if len(state.current_mapping) == len(state.queue):
            # can be empty (set(), tuple()) if there are no candidates.
            new_candidates, label = get_candidates(params, state)
            assert new_candidates.issubset(state.frontier)
            state.queue.append(new_candidates)
            state.current_labels.append(label)

        # execute step back
        if not state.queue[-1]:
            state.queue.pop()
            if state.current_labels:
                state.current_labels.pop()
            atom_id = None
            if state.current_mapping:
                atom_id, _ = state.current_mapping.popitem()
            if atom_id is not None:
                not_frontier_anymore = {
                    nbr
                    for nbr in params.nbrhds[atom_id]
                    if not any(
                        nbr2 in state.current_mapping for nbr2 in params.nbrhds[nbr]
                    )
                }
                state.frontier -= not_frontier_anymore
                state.frontier.add(atom_id)

        # execute step forward
        if (
            state.queue
            and state.queue[-1]
            and len(state.queue) - len(state.current_mapping) == 1
        ):
            next_atom = state.queue[-1].pop()
            state.current_mapping[next_atom] = len(state.current_mapping) + 1
            state.frontier |= params.nbrhds[next_atom]
            state.frontier -= state.current_mapping.keys()

            # update best
            if len(state.current_labels) > len(state.best_labels) or (
                len(state.current_labels) == len(state.best_labels)
                and better_than(state.current_labels, state.best_labels)
            ):
                state.best_labels = state.current_labels.copy()
                state.best_mapping = state.current_mapping.copy()

    return state.best_mapping