from __future__ import annotations

import heapq
import itertools
from collections import defaultdict
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from stereomolgraph.algorithms.color_refine import color_refine_mg

if TYPE_CHECKING:
    from stereomolgraph.graphs import MolGraph
    from stereomolgraph.graphs.mg import AtomId


@dataclass(order=True)
class _QueueState:
    priority: tuple[int, tuple[int, ...], int]
    serial: int
    colors: np.ndarray = field(compare=False)
    depth: int = field(compare=False)


def _normalize_colors(colors: np.ndarray) -> np.ndarray:
    unique = np.unique(colors)
    normalized = np.searchsorted(unique, colors)
    return normalized.astype(np.int64, copy=False)


def _color_classes(colors: np.ndarray) -> list[np.ndarray]:
    grouped: dict[int, list[int]] = defaultdict(list)
    for idx, color in enumerate(colors):
        grouped[int(color)].append(idx)
    return [np.array(grouped[color], dtype=np.int64) for color in sorted(grouped)]


def _state_priority(colors: np.ndarray, depth: int) -> tuple[int, tuple[int, ...], int]:
    class_sizes = tuple(
        sorted((int(class_.size) for class_ in _color_classes(colors)), reverse=True)
    )
    return (-len(class_sizes), class_sizes, -depth)


def _adjacency_matrix(graph: MolGraph, atoms: tuple[AtomId, ...]) -> np.ndarray:
    atom_idx = {atom: idx for idx, atom in enumerate(atoms)}
    adjacency = np.zeros((len(atoms), len(atoms)), dtype=np.int8)

    for bond in graph.bonds:
        atom1, atom2 = tuple(bond)
        idx1 = atom_idx[atom1]
        idx2 = atom_idx[atom2]
        adjacency[idx1, idx2] = 1
        adjacency[idx2, idx1] = 1

    return adjacency


def _canonical_code(
    order: tuple[int, ...],
    atom_types: np.ndarray,
    adjacency: np.ndarray,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    n_atoms = len(order)
    type_code = tuple(int(atom_types[idx]) for idx in order)

    if n_atoms < 2:
        return type_code, tuple()

    permuted = adjacency[np.ix_(order, order)]
    upper = permuted[np.triu_indices(n_atoms, k=1)]
    bond_code = tuple(int(value) for value in upper.tolist())

    return type_code, bond_code


def _individualize_and_refine(
    graph: MolGraph,
    colors: np.ndarray,
    atom_idx: int,
) -> np.ndarray:
    individualized = colors.copy()
    individualized[atom_idx] = int(np.max(individualized)) + 1
    refined = color_refine_mg(graph, atom_labels=individualized)
    return _normalize_colors(refined)


def canonical_atom_order(graph: MolGraph) -> tuple[AtomId, ...]:
    """Return a canonical atom enumeration for a MolGraph.

    The algorithm starts with color refinement (1-WL), then performs a
    best-first individualization/refinement search over remaining ambiguous
    color classes. All branches are explored via backtracking, ensuring that
    the globally best canonical labeling is returned.
    """

    atoms = tuple(graph.atoms)
    n_atoms = len(atoms)
    if n_atoms == 0:
        return tuple()

    atom_types = np.array(
        [int(graph.get_atom_type(atom)) for atom in atoms],
        dtype=np.int64,
    )
    adjacency = _adjacency_matrix(graph, atoms)

    initial_colors = _normalize_colors(color_refine_mg(graph))

    serial_counter = itertools.count()
    frontier: list[_QueueState] = []
    seen: set[tuple[int, ...]] = set()

    heapq.heappush(
        frontier,
        _QueueState(
            priority=_state_priority(initial_colors, depth=0),
            serial=next(serial_counter),
            colors=initial_colors,
            depth=0,
        ),
    )

    best_code: None | tuple[tuple[int, ...], tuple[int, ...]] = None
    best_order: None | tuple[int, ...] = None

    while frontier:
        state = heapq.heappop(frontier)
        state_key = tuple(int(color) for color in state.colors.tolist())
        if state_key in seen:
            continue
        seen.add(state_key)

        classes = _color_classes(state.colors)

        if all(class_.size == 1 for class_ in classes):
            order = tuple(int(class_[0]) for class_ in classes)
            code = _canonical_code(order, atom_types, adjacency)
            if best_code is None or code > best_code:
                best_code = code
                best_order = order
            continue

        branch_class = min(
            (class_ for class_ in classes if class_.size > 1),
            key=lambda class_: int(class_.size),
        )

        children: list[tuple[tuple[int, ...], np.ndarray]] = []
        for atom_idx in branch_class.tolist():
            child_colors = _individualize_and_refine(
                graph,
                state.colors,
                atom_idx=int(atom_idx),
            )
            child_sizes = tuple(
                sorted(
                    (
                        int(child_class.size)
                        for child_class in _color_classes(child_colors)
                    ),
                    reverse=True,
                )
            )
            children.append((child_sizes, child_colors))

        children.sort(key=lambda item: item[0])

        for _child_sizes, child_colors in children:
            child_key = tuple(int(color) for color in child_colors.tolist())
            if child_key in seen:
                continue
            child_depth = state.depth + 1
            heapq.heappush(
                frontier,
                _QueueState(
                    priority=_state_priority(child_colors, depth=child_depth),
                    serial=next(serial_counter),
                    colors=child_colors,
                    depth=child_depth,
                ),
            )

    assert best_order is not None
    return tuple(atoms[idx] for idx in best_order)


def canonical_relabel_mapping(graph: MolGraph) -> dict[AtomId, int]:
    """Return mapping from current atom ids to canonical atom ids."""
    return {atom: idx for idx, atom in enumerate(canonical_atom_order(graph))}
