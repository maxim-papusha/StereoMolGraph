import pytest

from stereomolgraph.experimental import canonical_atom_order
from stereomolgraph.graphs import MolGraph


def _cycle_graph(atom_ids: tuple[int, ...]) -> MolGraph:
    graph = MolGraph()
    for atom in atom_ids:
        graph.add_atom(atom, "C")

    for atom1, atom2 in zip(atom_ids, atom_ids[1:] + atom_ids[:1]):
        graph.add_bond(atom1, atom2)

    return graph


def _path_graph(atom_ids: tuple[int, ...], atom_types: tuple[str, ...]) -> MolGraph:
    graph = MolGraph()
    for atom, atom_type in zip(atom_ids, atom_types):
        graph.add_atom(atom, atom_type)

    for atom1, atom2 in zip(atom_ids, atom_ids[1:]):
        graph.add_bond(atom1, atom2)

    return graph


def _canonical_signature(graph: MolGraph) -> tuple[tuple[int, ...], tuple[int, ...]]:
    order = canonical_atom_order(graph)

    atom_types = tuple(int(graph.get_atom_type(atom)) for atom in order)
    bond_bits = []
    for idx, atom1 in enumerate(order):
        neighbors = graph.bonded_to(atom1)
        for atom2 in order[idx + 1 :]:
            bond_bits.append(1 if atom2 in neighbors else 0)

    return atom_types, tuple(bond_bits)


def test_canonical_order_contains_all_atoms() -> None:
    graph = _cycle_graph((10, 20, 30, 40, 50, 60))

    order = canonical_atom_order(graph)

    assert len(order) == graph.n_atoms
    assert set(order) == set(graph.atoms)


@pytest.mark.parametrize(
    "mapping",
    [
        {10: 101, 20: 17, 30: 88, 40: 3, 50: 55, 60: 209},
        {10: 209, 20: 55, 30: 3, 40: 88, 50: 17, 60: 101},
        {10: 55, 20: 209, 30: 101, 40: 17, 50: 88, 60: 3},
    ],
)
def test_canonical_signature_is_invariant_to_relabeling(
    mapping: dict[int, int],
) -> None:
    graph = _cycle_graph((10, 20, 30, 40, 50, 60))
    expected_signature = _canonical_signature(graph)

    relabeled = graph.relabel_atoms(mapping, copy=True)

    assert _canonical_signature(relabeled) == expected_signature


def test_atom_types_are_part_of_canonical_ordering() -> None:
    carbon_path = _path_graph((10, 20, 30, 40), ("C", "C", "C", "C"))
    hetero_path = _path_graph((100, 200, 300, 400), ("C", "N", "C", "C"))

    assert _canonical_signature(carbon_path) != _canonical_signature(hetero_path)
