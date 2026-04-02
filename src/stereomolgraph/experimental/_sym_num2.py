from __future__ import annotations

import itertools
from collections import Counter
from typing import cast

from stereomolgraph import Bond, StereoMolGraph
from stereomolgraph.algorithms._color_refine import color_refine_smg
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.experimental.new_stereodescriptors import (
    NonRotatableBond13,
    NonRotatableBond23,
    NonRotatableBond33,
)
from stereomolgraph.stereodescriptors import OInt, Tetrahedral

BondAutoCls = int  # bond automorphism class
MappingId = int
NonRotatableBond33Atoms = tuple[OInt, OInt, OInt, int, int, OInt, OInt, OInt]
NonRotatableBond = NonRotatableBond33 | NonRotatableBond23 | NonRotatableBond13


def bond_equivalence_classes(graph: StereoMolGraph) -> dict[Bond, int]:
    """
    Returns a mapping from bonds to their equivalence class.
    Bonds are in the same equivalence class if they can
    """
    if any(stereo.parity is None for stereo in graph.stereo.values()):
        raise NotImplementedError(
            "all stereocenters have to be defined to calculate the symmetry number"
        )
    colorings = color_refine_smg(graph)
    mappings = vf2pp_all_isomorphisms(
        graph,
        graph,
        atom_labels=(colorings, colorings),  # type: ignore[arg-type]
        stereo=True,
    )

    bond_auto_set: dict[Bond, set[Bond]] = {b: {b} for b in graph.bonds}

    for mapping in mappings:
        for b1 in graph.bonds:
            b2 = Bond({mapping[a] for a in b1})
            if b2 in bond_auto_set[b1]:
                continue
            else:
                bond_auto_set[b1].update(bond_auto_set[b2])
                bond_auto_set[b2] = bond_auto_set[b1]

    bond_auto_int: dict[Bond, int] = {
        bond: hash(frozenset(bond_set)) for bond, bond_set in bond_auto_set.items()
    }

    return bond_auto_int


def _require_complete_stereo(graph: StereoMolGraph) -> None:
    if any(stereo.parity is None for stereo in graph.stereo.values()):
        raise NotImplementedError(
            "all stereocenters have to be defined to calculate the symmetry number"
        )


def _self_mappings(graph: StereoMolGraph) -> tuple[dict[int, int], ...]:
    return tuple(
        vf2pp_all_isomorphisms(
            graph,
            graph,
            stereo=True,
        )
    )


def _automorphism_classes(
    graph: StereoMolGraph,
    mappings: tuple[dict[int, int], ...],
) -> tuple[dict[int, int], dict[BondAutoCls, set[Bond]]]:
    atom_auto_set: dict[int, set[int]] = {a: {a} for a in graph.atoms}
    bond_auto_set: dict[Bond, set[Bond]] = {b: {b} for b in graph.bonds}

    for mapping in mappings:
        for a1, a2 in mapping.items():
            if a2 in atom_auto_set[a1]:
                continue
            atom_auto_set[a1].update(atom_auto_set[a2])
            atom_auto_set[a2] = atom_auto_set[a1]

        for b1 in graph.bonds:
            b2 = Bond({mapping[a] for a in b1})
            if b2 in bond_auto_set[b1]:
                continue
            bond_auto_set[b1].update(bond_auto_set[b2])
            bond_auto_set[b2] = bond_auto_set[b1]

    atom_auto_int: dict[int, int] = {
        atom: min(atom_set) for atom, atom_set in atom_auto_set.items()
    }
    bond_auto_int: dict[Bond, int] = {
        bond: hash(frozenset(bond_set)) for bond, bond_set in bond_auto_set.items()
    }
    int_bond_auto_set: dict[BondAutoCls, set[Bond]] = {
        auto_cls_int: bond_auto_set[bond]
        for bond, auto_cls_int in bond_auto_int.items()
    }
    return atom_auto_int, int_bond_auto_set


def _bond_state(
    graph: StereoMolGraph,
    atom_auto_int: dict[int, int],
    bond: Bond,
) -> None | set[NonRotatableBond33]:
    if graph.get_bond_stereo(bond):
        return None

    a1, a2 = tuple(bond)
    nbrs1 = set(graph.bonded_to(a1))
    nbrs1.remove(a2)
    nbrs2 = set(graph.bonded_to(a2))
    nbrs2.remove(a1)
    if len(nbrs1) != 3 or len(nbrs2) != 3:
        return None

    left_auto_classes = {atom_auto_int[nbr] for nbr in nbrs1}
    right_auto_classes = {atom_auto_int[nbr] for nbr in nbrs2}
    nrb_set: set[NonRotatableBond33] = set()

    if len(left_auto_classes) == 2 and len(right_auto_classes) == 2:
        unique_left_int = next(
            auto_cls
            for auto_cls, count in Counter(atom_auto_int[nbr] for nbr in nbrs1).items()
            if count == 1
        )
        unique_right_int = next(
            auto_cls
            for auto_cls, count in Counter(atom_auto_int[nbr] for nbr in nbrs2).items()
            if count == 1
        )
        unique_left = next(
            nbr for nbr in nbrs1 if atom_auto_int[nbr] == unique_left_int
        )
        unique_right = next(
            nbr for nbr in nbrs2 if atom_auto_int[nbr] == unique_right_int
        )

        s1 = graph.get_atom_stereo(a1)
        s2 = graph.get_atom_stereo(a2)
        assert isinstance(s1, Tetrahedral)
        assert isinstance(s2, Tetrahedral)

        left = next(
            (stereo_atoms[2:5] if s1.parity == 1 else stereo_atoms[2:5][::-1])
            for stereo_atoms in s1._perm_atoms()
            if stereo_atoms[1] == a2 and stereo_atoms[2] == unique_left
        )
        right = next(
            (stereo_atoms[2:5] if s2.parity == -1 else stereo_atoms[2:5][::-1])
            for stereo_atoms in s2._perm_atoms()
            if stereo_atoms[1] == a1 and stereo_atoms[2] == unique_right
        )
        nrb_set.add(NonRotatableBond33(atoms=(*left, a1, a2, *right), parity=1))

    elif len(left_auto_classes) == 1 or len(right_auto_classes) == 1:
        s1 = graph.get_atom_stereo(a1)
        s2 = graph.get_atom_stereo(a2)
        assert isinstance(s1, Tetrahedral)
        assert isinstance(s2, Tetrahedral)

        left = next(
            (stereo_atoms[2:5] if s1.parity == 1 else stereo_atoms[2:5][::-1])
            for stereo_atoms in s1._perm_atoms()
            if stereo_atoms[1] == a2
        )
        right = next(
            (stereo_atoms[2:5] if s2.parity == -1 else stereo_atoms[2:5][::-1])
            for stereo_atoms in s2._perm_atoms()
            if stereo_atoms[1] == a1
        )
        nrb_set.add(NonRotatableBond33(atoms=(*left, a1, a2, *right), parity=1))
    else:
        return None

    return nrb_set


def _bond_states(
    graph: StereoMolGraph,
    atom_auto_int: dict[int, int],
) -> dict[Bond, set[NonRotatableBond33]]:
    bond_states: dict[Bond, set[NonRotatableBond33]] = {}
    for bond in graph.bonds:
        state = _bond_state(graph, atom_auto_int, bond)
        if state is not None:
            bond_states[bond] = state
    return bond_states


def _map_non_rotatable_bond(
    non_rot_bond: NonRotatableBond33,
    mapping: dict[int, int],
) -> NonRotatableBond33:
    mapped_atoms = cast(
        NonRotatableBond33Atoms,
        tuple(atom if atom is None else mapping[atom] for atom in non_rot_bond.atoms),
    )
    return NonRotatableBond33(atoms=mapped_atoms, parity=1)


def _bond_state_counter(
    int_bond_auto_set: dict[BondAutoCls, set[Bond]],
    bond_states: dict[Bond, set[NonRotatableBond33]],
    ignored_bond: None | Bond = None,
) -> dict[BondAutoCls, dict[frozenset[NonRotatableBond33], set[MappingId]]]:
    counter: dict[BondAutoCls, dict[frozenset[NonRotatableBond33], set[MappingId]]] = {}

    for bond_auto_cls, bonds in int_bond_auto_set.items():
        counter[bond_auto_cls] = {}
        bond_to_state = {
            bond: bond_states[bond]
            for bond in bonds
            if bond in bond_states and bond != ignored_bond
        }

        for comb in itertools.product(*bond_to_state.values()):
            counter[bond_auto_cls][frozenset(comb)] = set()

    return counter


def _populate_mapping_ids(
    counter: dict[BondAutoCls, dict[frozenset[NonRotatableBond33], set[MappingId]]],
    mappings: tuple[dict[int, int], ...],
) -> None:
    for mapping_id, mapping in enumerate(mappings):
        for bond_state_dict in counter.values():
            for bond_state, mapping_ids in bond_state_dict.items():
                if all(
                    _map_non_rotatable_bond(non_rot_bond, mapping) in bond_state
                    for non_rot_bond in bond_state
                ):
                    mapping_ids.add(mapping_id)


def _max_mapping_intersection(
    counter: dict[BondAutoCls, dict[frozenset[NonRotatableBond33], set[MappingId]]],
    mapping_count: int,
) -> int:
    if not counter:
        return 1

    mapping_sets_by_cls = [
        list(bond_state_dict.values())
        for bond_state_dict in counter.values()
        if bond_state_dict
    ]

    sym_number = 1
    for combo in itertools.product(*mapping_sets_by_cls):
        intersect: set[MappingId] = set.intersection(*combo) if combo else set()
        sym_number = max(sym_number, len(intersect))
        if sym_number == mapping_count:
            break
    return sym_number


def _prepare_external_symmetry_data(
    graph: StereoMolGraph,
) -> tuple[
    tuple[dict[int, int], ...],
    dict[BondAutoCls, set[Bond]],
    dict[Bond, set[NonRotatableBond33]],
]:
    _require_complete_stereo(graph)
    mappings = _self_mappings(graph)
    atom_auto_int, int_bond_auto_set = _automorphism_classes(graph, mappings)
    bond_states = _bond_states(graph, atom_auto_int)
    return mappings, int_bond_auto_set, bond_states


def _external_symmetry_number_from_data(
    mappings: tuple[dict[int, int], ...],
    int_bond_auto_set: dict[BondAutoCls, set[Bond]],
    bond_states: dict[Bond, set[NonRotatableBond33]],
    ignored_bond: None | Bond = None,
) -> int:
    counter = _bond_state_counter(
        int_bond_auto_set,
        bond_states,
        ignored_bond=ignored_bond,
    )
    _populate_mapping_ids(counter, mappings)
    return _max_mapping_intersection(counter, len(mappings))


def external_symmetry_number(graph: StereoMolGraph) -> int:
    """
    TODO: add paper reference.
    """
    mappings, int_bond_auto_set, bond_states = _prepare_external_symmetry_data(graph)
    return _external_symmetry_number_from_data(
        mappings,
        int_bond_auto_set,
        bond_states,
    )


def external_symmetry_number_bond(graph: StereoMolGraph, bond: Bond) -> int:
    """
    External symmetry number with `bond` treated as the single rotatable bond.
    All other modeled bonds remain constrained.
    """
    bond = Bond(bond)
    if bond not in graph.bonds:
        raise ValueError("bond must belong to the graph")

    mappings, int_bond_auto_set, bond_states = _prepare_external_symmetry_data(graph)
    return _external_symmetry_number_from_data(
        mappings,
        int_bond_auto_set,
        bond_states,
        ignored_bond=bond,
    )


def external_symmetry_number_per_rotatable_bond(
    graph: StereoMolGraph,
) -> dict[Bond, int]:
    """
    External symmetry number for each modeled rotatable bond, treating that
    bond as rotatable and all other modeled bonds as constrained.
    """
    mappings, int_bond_auto_set, bond_states = _prepare_external_symmetry_data(graph)
    return {
        bond: _external_symmetry_number_from_data(
            mappings,
            int_bond_auto_set,
            bond_states,
            ignored_bond=bond,
        )
        for bond in sorted(
            bond_states,
            key=lambda current_bond: tuple(sorted(current_bond)),
        )
    }


def get_bond_automorphims_classes(smg: StereoMolGraph) -> set[frozenset[Bond]]:
    """
    Bonds are only considered if both of its atoms have other substituents
    and therefore can have a internal rotation.
    """
    mappings = vf2pp_all_isomorphisms(smg, smg, stereo=True)
    bonds = {
        bond for bond in smg.bonds if all(len(smg.bonded_to(atom)) > 2 for atom in bond)
    }

    parent: dict[Bond, Bond] = {bond: bond for bond in bonds}
    rank: dict[Bond, int] = {bond: 0 for bond in bonds}

    def find(bond: Bond) -> Bond:
        root = bond
        while parent[root] != root:
            root = parent[root]
        while parent[bond] != bond:
            next_bond = parent[bond]
            parent[bond] = root
            bond = next_bond
        return root

    def union(bond1: Bond, bond2: Bond) -> None:
        root1 = find(bond1)
        root2 = find(bond2)
        if root1 == root2:
            return
        if rank[root1] < rank[root2]:
            root1, root2 = root2, root1
        parent[root2] = root1
        if rank[root1] == rank[root2]:
            rank[root1] += 1

    for mapping in mappings:
        for bond in bonds:
            mapped_bond = Bond(mapping[atom] for atom in bond)
            if mapped_bond in parent:
                union(bond, mapped_bond)

    classes: dict[Bond, set[Bond]] = {}
    for bond in bonds:
        root = find(bond)
        classes.setdefault(root, set()).add(bond)

    return {frozenset(bond_class) for bond_class in classes.values()}


def bond_symmetry_number(graph: StereoMolGraph, bond: Bond) -> int:
    mappings = vf2pp_all_isomorphisms(graph, graph, stereo=True)

    a1, a2 = bond

    nbrs1 = tuple(a for a in graph.bonded_to(a1) if a != a2)
    nbrs2 = tuple(a for a in graph.bonded_to(a2) if a != a1)

    if len(nbrs1) == 0 or len(nbrs2) == 0:
        return 1

    if len(nbrs1) > len(nbrs2):
        a1, a2 = a2, a1
        nbrs1, nbrs2 = nbrs2, nbrs1

    if len(nbrs1) == 3 and len(nbrs2) == 3:
        s = NonRotatableBond33(atoms=(*nbrs1, a1, a2, *nbrs2), parity=1)
    elif len(nbrs1) == 2 and len(nbrs2) == 3:
        s = NonRotatableBond23(atoms=(*nbrs1, a1, a2, *nbrs2), parity=1)
    elif len(nbrs1) == 1 and len(nbrs2) == 3:
        s = NonRotatableBond13(atoms=(*nbrs1, a1, a2, *nbrs2), parity=1)
    else:
        raise NotImplementedError

    # diff nrbs1
    unique_reorderings: set[NonRotatableBond] = set()
    bond_class = s.__class__

    for mapping in mappings:
        if mapping[a1] != a1 or mapping[a2] != a2:
            continue
        map_nrbrs1 = tuple(mapping[nbr] for nbr in nbrs1)
        map_nrbrs2 = tuple(mapping[nbr] for nbr in nbrs2)

        if (map_nrbr not in nbrs1 for map_nrbr in map_nrbrs1) or any(
            map_nrbr not in nbrs2 for map_nrbr in map_nrbrs2
        ):
            continue

        reordering = bond_class(atoms=(*map_nrbrs1, a1, a2, *map_nrbrs2), parity=1)
        unique_reorderings.add(reordering)

    return len(unique_reorderings)
