from __future__ import annotations

import itertools
from collections import Counter

from stereomolgraph import Bond, StereoMolGraph
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.experimental.bond_sym import (
    AtropBond,
    NonRotatableBond13,
    NonRotatableBond23,
    NonRotatableBond33,
    PlanarBond,
)
from stereomolgraph.stereodescriptors import Tetrahedral

NonRotatableBond = (
    NonRotatableBond33
    | NonRotatableBond23
    | NonRotatableBond13
    | PlanarBond
    | AtropBond
)
GeneratedBondStereo = NonRotatableBond33 | NonRotatableBond23 | NonRotatableBond13

BondAutoCls = int  # bond automorphism class
MappingId = int


def external_symmetry_number(graph: StereoMolGraph) -> int:
    """
    0. Get automorphic mappings
    1. Group bonds in equivalence classes
     2. Check for kind of bond. Ignore if it is already assigned as
         PlanarBond or AtropBond.
       Depending on number of left and right substituents of the bond use
       3 substituents, 3 substituents: NonRotatableBond33
       2 substituents, 3 substituents: NonRotatableBond23
       1 substituent, 3 substituents: NonRotatableBond13
    3. Assign the most rotationally symetric state for one of each automorphism class
    4. Assign remaining bonds of each equivalence class with mapping from other.
    5. Reevaluate mappings including new bond stereodescriptors.
    """
    if any(stereo.parity is None for stereo in graph.stereo.values()):
        raise NotImplementedError(
            "all stereocenters have to be defined to calculate the symmetry number"
        )

    mappings = tuple(
        vf2pp_all_isomorphisms(
            graph,
            graph,
            stereo=True,
        )
    )

    atom_auto_set: dict[int, set[int]] = {a: {a} for a in graph.atoms}
    bond_auto_set: dict[Bond, set[Bond]] = {b: {b} for b in graph.bonds}

    # Step 1: build automorphism classes for atoms and bonds.
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
    int_bond_auto_set: dict[int, set[Bond]] = {
        auto_cls_int: bond_auto_set[bond]
        for bond, auto_cls_int in bond_auto_int.items()
    }

    def _map_optional_atom(atom: int | None, mapping: dict[int, int]) -> int | None:
        if atom is None:
            return None
        return mapping[atom]

    def _ordered_three_neighbors(
        atom: int,
        bonded_atom: int,
        *,
        is_left: bool,
        unique_neighbor: int | None = None,
    ) -> tuple[int | None, int | None, int | None]:
        stereo = graph.get_atom_stereo(atom)
        if isinstance(stereo, Tetrahedral):
            for stereo_atoms in stereo._perm_atoms():
                if stereo_atoms[1] != bonded_atom:
                    continue
                if unique_neighbor is not None and stereo_atoms[2] != unique_neighbor:
                    continue

                ordered = stereo_atoms[2:5]
                if (is_left and stereo.parity == 1) or (
                    not is_left and stereo.parity == -1
                ):
                    return ordered[0], ordered[1], ordered[2]

                reversed_ordered = ordered[::-1]
                return (
                    reversed_ordered[0],
                    reversed_ordered[1],
                    reversed_ordered[2],
                )

        neighbors = tuple(sorted(n for n in graph.bonded_to(atom) if n != bonded_atom))
        return neighbors[0], neighbors[1], neighbors[2]

    def _build_non_rotatable_state(bond: Bond) -> GeneratedBondStereo | None:
        bond_stereo = graph.get_bond_stereo(bond)
        if isinstance(bond_stereo, (PlanarBond, AtropBond)):
            return None
        if bond_stereo is not None:
            return None

        a1, a2 = sorted(bond)
        nbrs1 = tuple(sorted(n for n in graph.bonded_to(a1) if n != a2))
        nbrs2 = tuple(sorted(n for n in graph.bonded_to(a2) if n != a1))

        if len(nbrs1) > len(nbrs2):
            a1, a2 = a2, a1
            nbrs1, nbrs2 = nbrs2, nbrs1

        side_counts = (len(nbrs1), len(nbrs2))
        if side_counts == (3, 3):
            left_counts = Counter(atom_auto_int[n] for n in nbrs1)
            right_counts = Counter(atom_auto_int[n] for n in nbrs2)
            unique_left = next(
                (n for n in nbrs1 if left_counts[atom_auto_int[n]] == 1),
                None,
            )
            unique_right = next(
                (n for n in nbrs2 if right_counts[atom_auto_int[n]] == 1),
                None,
            )

            left3 = _ordered_three_neighbors(
                a1,
                a2,
                is_left=True,
                unique_neighbor=unique_left,
            )
            right3 = _ordered_three_neighbors(
                a2,
                a1,
                is_left=False,
                unique_neighbor=unique_right,
            )

            return NonRotatableBond33(
                atoms=(*left3, a1, a2, *right3),
                parity=1,
            )

        if side_counts == (2, 3):
            right_counts = Counter(atom_auto_int[n] for n in nbrs2)
            unique_right = next(
                (n for n in nbrs2 if right_counts[atom_auto_int[n]] == 1),
                None,
            )
            right3 = _ordered_three_neighbors(
                a2,
                a1,
                is_left=False,
                unique_neighbor=unique_right,
            )
            left2 = (nbrs1[0], nbrs1[1])

            return NonRotatableBond23(
                atoms=(*left2, a1, a2, *right3),
                parity=1,
            )

        if side_counts == (1, 3):
            right_counts = Counter(atom_auto_int[n] for n in nbrs2)
            unique_right = next(
                (n for n in nbrs2 if right_counts[atom_auto_int[n]] == 1),
                None,
            )
            right3 = _ordered_three_neighbors(
                a2,
                a1,
                is_left=False,
                unique_neighbor=unique_right,
            )

            return NonRotatableBond13(
                atoms=(nbrs1[0], a1, a2, *right3),
                parity=1,
            )

        return None

    def _map_non_rotatable_state(
        state: GeneratedBondStereo,
        mapping: dict[int, int],
    ) -> GeneratedBondStereo:
        if isinstance(state, NonRotatableBond33):
            atoms = state.atoms
            mapped = (
                _map_optional_atom(atoms[0], mapping),
                _map_optional_atom(atoms[1], mapping),
                _map_optional_atom(atoms[2], mapping),
                mapping[atoms[3]],
                mapping[atoms[4]],
                _map_optional_atom(atoms[5], mapping),
                _map_optional_atom(atoms[6], mapping),
                _map_optional_atom(atoms[7], mapping),
            )
            return NonRotatableBond33(atoms=mapped, parity=1)

        if isinstance(state, NonRotatableBond23):
            atoms = state.atoms
            mapped = (
                _map_optional_atom(atoms[0], mapping),
                _map_optional_atom(atoms[1], mapping),
                mapping[atoms[2]],
                mapping[atoms[3]],
                _map_optional_atom(atoms[4], mapping),
                _map_optional_atom(atoms[5], mapping),
                _map_optional_atom(atoms[6], mapping),
            )
            return NonRotatableBond23(atoms=mapped, parity=1)

        atoms = state.atoms
        mapped = (
            _map_optional_atom(atoms[0], mapping),
            mapping[atoms[1]],
            mapping[atoms[2]],
            _map_optional_atom(atoms[3], mapping),
            _map_optional_atom(atoms[4], mapping),
            _map_optional_atom(atoms[5], mapping),
        )
        return NonRotatableBond13(atoms=mapped, parity=1)

    counter: dict[
        BondAutoCls,
        dict[frozenset[GeneratedBondStereo], set[MappingId]],
    ] = {}

    for bond_auto_cls, bonds in int_bond_auto_set.items():
        counter[bond_auto_cls] = {}
        states_by_bond: dict[Bond, GeneratedBondStereo] = {}

        # Step 2: classify bond type and skip already assigned Planar/Atrop bonds.
        for bond in bonds:
            state = _build_non_rotatable_state(bond)
            if state is not None:
                states_by_bond[bond] = state

        if not states_by_bond:
            continue

        representative = min(tuple(sorted(bond)) for bond in states_by_bond)
        representative_bond = Bond(representative)
        representative_state = states_by_bond[representative_bond]

        # Step 3: assign one most symmetric state for a representative bond.
        # Step 4: assign states to equivalent bonds by mapping from the representative.
        bond_to_state: dict[Bond, set[GeneratedBondStereo]] = {
            bond: set() for bond in states_by_bond
        }
        bond_to_state[representative_bond].add(representative_state)

        for mapping in mappings:
            mapped_bond = Bond(mapping[a] for a in representative_bond)
            if mapped_bond not in bond_to_state:
                continue

            mapped_state = _map_non_rotatable_state(representative_state, mapping)
            if mapped_state.bond == mapped_bond:
                bond_to_state[mapped_bond].add(mapped_state)

        # Defensive fallback: keep at least one state for each bond.
        for bond, states in bond_to_state.items():
            if not states:
                states.add(states_by_bond[bond])

        for comb in itertools.product(*bond_to_state.values()):
            counter[bond_auto_cls][frozenset(comb)] = set()

    # Step 5: reevaluate mappings by requiring mapped bond stereodescriptors
    # to remain inside the selected state set.
    for mapping_id, mapping in enumerate(mappings):
        for bond_state_dict in counter.values():
            for bond_state, mapping_ids in bond_state_dict.items():
                if all(
                    _map_non_rotatable_state(non_rot_bond, mapping) in bond_state
                    for non_rot_bond in bond_state
                ):
                    mapping_ids.add(mapping_id)

    mapping_sets_by_cls = [
        list(bond_state_dict.values())
        for bond_state_dict in counter.values()
        if bond_state_dict
    ]
    if not mapping_sets_by_cls:
        return 1

    sym_number = 1
    for combo in itertools.product(*mapping_sets_by_cls):
        intersection = set.intersection(*combo) if combo else set()
        sym_number = max(sym_number, len(intersection))
        if sym_number == len(mappings):
            break

    return sym_number


def old_external_symmetry_number(graph: StereoMolGraph) -> int:

    if any(stereo.parity is None for stereo in graph.stereo.values()):
        raise NotImplementedError(
            "all stereocenters have to be defined to calculate the symmetry number"
        )
    # colorings = color_refine_smg(graph)
    mappings = tuple(
        vf2pp_all_isomorphisms(
            graph,
            graph,
            stereo=True,
            # atom_labels=(colorings, colorings),
        )
    )

    atom_auto_set: dict[int, set[int]] = {a: {a} for a in graph.atoms}

    bond_auto_set: dict[Bond, set[Bond]] = {b: {b} for b in graph.bonds}

    for mapping in mappings:
        for a1, a2 in mapping.items():
            if a2 in atom_auto_set[a1]:
                continue
            else:
                atom_auto_set[a1].update(atom_auto_set[a2])
                atom_auto_set[a2] = atom_auto_set[a1]

        for b1 in graph.bonds:
            b2 = Bond({mapping[a] for a in b1})
            if b2 in bond_auto_set[b1]:
                continue
            else:
                bond_auto_set[b1].update(bond_auto_set[b2])
                bond_auto_set[b2] = bond_auto_set[b1]

    atom_auto_int: dict[int, int] = {
        a: min(atom_set) for a, atom_set in atom_auto_set.items()
    }
    bond_auto_int: dict[Bond, int] = {
        bond: hash(frozenset(bond_set)) for bond, bond_set in bond_auto_set.items()
    }

    # int_atom_auto_set: dict[int, set[int]] = {auto_cls_int: atom_auto_set[a]
    #                                          for a, auto_cls_int
    #                                          in atom_auto_int.items()}
    int_bond_auto_set: dict[int, set[Bond]] = {
        auto_cls_int: bond_auto_set[b] for b, auto_cls_int in bond_auto_int.items()
    }

    counter: dict[BondAutoCls, dict[frozenset[NonRotatableBond], set[MappingId]]] = {}

    for bond_auto_cls, bonds in int_bond_auto_set.items():
        counter[bond_auto_cls] = {}

        bond_to_state: dict[Bond, set[NonRotatableBond]] = {}

        for bond in bonds:
            if graph.get_bond_stereo(bond):
                continue
            a1, a2 = tuple(bond)
            nbrs1 = set(graph.bonded_to(a1))
            nbrs1.remove(a2)
            nbrs2 = set(graph.bonded_to(a2))
            nbrs2.remove(a1)
            if len(nbrs1) != 3 or len(nbrs2) != 3:
                continue

            nrb_set: set[NonRotatableBond] = set()
            if (
                len({atom_auto_int[nbr] for nbr in nbrs1}) == 2
                and len({atom_auto_int[nbr] for nbr in nbrs2}) == 2
            ):
                unique_left_int = next(
                    v
                    for v, c in Counter(atom_auto_int[x] for x in nbrs1).items()
                    if c == 1
                )
                unique_right_int = next(
                    v
                    for v, c in Counter(atom_auto_int[x] for x in nbrs2).items()
                    if c == 1
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

                nrb_set.add(NonRotatableBond(atoms=(*left, a1, a2, *right), parity=0))

            elif (
                len({atom_auto_int[nbr] for nbr in nbrs1}) == 1
                or len({atom_auto_int[nbr] for nbr in nbrs2}) == 1
            ):
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

                nrb_set.add(NonRotatableBond(atoms=(*left, a1, a2, *right), parity=0))

                # for perm_nbrs1 in itertools.permutations(nbrs1):
                #    for perm_nbrs2 in itertools.permutations(nbrs2):
                #        nrb = NonRotatableBond.from_bond_and_nrbs(perm_nbrs1, a1, a2, perm_nbrs2)
                #        nrb_set.add(nrb)
            else:
                continue
                raise Exception(
                    len({atom_auto_int[nbr] for nbr in nbrs1}),
                    len({atom_auto_int[nbr] for nbr in nbrs2}),
                )

            bond_to_state[bond] = nrb_set

        for _, comb in enumerate(itertools.product(*bond_to_state.values())):
            counter[bond_auto_cls][frozenset(comb)] = set()

    # reevalutates the automorphic mappings and removes internal rotations
    for mapping_id, mapping in enumerate(mappings):
        for bond_auto_cls, bond_state_dict in counter.items():
            for bond_state, mapping_ids in bond_state_dict.items():
                if all(
                    NonRotatableBond(
                        tuple(mapping[a] for a in non_rot_bond.atoms), parity=0
                    )
                    in bond_state
                    for non_rot_bond in bond_state
                ):
                    mapping_ids.add(mapping_id)

    # pick one bond state per auto class, maximize the intersection size of mapping ids
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
        if sym_number == len(mappings):
            break

    return sym_number
