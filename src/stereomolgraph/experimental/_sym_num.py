from __future__ import annotations

import itertools
from collections import Counter, deque
from typing import Literal

from typing_extensions import Self

from stereomolgraph import Bond, StereoMolGraph
from stereomolgraph.algorithms._color_refine import color_refine_smg
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.stereodescriptors import OInt, Tetrahedral, _StereoMixin

BondAutoCls = int  # bond automorphism class
MappingId = int


def topological_symmetry_number(graph: StereoMolGraph) -> int:
    """
    Calculated from the number of graph isomorphisms which conserve the
    stereo information.
    symmetry_number = internal_symmetry_number * rotational_symmetry_number
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
    return deque(enumerate(mappings, 1), maxlen=1)[0][0]


class NonRotatableBond(
    _StereoMixin[
        tuple[OInt, OInt, OInt, int, int, OInt, OInt, OInt], None | Literal[0]
    ],
):
    r"""
    Represents a bond that cannot freely rotate

             0    5
             |    |
        1  ▷ 3 - 4 ◁ 6
            ◀     ▶
           2        7
    """

    parity = 0
    inversion = None
    _bond: Bond
    PERMUTATION_GROUP = (
        (0, 1, 2, 3, 4, 5, 6, 7),
        (5, 7, 6, 4, 3, 0, 2, 1),
        (1, 2, 0, 3, 4, 6, 7, 5),
        (6, 5, 7, 4, 3, 1, 0, 2),
        (2, 0, 1, 3, 4, 7, 5, 6),
        (7, 6, 5, 4, 3, 2, 1, 0),
    )

    def get_isomers(self) -> set[Self]:
        return {self}

    @property
    def bond(self) -> Bond:
        bond = frozenset(self.atoms[3:5])
        assert len(bond) == 2
        return bond


class NonRotatableBond23: ...


def external_symmetry_number(graph: StereoMolGraph) -> int:
    """
    TODO: add paper reference.
    """
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
