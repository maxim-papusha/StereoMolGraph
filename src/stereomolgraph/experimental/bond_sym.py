from __future__ import annotations

from typing import Literal

from typing_extensions import Self

from stereomolgraph import Bond, StereoMolGraph
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.stereodescriptors import AtropBond, OInt, PlanarBond, _StereoMixin


class NonRotatableBond33(
    _StereoMixin[
        tuple[OInt, OInt, OInt, int, int, OInt, OInt, OInt], None | Literal[1, -1]
    ],
):
    r"""
    Represents a bond that cannot freely rotate::
            parity = 1
             0    5
             |    |
        1  ▷ 3 - 4 ◁ 6
            ◀     ▶
           2        7

            parity = -1
             0    5
             |    |
        1  ▷ 3 - 4 ◁ 7
            ◀     ▶
           2        6
    """

    parity = 0
    inversion = (0, 2, 1, 3, 4, 5, 7, 6)
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


class NonRotatableBond23(
    _StereoMixin[tuple[OInt, OInt, int, int, OInt, OInt, OInt], None | Literal[1, -1]],
):
    r"""
    Represents a bond that cannot freely rotate::
            parity = 1
           0     4
            \    |
             2 - 3 ◁ 5
            /     ▶
           1        6

            parity = -1
           0     4
            \    |
             2 - 3 ◁ 6
            /     ▶
           1        5
    """

    parity = 0
    inversion = (0, 1, 2, 3, 4, 6, 5)
    _bond: Bond
    PERMUTATION_GROUP = ((0, 1, 2, 3, 4, 5, 6),)

    def get_isomers(self) -> set[Self]:
        return {self}

    @property
    def bond(self) -> Bond:
        bond = frozenset(self.atoms[2:4])
        assert len(bond) == 2
        return bond


class NonRotatableBond13(
    _StereoMixin[tuple[OInt, int, int, OInt, OInt, OInt], None | Literal[1, -1]],
):
    r"""
    Represents a bond that cannot freely rotate::
            parity = 1
          0     3
           \    |
            1 - 2 ◁ 4
                  ▶
                   5

            parity = -1
          0     3
           \    |
            1 - 2 ◁ 5
                  ▶
                   4
    """

    parity = 0
    inversion = (0, 1, 2, 3, 5, 4)
    _bond: Bond
    PERMUTATION_GROUP = ((0, 1, 2, 3, 4, 5),)

    def get_isomers(self) -> set[Self]:
        return {self}

    @property
    def bond(self) -> Bond:
        bond = frozenset(self.atoms[1:3])
        assert len(bond) == 2
        return bond


NonRotatableBond = (
    NonRotatableBond33
    | NonRotatableBond23
    | NonRotatableBond13
    | PlanarBond
    | AtropBond
)


def get_bond_automorphims_classes(
    smg: StereoMolGraph, only_rotatable: bool = True
) -> set[frozenset[Bond]]:
    """
    Bonds are only considered if both of its atoms have other substituents
    and therefore can have a internal rotation.
    """
    mappings = vf2pp_all_isomorphisms(smg, smg, stereo=True)
    if only_rotatable:
        bonds = {
            bond
            for bond in smg.bonds
            if all(len(smg.bonded_to(atom)) > 2 for atom in bond)
            and not isinstance(smg.get_bond_stereo(bond), NonRotatableBond)
        }
    else:
        bonds = {bond for bond in smg.bonds}

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

    unique_reorderings: set[NonRotatableBond] = set()
    bond_class = s.__class__

    for mapping in mappings:
        if mapping[a1] != a1 or mapping[a2] != a2:
            continue
        map_nrbrs1 = tuple(mapping[nbr] for nbr in nbrs1)
        map_nrbrs2 = tuple(mapping[nbr] for nbr in nbrs2)

        if any(map_nrbr not in nbrs1 for map_nrbr in map_nrbrs1) or any(
            map_nrbr not in nbrs2 for map_nrbr in map_nrbrs2
        ):
            continue

        reordering = bond_class(atoms=(*map_nrbrs1, a1, a2, *map_nrbrs2), parity=1)
        unique_reorderings.add(reordering)

    return len(unique_reorderings)
