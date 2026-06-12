from __future__ import annotations

from collections import deque
from collections.abc import Callable, Hashable, Iterable, Mapping
from typing import TypeVar

from stereomolgraph import AtomId, Bond, StereoMolGraph
from stereomolgraph.algorithms.circular import color_refine_smg
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.stereodescriptors import (
    HinderedBond,
    HinderedBond12,
    HinderedBond13,
    HinderedBond23,
    HinderedBond33,
    OInt,
    PlanarBond,
    Tetrahedral,
)

ItemT = TypeVar("ItemT", bound=Hashable)


def get_automorphism_classes(
    items: Iterable[ItemT],
    mappings: Iterable[Mapping[AtomId, AtomId]],
    map_item: Callable[[ItemT, Mapping[AtomId, AtomId]], ItemT],
) -> dict[ItemT, set[ItemT]]:
    items = tuple(items)

    parent: dict[ItemT, ItemT] = {item: item for item in items}
    rank: dict[ItemT, int] = {item: 0 for item in items}

    def find(item: ItemT) -> ItemT:
        root = item
        while parent[root] != root:
            root = parent[root]
        while parent[item] != item:
            next_item = parent[item]
            parent[item] = root
            item = next_item
        return root

    def union(item1: ItemT, item2: ItemT) -> None:
        root1 = find(item1)
        root2 = find(item2)
        if root1 == root2:
            return
        if rank[root1] < rank[root2]:
            root1, root2 = root2, root1
        parent[root2] = root1
        if rank[root1] == rank[root2]:
            rank[root1] += 1

    for mapping in mappings:
        for item in items:
            mapped_item = map_item(item, mapping)
            if mapped_item in parent:
                union(item, mapped_item)

    classes: dict[ItemT, set[ItemT]] = {}
    for item in items:
        root = find(item)
        classes.setdefault(root, set()).add(item)

    return {item: set(classes[find(item)]) for item in items}


def atom_automorphism_classes(
    smg: StereoMolGraph,
    mappings: Iterable[Mapping[AtomId, AtomId]] | None = None,
) -> dict[AtomId, set[AtomId]]:
    if mappings is None:
        mappings = vf2pp_all_isomorphisms(smg, smg, stereo=True)

    return get_automorphism_classes(
        smg.atoms,
        mappings,
        lambda atom, mapping: mapping[atom],
    )


def bond_automorphism_classes(
    smg: StereoMolGraph,
    mappings: Iterable[Mapping[AtomId, AtomId]] | None = None,
) -> dict[Bond, set[Bond]]:
    """
    Bonds are only considered if both of its atoms have other substituents
    and therefore can have a internal rotation.
    """
    if mappings is None:
        mappings = vf2pp_all_isomorphisms(smg, smg, stereo=True)

    classes = get_automorphism_classes(
        {bond for bond in smg.bonds},
        mappings,
        lambda bond, mapping: Bond(mapping[atom] for atom in bond),
    )

    return classes


def topological_symmetry_number(graph: StereoMolGraph, atom_labels=None) -> int:
    """
    Calculated from the number of graph isomorphisms which conserve the
    stereo information.
    symmetry_number = internal_symmetry_number * rotational_symmetry_number
    """
    if atom_labels is None:
        atom_labels = color_refine_smg(graph, atom_labels=atom_labels)

    mappings = vf2pp_all_isomorphisms(
        graph, graph, stereo=True, atom_labels=(atom_labels, atom_labels)
    )
    return deque(enumerate(mappings, 1), maxlen=1)[0][0]


def bond_symmetry_number(
    graph: StereoMolGraph,
    bond: Bond,
    mappings: Iterable[Mapping[int, int]] | None = None,
) -> int:
    if Bond(bond) not in graph.bonds:
        raise ValueError("Bond has to be part of the graph")

    if mappings is None:
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
        s = HinderedBond33(atoms=(*nbrs1, a1, a2, *nbrs2), parity=1)
    elif len(nbrs1) == 2 and len(nbrs2) == 3:
        s = HinderedBond23(atoms=(*nbrs1, a1, a2, *nbrs2), parity=1)
    elif len(nbrs1) == 1 and len(nbrs2) == 3:
        s = HinderedBond13(atoms=(*nbrs1, a1, a2, *nbrs2), parity=1)
    elif len(nbrs1) == 2 and len(nbrs2) == 2:
        s = PlanarBond(atoms=(*nbrs1, a1, a2, *nbrs2), parity=0)
    elif len(nbrs1) == 1 and len(nbrs2) == 2:
        s = HinderedBond12(atoms=(*nbrs1, a1, a2, *nbrs2), parity=0)
    else:
        raise NotImplementedError(
            "Bonds with more than 3 substituents are not supported"
        )
    unique_reorderings: set[HinderedBond] = set()
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


def ext_sym_num(
    graph: StereoMolGraph, mappings: Iterable[dict[int, int]] | None = None
) -> int:
    """Calculate the upper bound of the external symmetry number for StereoMolGraph"""

    if any(stereo.parity is None for stereo in graph.stereo.values()):
        raise ValueError("Stereodescriptor parity has to be assigned")

    if mappings is None:
        mappings = vf2pp_all_isomorphisms(g1=graph, g2=graph, stereo=True)
    mappings: tuple[dict[int | None, int | None], ...] = tuple(mappings)
    for mapping in mappings:
        mapping[None] = None

    # atom_eq_classes: dict[AtomId | None, set[AtomId | None]]
    atom_eq_classes = atom_automorphism_classes(graph, mappings=mappings)
    atom_eq_classes[None] = {None}
    bond_eq_classes = bond_automorphism_classes(graph, mappings=mappings)

    def group_by_eq(atoms: Iterable[OInt]) -> list[set[OInt]]:
        groups: list[set[OInt]] = []
        for a in atoms:
            for g in groups:
                if any(a in atom_eq_classes[b] for b in g):
                    g.add(a)
                    break
            else:
                groups.append({a})
        return groups

    def ordered_tetrahedral_neighbors(
        stereo: Tetrahedral,
        bonded_atom: AtomId,
        *,
        is_left: bool,
        unique_neighbor: OInt = None,
    ) -> tuple[OInt, OInt, OInt]:
        for permuted_atoms in stereo._perm_atoms():
            if permuted_atoms[1] != bonded_atom:
                continue
            if unique_neighbor is not None and permuted_atoms[2] != unique_neighbor:
                continue

            ordered_neighbors = (
                permuted_atoms[2],
                permuted_atoms[3],
                permuted_atoms[4],
            )
            if (is_left and stereo.parity == 1) or (
                not is_left and stereo.parity == -1
            ):
                return ordered_neighbors

            return (
                ordered_neighbors[2],
                ordered_neighbors[1],
                ordered_neighbors[0],
            )

        raise ValueError(
            f"Could not order tetrahedral neighbors for atom {stereo.central_atom}"
        )

    def singled_out_neighbor(groups: list[set[OInt]]) -> OInt:
        if len(groups) != 2:
            return None

        singleton_group = next((group for group in groups if len(group) == 1), None)
        if singleton_group is None:
            return None

        return next(iter(singleton_group))

    hindered_bonds: dict[Bond, HinderedBond] = {}

    for eq_cls in {frozenset(eq_cls) for eq_cls in bond_eq_classes.values()}:
        eq_bonds = tuple(eq_cls)
        eq_cls_iterator = iter(eq_bonds)
        a1, a2 = next(eq_cls_iterator)

        if graph.get_bond_stereo({a1, a2}) is not None:
            continue

        nbrs1 = tuple(a for a in graph.bonded_to(a1) if a != a2)
        nbrs2 = tuple(a for a in graph.bonded_to(a2) if a != a1)
        if len(nbrs1) > len(nbrs2):
            a1, a2 = a2, a1
            nbrs1, nbrs2 = nbrs2, nbrs1

        if len(nbrs1) < 2:
            continue

        stereo1 = graph.get_atom_stereo(a1)
        stereo2 = graph.get_atom_stereo(a2)
        hb: HinderedBond | None = None

        if isinstance(stereo1, Tetrahedral) and isinstance(stereo2, Tetrahedral):
            eq_cls1 = group_by_eq(set(stereo1.atoms[1:5]) - {a2})
            eq_cls2 = group_by_eq(set(stereo2.atoms[1:5]) - {a1})

            if len(eq_cls2) > len(eq_cls1):
                a1, a2 = a2, a1
                stereo1, stereo2 = stereo2, stereo1
                eq_cls1, eq_cls2 = eq_cls2, eq_cls1

            assert stereo1.parity is not None and stereo2.parity is not None
            parity = stereo1.parity * stereo2.parity
            pattern1 = tuple(sorted((len(group) for group in eq_cls1), reverse=True))
            pattern2 = tuple(sorted((len(group) for group in eq_cls2), reverse=True))

            if (pattern1, pattern2) in {
                ((3,), (3,)),
                ((1, 1, 1), (3,)),
                ((1, 1, 1), (2, 1)),
                ((1, 1, 1), (1, 1, 1)),
            }:
                left3 = ordered_tetrahedral_neighbors(
                    stereo1,
                    a2,
                    is_left=True,
                )
                right3 = ordered_tetrahedral_neighbors(
                    stereo2,
                    a1,
                    is_left=False,
                )
                hb = HinderedBond33(atoms=(*left3, a1, a2, *right3), parity=parity)

            elif (pattern1, pattern2) == ((2, 1), (2, 1)):
                left3 = ordered_tetrahedral_neighbors(
                    stereo1,
                    a2,
                    is_left=True,
                    unique_neighbor=singled_out_neighbor(eq_cls1),
                )
                right3 = ordered_tetrahedral_neighbors(
                    stereo2,
                    a1,
                    is_left=False,
                    unique_neighbor=singled_out_neighbor(eq_cls2),
                )
                hb = HinderedBond33(atoms=(*left3, a1, a2, *right3), parity=parity)

            elif (pattern1, pattern2) == ((2, 1), (3,)):
                left3 = ordered_tetrahedral_neighbors(
                    stereo1,
                    a2,
                    is_left=True,
                    unique_neighbor=singled_out_neighbor(eq_cls1),
                )
                right3 = ordered_tetrahedral_neighbors(
                    stereo2,
                    a1,
                    is_left=False,
                )
                hb = HinderedBond33(atoms=(*left3, a1, a2, *right3), parity=parity)

        elif len(nbrs1) == 2 and isinstance(stereo2, Tetrahedral):
            if len(atom_eq_classes[nbrs1[0]]) == 2 == len(atom_eq_classes[nbrs1[1]]):
                eq_cls_t = sorted(
                    group_by_eq(set(stereo2.atoms[1:6]) - {a2}),
                    key=len,
                    reverse=True,
                )
                hb_atoms2 = next(
                    (a1, *reversed(p[2:6]))
                    for p in stereo2._perm_atoms()
                    if p[1] == a1 and eq_cls_t[0] == set(p[4:6])
                )
            else:
                hb_atoms2 = next(
                    (*reversed(p[2:6]), a1) for p in stereo2._perm_atoms() if p[1] == a1
                )
            hb = HinderedBond23(
                atoms=(*nbrs1, a1, a2, *hb_atoms2),
                parity=stereo2.parity,
            )

        if hb is None:
            continue

        hindered_bonds[frozenset({a1, a2})] = hb
        for eq_bond in eq_cls_iterator:
            for mapping in mappings:
                bond = frozenset({mapping[a1], mapping[a2]})
                if bond == eq_bond:
                    hindered_bonds[bond] = hb.__class__(
                        atoms=tuple(mapping[a] for a in hb.atoms),
                        parity=hb.parity,
                    )

    ext_sym_num = 0
    for mapping in mappings:
        if all(
            hb.__class__(atoms=tuple(mapping[a] for a in hb.atoms), parity=hb.parity)
            in hindered_bonds.values()
            for hb in hindered_bonds.values()
        ):
            ext_sym_num += 1

    return ext_sym_num
