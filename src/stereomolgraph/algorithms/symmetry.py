from __future__ import annotations

from collections import deque
from collections.abc import Callable, Hashable, Iterable, Mapping
from typing import TypeVar

from stereomolgraph import AtomId, Bond, StereoMolGraph
from stereomolgraph.algorithms.circular import color_refine_smg
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.stereodescriptors import (
    OInt,
    PlanarBond,
    RigidBond,
    RigidBond12,
    RigidBond13,
    RigidBond23,
    RigidBond33,
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


def symmetry_number(graph: StereoMolGraph, atom_labels=None) -> int:
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
    # most efficient way to get the len of a generator
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
        s = RigidBond33(atoms=(*nbrs1, a1, a2, *nbrs2), parity=1)
    elif len(nbrs1) == 2 and len(nbrs2) == 3:
        s = RigidBond23(atoms=(*nbrs1, a1, a2, *nbrs2), parity=1)
    elif len(nbrs1) == 1 and len(nbrs2) == 3:
        s = RigidBond13(atoms=(*nbrs1, a1, a2, *nbrs2), parity=1)
    elif len(nbrs1) == 2 and len(nbrs2) == 2:
        s = PlanarBond(atoms=(*nbrs1, a1, a2, *nbrs2), parity=0)
    elif len(nbrs1) == 1 and len(nbrs2) == 2:
        s = RigidBond12(atoms=(*nbrs1, a1, a2, *nbrs2), parity=0)
    else:
        raise NotImplementedError(
            "Bonds with more than 3 substituents are not supported"
        )
    unique_reorderings: set[RigidBond] = set()
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


def external_symmetry_number(
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

    bond_eq_classes = bond_automorphism_classes(graph, mappings=mappings)

    def ordered_tetrahedral_neighbors(
        stereo: Tetrahedral,
        bonded_atom: AtomId,
        *,
        is_left: bool,
    ) -> tuple[OInt, OInt, OInt]:
        for permuted_atoms in stereo._perm_atoms():
            if permuted_atoms[1] != bonded_atom:
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

    rigid_bonds: dict[Bond, RigidBond] = {}

    # 1. Loop over bond equivalence classes
    for eq_cls in {frozenset(eq_cls) for eq_cls in bond_eq_classes.values()}:
        eq_bonds = tuple(eq_cls)
        eq_cls_iterator = iter(eq_bonds)

        # a. Atoms a1, a2 in the first bond of the equivalence class
        a1, a2 = next(eq_cls_iterator)

        # b. If there is already bond stereochemistry, skip this equivalence class
        if graph.get_bond_stereo({a1, a2}) is not None:
            continue

        # c. Determine the neighbors external to the bond
        nbrs1 = tuple(a for a in graph.bonded_to(a1) if a != a2)
        nbrs2 = tuple(a for a in graph.bonded_to(a2) if a != a1)
        if len(nbrs1) > len(nbrs2):
            a1, a2 = a2, a1
            nbrs1, nbrs2 = nbrs2, nbrs1

        if len(nbrs1) < 2:
            continue

        # d. Get the stereochemistry of the atoms
        stereo1 = graph.get_atom_stereo(a1)
        stereo2 = graph.get_atom_stereo(a2)
        hb: RigidBond | None = None

        # e. Assign RigidBond33 for two tetrahedral atoms
        if isinstance(stereo1, Tetrahedral) and isinstance(stereo2, Tetrahedral):
            assert stereo1.parity is not None and stereo2.parity is not None
            parity = stereo1.parity * stereo2.parity

            # i. Form RigidBond33, ordering external neighbors consistently
            #    with the tetrahedral stereochemistry
            left3 = ordered_tetrahedral_neighbors(stereo1, a2, is_left=True)
            right3 = ordered_tetrahedral_neighbors(stereo2, a1, is_left=False)
            hb = RigidBond33(atoms=(*left3, a1, a2, *right3), parity=parity)

        # f. Assign RigidBond23 for an sp2 and a tetrahedral atom
        elif len(nbrs1) == 2 and isinstance(stereo2, Tetrahedral):
            # i. Get the right neighbors in an order consistent with their parity
            # (Reverse the order of the external neighbors since the parity is opposite)
            right5 = next(
                (a1, a2, *reversed(p[2:])) for p in stereo2._perm_atoms() if p[1] == a1
            )
            # i. If both external neighbors of a1 have orbits of length 2
            hb = RigidBond23(atoms=(*nbrs1, *right5), parity=stereo2.parity)

        if hb is None:
            continue

        rigid_bonds[frozenset({a1, a2})] = hb
        for eq_bond in eq_cls_iterator:
            for mapping in mappings:
                bond = frozenset({mapping[a1], mapping[a2]})
                if bond == eq_bond:
                    rigid_bonds[bond] = hb.__class__(
                        atoms=tuple(mapping[a] for a in hb.atoms),
                        parity=hb.parity,
                    )

    ext_sym_num = 0
    for mapping in mappings:
        if all(
            hb.__class__(atoms=tuple(mapping[a] for a in hb.atoms), parity=hb.parity)
            in rigid_bonds.values()
            for hb in rigid_bonds.values()
        ):
            ext_sym_num += 1
    if ext_sym_num == 0:
        raise RuntimeError("External symmetry number is zero, this should not happen")
    return ext_sym_num
