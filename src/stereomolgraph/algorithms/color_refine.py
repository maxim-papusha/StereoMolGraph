from __future__ import annotations

import itertools
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable, Collection, Iterator
    from typing import Literal, TypeVar

    from stereomolgraph.graphs import (
        CondensedReactionGraph,
        MolGraph,
        StereoCondensedReactionGraph,
        StereoMolGraph,
    )
    from stereomolgraph.stereodescriptors import AtomStereo, BondStereo, Stereo

    S = TypeVar("S", bound=Stereo)
    AS = TypeVar("AS", bound=AtomStereo)
    BS = TypeVar("BS", bound=BondStereo)
    N = TypeVar("N", bound=int)


def numpy_int_tuple_hash(
    arr: np.ndarray[tuple[int, ...], np.dtype[np.int64]],
    out: None | np.ndarray[tuple[Literal[1], ...], np.dtype[np.int64]] = None,
) -> np.ndarray:
    """
    Mimics the python SipHash hashing function for tuples of integers
    with numpy int64 arrays.

    def SipHash(arr_slice):
        h = 0x345678
        mult = 1000003
        length = len(arr_slice)
        for i in range(1, length + 1):
            h = (h ^ arr_slice[i-1]) * mult
            mult += 82520 + 2 * (length - i)
        return h + 97531
    """
    # overflow is an expected behavior in this case
    with np.errstate(over="ignore"):
        arr_shape = arr.shape
        length = arr_shape[-1]
        if out is None:
            output = np.full(arr_shape[:-1], 0x345678, dtype=np.int64)
        else:
            output = out
            output.fill(0x345678)

        n = 82518 + 2 * length
        m = range(n, n - 2 * (length - 1), -2)
        mults = itertools.accumulate(m, initial=1000003)

        for idx, mult in enumerate(mults):
            output ^= arr[..., idx]
            output *= mult

        output += 97531
        return output


def numpy_int_multiset_hash(
    arr: np.ndarray[tuple[int, ...], np.dtype[np.int64]],
    out: None | np.ndarray[tuple[Literal[1], ...], np.dtype[np.int64]] = None,
) -> np.ndarray:
    """
    Hash function for a multiset (order-independent with duplicates) of integers.
    Works by sorting the elements and then applying the tuple hashing function.
    """
    sorted_arr = np.sort(arr, axis=-1)
    return numpy_int_tuple_hash(sorted_arr, out)


def label_hash(
    mg: MolGraph,
    atom_labels: Collection[str] = ("atom_type",),
) -> np.ndarray[tuple[int], np.dtype[np.int64]]:
    """Generates a hash for each atom based on choosen attributes.

    :param mg: MolGraph object containing the atoms.
    :param atom_labels: Iterable of attribute names to use for hashing.
    """
    if len(atom_labels) == 1:
        atom_label = next(iter(atom_labels))
        atom_hash = [
            hash(mg.get_atom_attribute(atom, atom_label)) for atom in mg.atoms
        ]
    else:
        atom_hash = [
            hash(
                frozenset(
                    (attr, mg.get_atom_attribute(atom, attr))
                    for attr in atom_labels
                )
            )
            for atom in mg.atoms
        ]
    return np.array(atom_hash, dtype=np.int64)


def morgan_generator(
    mg: MolGraph,
    atom_labels: None | np.ndarray[tuple[int], np.dtype[np.int64]] = None,
) -> Iterator[np.ndarray[tuple[int], np.dtype[np.int64]]]:
    """Color refinement algorithm for MolGraph.

    This algorithm refines the atom coloring based on their connectivity.
    Identical to the Weisfeiler-Lehman (1-WL) algorithm.

    :param mg: MolGraph object containing the atoms and their connectivity.
    :param max_iter: Maximum number of iterations for refinement.
        Default is None, which means it will run until convergence."""
    n_atoms = len(mg.atoms)
    if atom_labels is not None:
        assert len(atom_labels) == n_atoms

    atom_hash = (
        label_hash(mg, ("atom_type",)) if atom_labels is None else atom_labels
    )
    yield atom_hash
    if n_atoms == 0:
        return

    arr_id_dict, id_arr_dict = {}, {}
    for id, atom in enumerate(mg.atoms):
        arr_id_dict[atom] = id
        id_arr_dict[id] = atom

    bonded_lst = [
        (id, [arr_id_dict[a] for a in mg.bonded_to(atom)])
        for atom, id in arr_id_dict.items()
    ]

    bonded_lst.sort(key=lambda x: len(x[1]))

    id_nbrs_tuple_list = []
    for n_nbrs, group in itertools.groupby(
        bonded_lst, key=lambda x: len(x[1])
    ):
        if n_nbrs == 0:
            continue  # Skip aggregation if no neighbors

        group = list(group)
        ids_lst, nbrs_lists = [], []
        for id, nbrs in group:
            ids_lst.append(id)
            nbrs_lists.append(nbrs)

        ids = np.array(ids_lst, dtype=np.int16)
        nbrs = np.array(nbrs_lists, dtype=np.int16)
        id_nbrs_tuple_list.append((ids, nbrs))

    for _ in itertools.repeat(None):
        for ids, nbrs in id_nbrs_tuple_list:
            # Compute the new hash for each atom based on its neighbors
            atom_hash[ids] = numpy_int_multiset_hash(atom_hash[nbrs])

        yield atom_hash


def stereo_morgan_generator(
    smg: StereoMolGraph,
    atom_labels: None | np.ndarray[tuple[int], np.dtype[np.int64]] = None,
) -> Iterator[np.ndarray[tuple[int], np.dtype[np.int64]]]:
    n_atoms = len(smg.atoms)
    if atom_labels is not None:
        assert len(atom_labels) == n_atoms

    atom_hash = (
        label_hash(smg, ("atom_type",)) if atom_labels is None else atom_labels
    )

    yield atom_hash

    if (n_atoms == 0
        or n_atoms == 1
        or all(len(c) == 1 for c in smg.connected_components())
    ):
        while True:
            yield atom_hash

    prev_atom_hash = np.zeros_like(atom_hash, dtype=np.int64)

    arr_id_dict, id_arr_dict = {}, {}
    stereo_hash_pointer = {}
    for id, atom in enumerate(smg.atoms):
        arr_id_dict[atom] = id
        id_arr_dict[id] = atom
        stereo_hash_pointer[id] = []  # arr_id: list[stereo_pointer]

    grouped_atom_stereo: dict = defaultdict(list)
    atoms_with_atom_stereo: set[int] = set()
    atoms_with_bond_stereo: set[int] = set()

    grouped_bond_stereo: dict = defaultdict(list)
    atoms_with_atom_stereo: set[int] = set()

    as_atoms = []
    as_perm_atoms = []
    as_nbr_atoms = []

    bs_atoms = []
    bs_nbr_atoms = []
    bs_perm_atoms = []

    # i: arrays to store intermediate values.
    # Avoids additional memory allocation.
    i_a_perm_nbrs = []
    i_a_perm = []
    i_a = []

    i_b_perm_nbrs = []
    i_b_perm = []
    i_b = []

    for atom, stereo in smg.atom_stereo.items():
        if stereo.parity is not None:
            nbr_atoms = (
                stereo.atoms
                if stereo.parity != -1
                else stereo._inverted_atoms()  # type: ignore
            )
            grouped_atom_stereo[stereo.PERMUTATION_GROUP].append(
                (atom, nbr_atoms)
            )

            atoms_with_atom_stereo.add(atom)

    for bond, stereo in smg.bond_stereo.items():
        if stereo.parity is not None:
            nbr_atoms = (
                stereo.atoms
                if stereo.parity != -1
                else stereo._inverted_atoms()  # type: ignore
            )
            grouped_bond_stereo[stereo.PERMUTATION_GROUP].append(
                (bond, nbr_atoms)
            )

            for atom in bond:
                atoms_with_bond_stereo.add(atom)

    atoms_without_atom_stereo = set(smg.atoms) - atoms_with_atom_stereo

    for atom in atoms_without_atom_stereo:
        fake_stereo_atoms = (atom, *smg.bonded_to(atom))
        perm_gen = itertools.permutations(range(1, len(fake_stereo_atoms)))
        perm_group = tuple((0, *perm) for perm in perm_gen)
        grouped_atom_stereo[perm_group].append((atom, fake_stereo_atoms))

    # atom_stereo
    for perm_group, atom_nbr_atoms_list_tup in grouped_atom_stereo.items():
        atom_arr_ids = np.array(
            [arr_id_dict[atom] for atom, _ in atom_nbr_atoms_list_tup],
            dtype=np.uint16,
        )
        as_atoms.append(atom_arr_ids)

        nbr_atoms = np.array(
            [
                [arr_id_dict[a] for a in nbr_lst]
                for _atom, nbr_lst in atom_nbr_atoms_list_tup
            ],
            dtype=np.uint16,
        )

        as_nbr_atoms.append(nbr_atoms)

        perm_group = np.array(perm_group, dtype=np.uint8)
        perm_atoms = nbr_atoms[..., perm_group]
        as_perm_atoms.append(perm_atoms)

        # intermediate arrays
        a_perm_nbrs = np.zeros(perm_atoms.shape, dtype=np.int64)
        i_a_perm_nbrs.append(a_perm_nbrs)
        i_a_perm.append(np.zeros(a_perm_nbrs.shape[0:2], dtype=np.int64))
        i_atom_stereo = np.zeros(a_perm_nbrs.shape[0:1], dtype=np.int64)
        i_a.append(i_atom_stereo)

        for stereo_id, atom_arr_id in enumerate(atom_arr_ids):
            stereo_hash_pointer[atom_arr_id].append(
                i_atom_stereo[stereo_id : stereo_id + 1]
            )
            # by reference

    # bond_stereo
    for perm_group, atom_nbr_atoms_list_tup in grouped_bond_stereo.items():
        atom_arr_ids = np.array(
            [
                [arr_id_dict[atom] for atom in bond]
                for bond, _ in atom_nbr_atoms_list_tup
            ],
            dtype=np.uint16,
        )
        bs_atoms.append(atom_arr_ids)

        nbr_atoms = np.array(
            [
                [arr_id_dict[a] for a in nbr_lst]
                for _atom, nbr_lst in atom_nbr_atoms_list_tup
            ],
            dtype=np.uint16,
        )

        bs_nbr_atoms.append(nbr_atoms)
        arr_perm_group = np.array(perm_group, dtype=np.uint8)
        perm_atoms = nbr_atoms[..., arr_perm_group]
        bs_perm_atoms.append(perm_atoms)

        # intermediate arrays
        b_perm_nbrs = np.zeros(perm_atoms.shape, dtype=np.int64)
        i_b_perm_nbrs.append(b_perm_nbrs)
        i_b_perm.append(np.zeros(b_perm_nbrs.shape[0:2], dtype=np.int64))
        i_bond_stereo = np.zeros(b_perm_nbrs.shape[0:1], dtype=np.int64)
        i_bond_stereo.fill(hash(perm_group))
        i_b.append(i_bond_stereo)

        for stereo_id, (atom_arr_id1, atom_arr_id2) in enumerate(atom_arr_ids):
            stereo_hash_pointer[atom_arr_id1].append(
                i_bond_stereo[stereo_id : stereo_id + 1]
            )
            stereo_hash_pointer[atom_arr_id2].append(
                i_bond_stereo[stereo_id : stereo_id + 1]
            )
        # by reference

    i_atoms_with_n_stereo = []  # atoms, i_stereo, group

    pntr = sorted(
        ((id, ptr) for id, ptr in stereo_hash_pointer.items() if ptr),
        key=lambda x: len(x[1]),
    )

    for key, group in itertools.groupby(pntr, key=lambda x: len(x[1])):
        group = list(group)

        ids = []
        pntr_groups = []
        for id, ptrs in group:
            ids.append(id)
            pntr_groups.append(ptrs)

        atoms = np.array(ids, dtype=np.int16)
        i_hash = np.empty((len(pntr_groups), key), dtype=np.int64)
        i_atoms_with_n_stereo.append((atoms, i_hash, pntr_groups))

    for count in itertools.count(1, 1):
        # atom stereo
        for perm_atoms, a_perm_nbrs, a_perm, a in zip(
            as_perm_atoms, i_a_perm_nbrs, i_a_perm, i_a
        ):
            numpy_int_tuple_hash(atom_hash[perm_atoms], out=a_perm)
            numpy_int_multiset_hash(a_perm, out=a)

        # bond stereo
        if count != 0:
            for perm_atoms, b_perm_nbrs, b_perm, b in zip(
                bs_perm_atoms, i_b_perm_nbrs, i_b_perm, i_b
            ):
                numpy_int_tuple_hash(prev_atom_hash[perm_atoms], out=b_perm)
                numpy_int_multiset_hash(b_perm, out=b)

        for (
            atoms,
            i_stereo,
            ptr_lsts,
        ) in i_atoms_with_n_stereo:  # atoms, i_stereo, group
            prev_atom_hash[:] = atom_hash[:]

            i_stereo[:] = np.asarray(
                [[ptr.item() for ptr in ptr_list] for ptr_list in ptr_lsts]
            )
            atom_hash[atoms] = numpy_int_multiset_hash(i_stereo)

        yield atom_hash


def _reaction_generator(
    graph: CondensedReactionGraph,
    generator: Callable,
    max_iter: int | None = None,
    atom_labels: None | np.ndarray[tuple[int], np.dtype[np.int64]] = None,
) -> Iterator[np.ndarray[tuple[int], np.dtype[np.int64]]]:
    r_colors = generator(
        graph.reactant(), atom_labels=atom_labels
    )
    p_colors = generator(
        graph.product(),  atom_labels=atom_labels
    )
    ts_colors = generator(graph._ts(), atom_labels=atom_labels)

    for _ in itertools.repeat(None) if max_iter is None else range(max_iter):
        stacked = np.stack(
            [next(c) for c in (r_colors, p_colors, ts_colors)],
            axis=-1,
            dtype=np.int64,
        )

        yield numpy_int_tuple_hash(stacked)


def reaction_morgan_generator(
    graph: CondensedReactionGraph,
    max_iter: int | None = None,
    atom_labels: None | np.ndarray[tuple[int], np.dtype[np.int64]] = None,
) -> Iterator[np.ndarray[tuple[int], np.dtype[np.int64]]]:
    return _reaction_generator(
        graph=graph,
        generator=morgan_generator,
        max_iter=max_iter,
        atom_labels=atom_labels,
    )


def stereo_reaction_morgan_generator(
    graph: StereoCondensedReactionGraph,
    max_iter: int | None = None,
    atom_labels: None | np.ndarray[tuple[int], np.dtype[np.int64]] = None,
) -> Iterator[np.ndarray[tuple[int], np.dtype[np.int64]]]:
    return _reaction_generator(
        graph=graph,
        generator=stereo_morgan_generator,
        max_iter=max_iter,
        atom_labels=atom_labels,
    )


def _color_refine(
    graph: MolGraph,
    generator: Callable,
    max_iter: int | None = None,
    atom_labels: None | np.ndarray[tuple[int], np.dtype[np.int64]] = None,
) -> np.ndarray[tuple[int], np.dtype[np.int64]]:
    sm_generator = generator(graph, atom_labels=atom_labels)

    n_atoms = graph.n_atoms

    atom_hash = next(sm_generator)
    n_atom_classes = np.unique(atom_hash).shape[0]

    counter = (
        itertools.count(1, 1) if max_iter is None else range(max_iter + 1)
    )
    for _ in counter:
        atom_hash = next(sm_generator)
        new_n_classes = np.unique(atom_hash).shape[0]
        if new_n_classes == n_atom_classes:
            break
        elif new_n_classes == n_atoms:
            break
        else:
            n_atom_classes = new_n_classes

    return atom_hash


def color_refine_mg(
    graph: MolGraph,
    max_iter: int | None = None,
    atom_labels: None | np.ndarray[tuple[int], np.dtype[np.int64]] = None,
) -> np.ndarray[tuple[int], np.dtype[np.int64]]:
    return _color_refine(
        graph=graph,
        generator=morgan_generator,
        max_iter=max_iter,
        atom_labels=atom_labels,
    )


def color_refine_smg(
    graph: StereoMolGraph,
    max_iter: int | None = None,
    atom_labels: None | np.ndarray[tuple[int], np.dtype[np.int64]] = None,
) -> np.ndarray[tuple[int], np.dtype[np.int64]]:
    return _color_refine(
        graph=graph,
        generator=stereo_morgan_generator,
        max_iter=max_iter,
        atom_labels=atom_labels,
    )


def color_refine_crg(
    graph: CondensedReactionGraph,
    max_iter: int | None = None,
    atom_labels: None | np.ndarray[tuple[int], np.dtype[np.int64]] = None,
) -> np.ndarray[tuple[int], np.dtype[np.int64]]:
    return _color_refine(
        graph=graph,
        generator=reaction_morgan_generator,
        max_iter=max_iter,
        atom_labels=atom_labels,
    )


def color_refine_scrg(
    graph: StereoCondensedReactionGraph,
    max_iter: int | None = None,
    atom_labels: None | np.ndarray[tuple[int], np.dtype[np.int64]] = None,
) -> np.ndarray[tuple[int], np.dtype[np.int64]]:
    return _color_refine(
        graph=graph,
        generator=stereo_reaction_morgan_generator,
        max_iter=max_iter,
        atom_labels=atom_labels,
    )
