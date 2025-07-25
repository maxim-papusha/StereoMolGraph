from __future__ import annotations

import itertools
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal, TypeVar

    from stereomolgraph.graphs.mg import (
        AtomId,
        MolGraph,
    )

def numpy_int_tuple_hash(
    arr: np.ndarray[tuple[int, ...], np.dtype[np.int64]],
    out: None|np.ndarray[tuple[Literal[1], ...], np.dtype[np.int64]] = None,
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

def label_hash(
    mg: MolGraph,
    atom_labels: Iterable[str] = ("atom_type",),
) -> dict[AtomId, int]:
    """Generates a hash for each atom based on its attributes.
    
    :param mg: MolGraph object containing the atoms.
    :param atom_labels: Iterable of attribute names to use for hashing.
    """
    atom_hash = {atom: hash(tuple(
        mg.get_atom_attribute(atom, attr ) for attr in atom_labels))
        for atom in mg.atoms
        }
    return atom_hash

def color_refine_mg(
    mg: MolGraph,
    max_iter: None|int = None,
    atom_labels: Iterable[str] = ("atom_type",),
) -> dict[AtomId, int]:
    """Color refinement algorithm for MolGraph.
    
    This algorithm refines the atom coloring based on their connectivity.
    Identical to the Weisfeiler-Lehman (1-WL) algorithm.

    :param mg: MolGraph object containing the atoms and their connectivity.
    :param max_iter: Maximum number of iterations for refinement.
        Default is None, which means it will run until convergence."""
    # np.uint32: chronological id of atom in order
    # np.int64: hash of atom label
    
    atom_label_hash = label_hash(mg, atom_labels)

    atom_hash: np.ndarray = np.array(
        [atom_label_hash[atom] for atom in mg.atoms], dtype=np.int64
    )
    
    n_atoms = np.int64(mg.n_atoms)


    id_arr = {atom: np.uint32(a_id) for a_id, atom in enumerate(mg.atoms)}
    d = {
        id_arr[atom]: np.array([id_arr[nbr] for nbr in mg.bonded_to(atom)])
        for atom in mg.atoms
    }

    grouped: defaultdict[int, dict[np.uint32, np.ndarray[tuple[int],
                                                        np.dtype[np.uint32]]]]
    grouped = defaultdict(dict)

    for atom, neighbors in d.items():
        grouped[len(neighbors)][atom] = neighbors

    # one item per group
    masks: list[np.ndarray[tuple[int], np.dtype[np.bool_]]] = []
    nbrs_id_arrs: list[np.ndarray[tuple[int, int], np.dtype[np.uint32]]] = []
    nbrs_hash_arrs: list[np.ndarray[tuple[int, int], np.dtype[np.int64]]] = []

    # groups with atoms of same number of neighbors
    for _n_nbrs, group in grouped.items():
        nbrs_id_arr = np.array([nbrs_id for nbrs_id in group.values()], dtype=np.uint32)
        nbrs_id_arrs.append(nbrs_id_arr)
        nbrs_hash_arrs.append(np.empty_like(nbrs_id_arr, dtype=np.int64))

        mask = np.zeros_like(atom_hash, dtype=np.bool_)
        
        #raise Exception(np.fromiter(group.keys(), dtype=np.int32))
        mask[[int(i) for i in group.keys()]] = True
        raise Exception(mask)
        masks.append(mask)
        assert mask.shape == atom_hash.shape
        assert nbrs_id_arr.shape == (mask.sum(), _n_nbrs)
        assert nbrs_hash_arrs[-1].shape == (mask.sum(), _n_nbrs)


    n_atom_classes: int = 0
    n = 0
    counter = itertools.repeat(None) if max_iter is None else range(max_iter)
    new_atom_hash = np.empty_like(atom_hash, dtype=np.int64)
    for _ in counter:
        n = -1
        for nbrs_id_arr, nbrs_hash_arr, mask in zip(
           nbrs_id_arrs, nbrs_hash_arrs, masks
        ):
            nbrs_hash_arr[:] = atom_hash[nbrs_id_arr]
            nbrs_hash_arr.sort(axis=-1)
            new_atom_hash[mask] = numpy_int_tuple_hash(nbrs_hash_arr, out=new_atom_hash[mask])

        new_n_classes = np.unique(new_atom_hash).shape[0]
        n = new_n_classes
        if new_n_classes == n_atom_classes:
            break
        elif new_n_classes == n_atoms:
            break
        else:
            n_atom_classes = new_n_classes
            atom_hash, new_atom_hash = new_atom_hash, atom_hash

    return {a: int(h) for a, h in zip(mg.atoms, atom_hash)}

