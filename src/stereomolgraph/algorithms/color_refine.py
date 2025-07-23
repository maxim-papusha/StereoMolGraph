from __future__ import annotations

import itertools
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal, Optional, TypeVar

    from stereomolgraph.graphs.mg import (
        AtomId,
        MolGraph,
    )

    N = TypeVar("N", bound=int)

def numpy_int_tuple_hash(
    arr: np.ndarray[tuple[int, ...], np.dtype[np.int64]],
    out: Optional[np.ndarray[tuple[Literal[1], ...], np.dtype[np.int64]]] = None,
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

def numpy_int_set_hash(): ...

def label_hash(
    mg: MolGraph,
    atom_labels: Optional[Iterable[str]] = ("atom_type",),
    bond_labels: Optional[Iterable[str]] = None,
) -> dict[AtomId, int]:
    if atom_labels == ("atom_type",) and bond_labels is None:
        atom_hash = {
            atom: hash(mg.get_atom_type(atom))
            for atom in mg.atoms
        }

    elif atom_labels is None and bond_labels is None:
        atom_hash = {atom: 0 for atom in mg.atoms}

    elif atom_labels:
        atom_labels = sorted(atom_labels)
        atom_labels.append("atom_type")
        bond_labels = sorted(bond_labels) if bond_labels else []
        bond_labels.append("reaction")
        atom_hash = {atom:
              hash((
            tuple([(atom_label, label_dict.get(atom_label, None))
                   for atom_label in atom_labels]),
            tuple(sorted([(tuple(sorted(
                mg.get_bond_attributes(atom, nbr, bond_labels).items()))
                           for nbr in mg.bonded_to(atom))])),
            ))
              for atom, label_dict in mg.atoms_with_attributes.items()
        }
    else:
        raise ValueError("Invalid combination of atom and bond labels.")
    return atom_hash

def color_refine_mg(
    mg: MolGraph,
    max_iter: Optional[int] = None,
    atom_labels: Optional[Iterable[str]] = ("atom_type",),
    bond_labels: Optional[Iterable[str]] = None,
) -> dict[AtomId, int]:
    atom_label_hash = label_hash(mg, atom_labels, bond_labels)

    atom_hash: np.ndarray = np.array(
        [atom_label_hash[atom] for atom in mg.atoms], dtype=np.int64
    )

    n_atoms = np.int64(mg.n_atoms)
    id_arr = {atom: a_id for a_id, atom in enumerate(mg.atoms)}
    d = {
        id_arr[atom]: {id_arr[nbr] for nbr in mg.bonded_to(atom)}
        for atom in mg.atoms
    }

    grouped: defaultdict[int, dict[int, set[int]]] = defaultdict(dict)
    for key, value in d.items():
        grouped[len(value)][key] = value

    masks: list[np.ndarray] = []
    data: list[np.ndarray] = []
    t_arrs: list[np.ndarray] = []
    t_hashs: list[np.ndarray] = []
    a_hashs: list[np.ndarray] = []

    for group in list(grouped.values()):
        mask = np.zeros_like(atom_hash, dtype=np.bool_)
        k = [int(i) for i in group.keys()]
        mask[k] = True
        group_values = [(k, *v) for k, v in group.items()]  # rename me

        n_neigh = len(group_values[0])
        perm = itertools.permutations(range(1, n_neigh))

        perm_with_zero = [(0,) + p for p in perm]

        g = np.array(group_values, dtype=np.int64)[:, perm_with_zero]

        masks.append(mask)
        data.append(g)
        t_arrs.append(np.empty_like(g, dtype=np.int64))
        t_hashs.append(np.empty(shape=t_arrs[-1].shape[:-1], dtype=np.int64))
        a_hashs.append(np.empty(shape=t_hashs[-1].shape[:-1], dtype=np.int64))

    n_atom_classes = None
    counter = itertools.repeat(None) if max_iter is None else range(max_iter)
    new_atom_hash = np.empty_like(atom_hash, dtype=np.int64)
    
    for _ in counter:
        for d, m, t_arr, t_hash, a_hash in zip(
            data, masks, t_arrs, t_hashs, a_hashs
        ):
            t_arr[:] = atom_hash[d]
            t_hash = numpy_int_tuple_hash(t_arr, out=t_hash)
            t_hash.sort(axis=-1) # defaults to quicksort
            a_hash = numpy_int_tuple_hash(t_hash, out=a_hash)
            new_atom_hash[m] = a_hash

        new_n_classes = np.unique(new_atom_hash).shape[0]

        if new_n_classes == n_atom_classes:
            break
        elif new_n_classes == n_atoms:
            break
        else:
            n_atom_classes = new_n_classes
            atom_hash, new_atom_hash = new_atom_hash, atom_hash

    return {a: int(h) for a, h in zip(mg.atoms, atom_hash)}

