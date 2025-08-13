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
        StereoMolGraph,
    )
    from stereomolgraph.stereodescriptors import Stereo, AtomStereo

    S = TypeVar("S", bound=Stereo)
    AS = TypeVar("AS", bound=AtomStereo)
    BS = TypeVar("BS", bound=BondStereo)
    N = TypeVar("N", bound=int)

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

def numpy_int_multiset_hash(
    arr: np.ndarray[tuple[int, ...], np.dtype[np.int64]],
    out: None|np.ndarray[tuple[Literal[1], ...], np.dtype[np.int64]] = None,
) -> np.ndarray:
    """
    Hash function for a multiset (order-independent with duplicates) of integers.
    Works by sorting the elements and then applying the tuple hashing function.
    """
    sorted_arr = np.sort(arr, axis=-1)
    return numpy_int_tuple_hash(sorted_arr, out)

def label_hash(
    mg: MolGraph,
    atom_labels: Iterable[str] = ("atom_type",),
) -> dict[AtomId, int]:
    """Generates a hash for each atom based on its attributes.
    
    :param mg: MolGraph object containing the atoms.
    :param atom_labels: Iterable of attribute names to use for hashing.
    """
    if len(atom_labels) == 1:
        atom_label = atom_labels[0]
        atom_hash = {atom: hash(mg.get_atom_attribute(atom, atom_label))
        for atom in mg.atoms}
    else:
        atom_hash = {atom: hash(frozenset(
        (attr, mg.get_atom_attribute(atom, attr)) for attr in atom_labels))
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
    
    atom_label_hash = label_hash(mg, atom_labels)

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

def color_refine_smg(
    smg: StereoMolGraph,
    max_iter: None|int = None,
    atom_labels: Iterable[str] = ("atom_type",),
) -> dict[AtomId, int]:

    atom_label_hash = label_hash(smg, atom_labels)
    atom_hash = np.array(
        [atom_label_hash[atom] for atom in mg.atoms], dtype=np.int64
    )

    arr_id = {atom: a_id for a_id, atom in enumerate(mg.atoms)}

    grouped_atom_stereo: dict[type[S], list[tuple[int, ...]]] = defaultdict(list)
    atoms_with_atom_stereo: set[int] = set()

    s_perm_groups = []
    s_atoms = []
    s_nbr_atoms = []

    s_nbrs_hash_i = []
    s_nbrs_hash_perm_i = []

    for atom, stereo in smg.get_atom_stereo.items():
        if stereo.parity is not None:
            atoms_with_atom_stereo.add(atom)
            nbr_atoms = stereo.atoms if parity != -1 else stereo.inverted_atoms
            grouped_atom_stereo[stereo.__class__.PERMUTATION_GROUP].append(nbr_atoms)

    
    atoms_with_atom_stereo = set(smg.atoms) - atoms_with_atom_stereo
    # TODO: for atoms without atom_stereo

    for perm_group, nbr_atoms_list in grouped_atom_stereo.items():
        
        s_perm_groups.append(np.array(perm_group, dtype=np.uint8))

        s_atoms.append(np.fromiter((s.central_atom for s in stereo_list),
                        dtype=np.uint16))

        nbr_atoms = np.array([[arr_id[a] for a in 
                        (s.atoms if s.parity != -1 else s._inverted_atoms)]  
                     for s in stereo_list ], dtype=np.uint16)

        s_nbr_atoms.append(nbr_atoms)

        nbrs_hash = np.empty_like(nbr_atoms, dtype=np.int64)
        s_nbrs_hash_i.append(nbrs_hash)

        nbrs_hash_perm = np.empty( (*(nbrs_hash.shape), len(perm_group)),
                                 dtype=np.int64)
        s_nbrs_hash_perm_i.append(nbrs_hash_perm)