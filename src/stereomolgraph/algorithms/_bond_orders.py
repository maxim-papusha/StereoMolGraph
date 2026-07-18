"""Bond-order assignment from atom connectivity.

Derived from xyz2mol (https://github.com/jensengroup/xyz2mol) by
Jan H. Jensen, University of Copenhagen, implementing the algorithm from

    Yeonjoon Kim and Woo Youn Kim
    "Universal Structure Conversion Method for Organic Molecules:
    From Atomic Connectivity to Three-Dimensional Geometry"
    Bull. Korean Chem. Soc. 2015, Vol. 36, 1769-1777
    DOI: 10.1002/bkcs.10334

Adapted to use ``MolGraph`` data structures and dict-based bond
representation instead of numpy adjacency matrices.

License: MIT — see end of file.
"""

from __future__ import annotations

import itertools
import warnings
from collections import defaultdict
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

    from stereomolgraph.graphs.mg import AtomId, Bond, MolGraph


@dataclass(slots=True)
class BondOrders:
    """Result of bond order assignment for a molecular graph.

    Attributes
    ----------
    bond_order:
        Bond order for each bond, keyed by ``Bond`` (``frozenset[AtomId, AtomId]``).
        Only bonds with order >= 1 are included.
    charges:
        Formal charge per atom, keyed by ``AtomId``. Defaults to 0.
    unpaired_electrons:
        Number of unpaired (radical) electrons per atom, keyed by ``AtomId``.
        Defaults to 0.
    """

    bond_order: dict[Bond, int]
    charges: dict[AtomId, int]
    unpaired_electrons: dict[AtomId, int]


atomic_valence: dict[int, list[int]] = defaultdict(
    list,
    {
        1: [1],
        5: [3, 4],
        6: [4, 2],
        7: [3, 4],
        8: [2, 1, 3],
        9: [1],
        13: [3, 4],
        14: [4],
        15: [3, 5],
        16: [2, 4, 6],
        17: [1],
        18: [0],
        32: [4],
        33: [5, 3],
        34: [2],
        35: [1],
        52: [2],
        53: [1],
        **{z: [20] for z in range(21, 31)},
        **{z: [20] for z in range(39, 49)},
        **{z: [20] for z in range(57, 81)},
        **{z: [20] for z in range(89, 104)},
    },
)

atomic_valence_electrons: dict[int, int] = {
    1: 1,
    5: 3,
    6: 4,
    7: 5,
    8: 6,
    9: 7,
    13: 3,
    14: 4,
    15: 5,
    16: 6,
    17: 7,
    18: 8,
    32: 4,
    33: 5,
    34: 6,
    35: 7,
    52: 6,
    53: 7,
    **{z: z - 18 for z in range(21, 31)},
    **{z: z - 36 for z in range(39, 49)},
    57: 3,
    **{z: z - 68 for z in range(71, 81)},
    **{z: 3 for z in range(58, 71)},
    **{89: 3, 90: 4, 91: 5, 92: 6, 93: 7, 94: 8},
    **{z: 3 for z in range(95, 104)},
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def connectivity2bond_orders(
    graph: MolGraph,
    allow_charged_fragments: bool = False,
    charge: int = 0,
) -> BondOrders:
    """Calculates Bond orders from atom connectivity.

    Bond orders can be assigned automatically using the algorithm from

    [Yeonjoon Kim and Woo Youn Kim "Universal Structure Conversion Method for
    Organic Molecules: From Atomic Connectivity to Three-Dimensional Geometry"
    Bull. Korean Chem. Soc. 2015, Vol. 36, 1769-1777](https://doi.org/10.1002/bkcs.10334)

    :param graph: Molecular graph whose bonds, atom types, and connectivity
                  are used as input.
    :param allow_charged_fragments: If false radicals are formed and if True
                                ions are preferred.
    :param charge: charge of the whole molecule, defaults to 0.
    :return: :class:`BondOrders` with bonds, charges, and unpaired electrons.
    """
    atom_ids = list(graph.atoms)
    num_atoms = len(atom_ids)
    atomic_nums = [int(graph.get_atom_type(aid)) for aid in atom_ids]

    charges: dict[int, int] = {}
    unpaired_electrons: dict[int, int] = {}

    # Pre-compute adjacency lookup (positional indices 0..num_atoms-1)
    aid_to_idx = {aid: i for i, aid in enumerate(atom_ids)}
    adj_lookup = [[False] * num_atoms for _ in range(num_atoms)]
    for bond in graph.bonds:
        a, b = bond
        i, j = aid_to_idx[a], aid_to_idx[b]
        adj_lookup[i][j] = adj_lookup[j][i] = True

    # Connectivity valence per atom
    conn_valence = [len(graph.neighbors[aid]) for aid in atom_ids]

    # Run the bond order assignment
    bond_orders = _AC2BO(
        atom_ids=atom_ids,
        atom_nrs=atomic_nums,
        charge=charge,
        adj_lookup=adj_lookup,
        conn_valence=conn_valence,
        allow_charged_fragments=allow_charged_fragments,
        allow_carbenes=False,
    )

    # Per-atom bond order totals
    bo_total = [0] * num_atoms
    for i, aid in enumerate(atom_ids):
        bo_total[i] = sum(
            bond_orders.get(frozenset({aid, atom_ids[j]}), 0)
            for j in range(num_atoms)
            if adj_lookup[i][j]
        )

    # Set atomic charges
    if allow_charged_fragments:
        mol_charge = 0
        for i, atomic_num in enumerate(atomic_nums):
            formal_charge = _get_atomic_charge(
                atomic_num,
                atomic_valence_electrons=atomic_valence_electrons[atomic_num],
                bo_total=bo_total[i],
            )
            mol_charge += formal_charge
            if atomic_num == 6:
                n_single = sum(
                    1
                    for j in range(num_atoms)
                    if adj_lookup[i][j]
                    and bond_orders.get(frozenset({atom_ids[i], atom_ids[j]}), 0) == 1
                )
                if n_single == 2 and bo_total[i] == 2:
                    mol_charge -= formal_charge
                    formal_charge = 0
                    unpaired_electrons[atom_ids[i]] = 2
                if n_single == 3 and mol_charge + 1 < charge:
                    mol_charge += 2
                    formal_charge = 1
            if formal_charge != 0:
                charges[atom_ids[i]] = formal_charge

    # Set atomic radicals
    for i, atomic_num in enumerate(atomic_nums):
        if unpaired_electrons.get(atom_ids[i], 0) > 0:
            continue
        formal_charge = _get_atomic_charge(
            atomic_num,
            atomic_valence_electrons=atomic_valence_electrons[atomic_num],
            bo_total=bo_total[i],
        )
        if abs(formal_charge) > 0:
            unpaired_electrons[atom_ids[i]] = abs(int(formal_charge))

    return BondOrders(
        bond_order={b: o for b, o in bond_orders.items() if o > 0},
        charges=charges,
        unpaired_electrons=unpaired_electrons,
    )


# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------


def _AC2BO(
    atom_ids: list[int],
    atom_nrs: list[int],
    charge: int,
    adj_lookup: list[list[bool]],
    conn_valence: list[int],
    allow_charged_fragments: bool = True,
    allow_carbenes: bool = False,
) -> dict[frozenset[int], int]:
    """Assign bond orders to a molecular graph.

    Returns a dict mapping ``frozenset[AtomId, AtomId]`` → bond order.
    """
    num_atoms = len(atom_ids)

    # Initial bond orders: every existing bond starts at 1
    bond_orders: dict[frozenset[int], int] = {}
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            if adj_lookup[i][j]:
                bond_orders[frozenset({atom_ids[i], atom_ids[j]})] = 1

    # Build allowed valences per atom (must be >= connectivity valence)
    possible_valences: list[list[int]] = []
    for i, (atomic_num, valence) in enumerate(zip(atom_nrs, conn_valence)):
        allowed = [x for x in atomic_valence[atomic_num] if x >= valence]
        if atomic_num == 6 and valence == 1:
            allowed = [x for x in allowed if x != 2]
        if atomic_num == 6 and not allow_carbenes and valence == 2:
            allowed = [x for x in allowed if x != 2]
        if atomic_num == 6 and valence == 2 and 3 not in allowed:
            allowed.append(3)
        if atomic_num == 16 and valence == 1:
            allowed = [1, 2]

        if not allowed:
            warnings.warn(
                f"Valence of atom {i} is {valence}, which exceeds allowed "
                f"maximum {max(atomic_valence[atomic_num])}. Continuing"
            )
        possible_valences.append(allowed)

    # Sort valence combos: prefer lower valences for O first, then N, C, P, S.
    # Lexicographic tuple ordering naturally gives this priority.
    # O valence 2 is the most stable (carbonyl, hydroxyl, ether) so it scores
    # lowest; O valence 3 (oxonium) scores next; O valence 1 (radical) scores
    # highest.  This prevents the algorithm from assigning the double bond to
    # the hydroxyl oxygen in carboxyl groups.
    _O_PREF = {2: 0, 3: 1, 1: 2}  # lower → preferred

    def _group_key(combo):
        groups: dict[int, list[int]] = {6: [], 7: [], 8: [], 15: [], 16: []}
        for v, an in zip(combo, atom_nrs):
            if an in groups:
                groups[an].append(v)
        return (
            tuple(_O_PREF.get(v, v) for v in groups[8]),  # O first (scored)
            tuple(groups[7]),  # N
            tuple(groups[6]),  # C
            tuple(groups[15]),  # P
            tuple(groups[16]),  # S
        )

    combos = list(itertools.product(*possible_valences))
    sorted_combos = sorted(combos, key=_group_key)

    best = dict(bond_orders)

    for valences in sorted_combos:
        UA, DU_from_AC = _get_UA(valences, conn_valence)

        if len(UA) == 0:
            if _BO_is_OK(
                bond_orders,
                adj_lookup,
                charge,
                DU_from_AC,
                atom_nrs,
                valences,
                atom_ids,
                allow_charged_fragments=allow_charged_fragments,
                allow_carbenes=allow_carbenes,
            ):
                return bond_orders
        else:
            UA_pairs_list = _get_UA_pairs(UA, adj_lookup, DU_from_AC)
            for UA_pairs in UA_pairs_list:
                candidate = _get_BO(
                    dict(best),
                    UA,
                    DU_from_AC,
                    valences,
                    UA_pairs,
                    adj_lookup,
                    atom_ids,
                )
                if _BO_is_OK(
                    candidate,
                    adj_lookup,
                    charge,
                    DU_from_AC,
                    atom_nrs,
                    valences,
                    atom_ids,
                    allow_charged_fragments=allow_charged_fragments,
                    allow_carbenes=allow_carbenes,
                ):
                    return candidate

                # Track best-so-far when the exact solution isn't found yet
                if (
                    sum(candidate.values()) >= sum(best.values())
                    and _valences_not_too_large(
                        candidate, valences, atom_ids, adj_lookup
                    )
                    and _charge_is_OK(
                        candidate,
                        adj_lookup,
                        charge,
                        DU_from_AC,
                        atom_nrs,
                        valences,
                        atom_ids,
                        allow_charged_fragments=allow_charged_fragments,
                        allow_carbenes=allow_carbenes,
                    )
                ):
                    best = candidate

    return best


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _get_UA(
    target_valences: Sequence[int], valence_list: list[int]
) -> tuple[list[int], list[int]]:
    """Find unsaturated atoms and their degree of unsaturation."""
    UA: list[int] = []
    DU: list[int] = []
    for i, (target, valence) in enumerate(zip(target_valences, valence_list)):
        if target - valence > 0:
            UA.append(i)
            DU.append(target - valence)
    return UA, DU


def _get_BO(
    bo: dict[frozenset[int], int],
    UA: Sequence[int],
    DU: Sequence[int],
    valences: Sequence[int],
    UA_pairs: tuple[tuple[int, int], ...],
    adj_lookup: list[list[bool]],
    atom_ids: list[int],
) -> dict[frozenset[int], int]:
    """Increment bond orders along UA pairs until all atoms are saturated."""
    bo = dict(bo)
    prev_DU: list[int] = []
    num_atoms = len(atom_ids)

    while prev_DU != list(DU):
        for i, j in UA_pairs:
            bond = frozenset({atom_ids[i], atom_ids[j]})
            bo[bond] = bo.get(bond, 0) + 1

        bo_total = [
            sum(
                bo.get(frozenset({atom_ids[i], atom_ids[j]}), 0)
                for j in range(num_atoms)
                if adj_lookup[i][j]
            )
            for i in range(num_atoms)
        ]
        prev_DU = list(DU)
        UA_new, DU_new = _get_UA(valences, bo_total)
        UA_pairs = _get_UA_pairs(UA_new, adj_lookup, DU_new)[0]

    return bo


def _valences_not_too_large(
    bo: dict[frozenset[int], int],
    valences: Sequence[int],
    atom_ids: list[int],
    adj_lookup: list[list[bool]],
) -> bool:
    """Check that no atom exceeds its maximum valence."""
    num_atoms = len(atom_ids)
    for i in range(num_atoms):
        total = sum(
            bo.get(frozenset({atom_ids[i], atom_ids[j]}), 0)
            for j in range(num_atoms)
            if adj_lookup[i][j]
        )
        if total > valences[i]:
            return False
    return True


def _BO_is_OK(
    bo: dict[frozenset[int], int],
    adj_lookup: list[list[bool]],
    charge: int,
    DU: list[int],
    atom_nrs: list[int],
    valences: Sequence[int],
    atom_ids: list[int],
    allow_charged_fragments: bool = True,
    allow_carbenes: bool = False,
) -> bool:
    """Check if the current bond order assignment is valid."""
    if not _valences_not_too_large(bo, valences, atom_ids, adj_lookup):
        return False

    extra_bonds = sum(v - 1 for v in bo.values())
    check_sum = extra_bonds == sum(DU)

    check_charge = _charge_is_OK(
        bo,
        adj_lookup,
        charge,
        DU,
        atom_nrs,
        valences,
        atom_ids,
        allow_charged_fragments,
        allow_carbenes=allow_carbenes,
    )
    return check_charge and check_sum


def _charge_is_OK(
    bo: dict[frozenset[int], int],
    adj_lookup: list[list[bool]],
    charge: int,
    DU: list[int],
    atom_nrs: list[int],
    valences: Sequence[int],
    atom_ids: list[int],
    allow_charged_fragments: bool = True,
    allow_carbenes: bool = False,
) -> bool:
    """Check whether computed atomic charges sum to the target charge."""
    q_tot = 0
    num_atoms = len(atom_ids)

    if allow_charged_fragments:
        for i, atomic_num in enumerate(atom_nrs):
            bo_total_i = sum(
                bo.get(frozenset({atom_ids[i], atom_ids[j]}), 0)
                for j in range(num_atoms)
                if adj_lookup[i][j]
            )
            q = _get_atomic_charge(
                atomic_num, atomic_valence_electrons[atomic_num], bo_total_i
            )
            q_tot += q
            if atomic_num == 6:
                n_single = sum(
                    1
                    for j in range(num_atoms)
                    if adj_lookup[i][j]
                    and bo.get(frozenset({atom_ids[i], atom_ids[j]}), 0) == 1
                )
                if not allow_carbenes and n_single == 2 and bo_total_i == 2:
                    q_tot += 1
                    q = 2
                if n_single == 3 and q_tot + 1 < charge:
                    q_tot += 2
                    q = 1

    return charge == q_tot


def _get_UA_pairs(
    UA: Sequence[int],
    adj_lookup: list[list[bool]],
    DU: Sequence[int] | None = None,
) -> list[tuple[()]] | list[tuple[tuple[int, int], ...]]:
    """Find a maximal set of non-overlapping unsaturated-atom pairs.

    Uses degree-ordered greedy matching (O(E log E)) — inspired by
    :func:`networkx.algorithms.matching.maximal_matching`.  A maximal
    matching is sufficient because :func:`_get_BO` iterates until all
    atoms are saturated.

    Atoms with degree of unsaturation > 1 are duplicated as virtual
    nodes so each unsaturation unit can be paired independently.
    """
    bonds = _get_bonds(UA, adj_lookup)

    if len(bonds) == 0:
        return [()]

    virtual_to_real: dict[int, int] = {}
    if DU is not None and any(d > 1 for d in DU):
        next_virtual = 10000
        real_to_virtual: dict[int, int] = {}
        for i, du in zip(UA, DU):
            if du > 1:
                v = next_virtual
                real_to_virtual[i] = v
                virtual_to_real[v] = i
                next_virtual += 1
        for i, j in list(bonds):
            if i in real_to_virtual:
                bonds.append((real_to_virtual[i], j))
            elif j in real_to_virtual:
                bonds.append((i, real_to_virtual[j]))

    # Degree-ordered greedy matching: prefer vertices with fewer options
    pairs = _maximal_matching(bonds)
    if not pairs:
        return [()]

    # Resolve any virtual nodes back to real atom indices
    if virtual_to_real:
        resolved = tuple(
            (virtual_to_real.get(i, i), virtual_to_real.get(j, j)) for i, j in pairs
        )
        return [resolved]

    return [tuple(pairs)]


def _maximal_matching(
    edges: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """Find a maximal matching by degree-ordered greedy selection.

    A matching is *maximal* when no edge can be added without sharing
    a vertex with an existing matched edge.  Sorting edges by endpoint
    degree (ascending) produces larger matchings than arbitrary order
    on typical molecular graphs.

    Runs in O(|E| log |E|) time.  Based on the approach used by
    :func:`networkx.algorithms.matching.maximal_matching`.
    """
    if not edges:
        return []

    # Vertex degree in the subgraph induced by `edges`
    degree: dict[int, int] = {}
    for u, v in edges:
        degree[u] = degree.get(u, 0) + 1
        degree[v] = degree.get(v, 0) + 1

    # Sort by combined degree — low-degree vertices get paired first
    edges_sorted = sorted(edges, key=lambda e: degree[e[0]] + degree[e[1]])

    matched: set[int] = set()
    result: list[tuple[int, int]] = []
    for u, v in edges_sorted:
        if u not in matched and v not in matched and u != v:
            matched.add(u)
            matched.add(v)
            result.append((u, v))

    return result


def _get_bonds(
    UA: Sequence[int],
    adj_lookup: list[list[bool]],
) -> list[tuple[int, int]]:
    """Return all edges between unsaturated atoms (as positional pairs)."""
    bonds: list[tuple[int, int]] = []
    for k, i in enumerate(UA):
        for j in UA[k + 1 :]:
            if adj_lookup[i][j]:
                bonds.append(tuple(sorted((i, j))))
    return bonds


def _get_atomic_charge(
    atomic_num: int, atomic_valence_electrons: int, bo_total: int
) -> int:
    """Compute formal charge from atomic number and bond order sum."""
    if atomic_num == 1:
        charge = 1 - bo_total
    elif atomic_num == 5:
        charge = 3 - bo_total
    elif atomic_num == 6 and bo_total == 2:
        charge = 0
    elif atomic_num == 13:
        charge = 3 - bo_total
    elif atomic_num == 15 and bo_total == 5:
        charge = 0
    elif atomic_num == 16 and bo_total == 6:
        charge = 0
    elif atomic_num == 16 and bo_total == 4:
        charge = 0
    elif atomic_num == 16 and bo_total == 5:
        charge = 1
    else:
        charge = atomic_valence_electrons - 8 + bo_total

    return charge


# MIT License

# Copyright (c) 2018 Jensen Group

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
