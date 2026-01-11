from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import numpy as np

from stereomolgraph.coords import (
    BondsFromDistance,
    angle_from_coords,
    are_planar,
    handedness,
)
from stereomolgraph.stereodescriptors import (
    AtomStereo,
    Octahedral,
    PlanarBond,
    SquarePlanar,
    Tetrahedral,
    TrigonalBipyramidal,
)

if TYPE_CHECKING:
    from typing import Literal, TypeVar

    from stereomolgraph import MolGraph, StereoMolGraph
    from stereomolgraph.coords import BondsFromDistance, Geometry

    N = TypeVar("N", bound=int, covariant=True)
    MG = TypeVar("MG", bound=MolGraph, covariant=True)
    SMG = TypeVar("SMG", bound=StereoMolGraph, covariant=True)
    NP_FLOAT = TypeVar(
        "NP_FLOAT", bound=np.dtype[np.floating], contravariant=True
    )

    THREE = Literal[3]


def connectivity_from_geometry(
    cls: type[MG],
    geo: Geometry,
    switching_function: BondsFromDistance = BondsFromDistance(),
) -> MG:
    """
    Creates a graph of a molecule from a Geometry and a switching Function.
    Uses the Default switching function if none are given.

    :param geo: Geometry
    :param switching_function: Function to determine if two atoms are
        connected
    :return: graph of molecule
    """

    connectivity_matrix = switching_function.array(geo.coords, geo.atom_types)
    return cls.from_atom_types_and_bond_order_matrix(
        geo.atom_types,
        connectivity_matrix,
    )


def stero_from_geometry(
    smg: SMG,
    geo: Geometry,
) -> SMG:
    for atom in range(geo.n_atoms):
        first_nbrs = smg.bonded_to(atom)
        atom_stereo_tup = (atom, *first_nbrs)
        atom_stereo = atom_stereo_from_coords(
            atom_stereo_tup, geo.coords.take(atom_stereo_tup, axis=0)
        )

        if atom_stereo is not None:
            smg.set_atom_stereo(atom_stereo)

        if len(first_nbrs) != 3:
            continue

        for nbr in first_nbrs:
            second_neighbors = smg.bonded_to(nbr).difference({atom})

            if len(second_neighbors) != 2:
                continue

            first_nbrs_reduced = first_nbrs.difference({nbr})
            stereo_atoms = (*first_nbrs_reduced, atom, nbr, *second_neighbors)
            assert len(stereo_atoms) == 6
            planar_bond = _planar_bond_from_coords(
                stereo_atoms, geo.coords.take(stereo_atoms, axis=0)
            )
            if planar_bond is not None:
                smg.set_bond_stereo(planar_bond)

    return smg


def atom_stereo_from_coords(
    atoms: tuple[int, ...],
    coords: np.ndarray[tuple[N, THREE], NP_FLOAT],
) -> None | AtomStereo:
    assert len(atoms) == coords.shape[0]
    assert coords.shape[1] == 3
    assert np.issubdtype(coords.dtype, np.floating)

    if len(atoms) == 5 and coords.shape[0] == 5:
        if are_planar(coords[[1, 2, 3, 4]]):
            return _square_planar_from_coords(atoms=atoms, coords=coords)  # type: ignore
        else:
            return _tetrahedral_from_coords(atoms=atoms, coords=coords)  # type: ignore

    elif len(atoms) == 6:
        return _trigonal_bipyramidal_from_coords(atoms=atoms, coords=coords)  # type: ignore

    elif len(atoms) == 7:
        return _octahedral_from_coords(atoms=atoms, coords=coords)  # type: ignore

    else:
        return None


def _tetrahedral_from_coords(
    atoms: tuple[int, int, int, int, int],
    coords: np.ndarray[tuple[Literal[5], THREE], NP_FLOAT],
) -> Tetrahedral:
    """
    Creates the representation of a Tetrahedral Stereochemistry
    from the coordinates of the atoms.

    :param atoms: Atoms of the stereochemistry
    :param coords: nAtomsx3 numpy array with cartesian coordinates
    """
    indeces = (1, 2, 3, 4)
    orientation = handedness(coords.take(indeces, axis=0))
    int_orientation = int(orientation)
    assert int_orientation in (1, -1), (
        f"Orientation {orientation} is not valid for Tetrahedral "
        "stereochemistry."
    )
    return Tetrahedral(atoms, int_orientation)


def _square_planar_from_coords(
    atoms: tuple[int, int, int, int, int],
    coords: np.ndarray[tuple[Literal[5], THREE], NP_FLOAT],
) -> None | SquarePlanar:
    if not are_planar(coords[[1, 2, 3, 4]]):
        return None

    distinct_orders = ((1, 2, 3, 4), (2, 1, 3, 4), (3, 2, 1, 4), (4, 2, 3, 1))

    max_order: tuple[int, int, int, int] = (1, 2, 3, 4)
    angle_sum = 0.0
    for order in distinct_orders:
        a = angle_from_coords(coords[[order[i] for i in (0, 1, 2)]])
        b = angle_from_coords(coords[[order[i] for i in (1, 2, 3)]])
        c = angle_from_coords(coords[[order[i] for i in (2, 3, 0)]])
        d = angle_from_coords(coords[[order[i] for i in (3, 0, 1)]])
        angle = a + b + c + d
        if angle > angle_sum:
            angle_sum, max_order = angle, order

    atoms = (
        atoms[0],
        atoms[max_order[0]],
        atoms[max_order[1]],
        atoms[max_order[2]],
        atoms[max_order[3]],
    )

    return SquarePlanar(atoms, parity=0)


def _trigonal_bipyramidal_from_coords(
    atoms: tuple[int, int, int, int, int, int],
    coords: np.ndarray[tuple[Literal[6], THREE], NP_FLOAT],
) -> None | TrigonalBipyramidal:
    """
    calculates the distance of the atom 5 from the plane defined by the
    first three atoms in Angstrom. The sign of the distance is determined
    by the side of the plane that atom 5 is on.
    """
    indices = (1, 2, 3, 4, 5)

    if np.any(are_planar(coords[[1, 2, 3, 4, 5]])):
        return None

    lst = np.array(
        [[i, 0, j] for i, j in itertools.combinations(indices, 2)],
        dtype=np.int8)

    # The atoms with the largest angle are the axial atoms
    angles = angle_from_coords(coords[lst])

    i, j = lst[angles.argmax()][[0, 2]]  # axial atoms
    i, j = int(i), int(j)

    equatorial = [a for a in indices if a not in (i, j)]
    i_rotation = -1 * handedness(coords.take([*equatorial, i], axis=0))
    j_rotation = handedness(coords.take([*equatorial, j], axis=0))

    assert int(i_rotation) == int(j_rotation)

    indeces_in_new_order = (i, j, *equatorial)

    orientation = int(i_rotation)
    tb_indeces = (0, *indeces_in_new_order)

    tb_atoms = tuple(atoms[i] for i in tb_indeces)

    assert len(tb_atoms) == 6
    assert orientation in (1, -1)
    return TrigonalBipyramidal(tb_atoms, orientation)


def _octahedral_from_coords(
    atoms: tuple[int, int, int, int, int, int, int],
    coords: np.ndarray[tuple[Literal[7], THREE], NP_FLOAT],
) -> None | Octahedral:
    indeces = (1, 2, 3, 4, 5, 6)
    planar_groups: list[set[int]] = []
    for p1, p2, p3, p4 in itertools.combinations(indeces, 4):
        points = [p1, p2, p3, p4]
        if are_planar(coords[points], threshold=0.3):
            planar_groups.append(set(points))

    if len(planar_groups) != 3:
        return None

    trans_atoms = planar_groups[0].intersection(planar_groups[1])

    cis_atoms0 = planar_groups[0].difference(trans_atoms)

    cis_atoms1 = planar_groups[1].difference(trans_atoms)
    assert len(trans_atoms) == 2
    assert len(cis_atoms0) == 2
    assert len(cis_atoms1) == 2
    assert cis_atoms0 | cis_atoms1 == planar_groups[2]

    a1, a2 = trans_atoms
    a3, a5 = cis_atoms0
    a4, a6 = cis_atoms1



    parity = int(handedness(coords[[a1, a3, a4, a2]]))
    parity2 = int(handedness(coords[[a1, a3, a6, a2]]))
    assert parity2 == -parity
    assert parity == 1 or parity == -1
    return Octahedral((atoms[0], atoms[a1], atoms[a2], atoms[a3], atoms[a4], atoms[a5], atoms[a6]), parity)


def _planar_bond_from_coords(
    atoms: tuple[int, int, int, int, int, int],
    coords: np.ndarray[tuple[Literal[6], Literal[3]], np.dtype[np.floating]],
) -> None | PlanarBond:
    if not are_planar(coords):
        return None

    a = (coords[0] - coords[1]) / np.linalg.norm(coords[0] - coords[1])
    b = (coords[4] - coords[5]) / np.linalg.norm(coords[4] - coords[5])
    result = int(np.sign(np.dot(a, b)))

    if result == -1:
        new_atoms = tuple(atoms[i] for i in (1, 0, 2, 3, 4, 5))
    elif result == 1:
        new_atoms = atoms
    elif result == 0:
        raise ValueError("atoms are tetrahedral")
    else:
        raise ValueError("something went wrong")
    assert len(new_atoms) == 6
    return PlanarBond(new_atoms, 0)
