from __future__ import annotations

from collections import deque
from itertools import combinations
from typing import TYPE_CHECKING, Protocol

import numpy as np

from stereomolgraph import (
    COVALENT_RADII,
    PERIODIC_TABLE,
    Element
)

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Sequence, Optional
    from os import PathLike
    from typing import Any, Literal


def are_planar(*points: np.ndarray, threshold: float = 0.5) -> bool:
    """Checks if all atoms are in one plane

    Checks if the all atoms are planar within a given threshold.
    The threshold is the maximal distance of an atom from the plane of
    three other atoms.

    :param points: coordinates of the atoms
    :type points: np.ndarray
    :param threshold: maximal distance of atom from plane [Angstrom]
    :type threshold: float
    :return: True if all atoms are planar
    :rtype: bool
    """
    if threshold < 0:
        raise ValueError("threshold has to be bigger than 0")
    if len(points) < 4:
        return True

    for p1, p2, p3, p4 in combinations(points, 4):
        d = deque([p1, p2, p3, p4])
        for _ in range(4):
            d.rotate()
            vec1 = p1 - p2
            vec2 = p3 - p2
            vec3 = p4 - p2
            normal = np.cross(vec1, vec2)
            norm_normal = normal / np.linalg.norm(normal)
            result = abs(np.dot(norm_normal, vec3))
            if result > threshold:
                return False
    return True


def handedness(
    p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray
) -> Literal[1, -1]:
    """
    Calculates the orientation of the atom 4 from the plane defined
    by the first three atoms from their coordinates.

    :param p1: coordinates of atom 1
    :type p1: np.ndarray
    :param p2: coordinates of atom 2
    :type p2: np.ndarray
    :param p3: coordinates of atom 3
    :type p3: np.ndarray
    :param p4: coordinates of atom 4
    :type p4: np.ndarray
    :raises ValueError: if atoms are planar
    :return: Tetrahedral stereo
    :rtype: Tetrahedral
    """
    vec1 = p1 - p2
    vec2 = p3 - p2
    vec3 = p4 - p2
    normal = np.cross(vec1, vec2)
    norm_normal = normal / np.linalg.norm(normal)
    result: Literal[-1, 0, 1] = int(np.sign(np.dot(norm_normal, vec3)))
    if result == 0:
        raise ValueError("atoms are planar")
    return result

def pairwise_distances(
    coords: np.ndarray[tuple[int, Literal[3]], np.dtype[np.float64]],
) -> np.ndarray[tuple[int, Literal[3]], np.dtype[np.float64]]:
    diff = coords[:, None] - coords[None, :]
    diff **= 2
    summed = np.sum(diff, axis=-1)
    np.sqrt(summed, out=summed)
    return summed


class BaseGeometry(Protocol):
    atom_types: tuple[Element, ...]
    coords: np.ndarray[tuple[int, Literal[3]], Any]


class Geometry:
    """
    Represents a molecular geometry, i.e. the coordinates and atom types.

    :param atom_types: tuple of Element objects
    :param coords: nAtomsx3 numpy array with cartesian coordinates
    """

    atom_types: tuple[Element, ...]
    coords: np.ndarray[tuple[int, Literal[3]], np.dtype[np.float64]]

    @property
    def n_atoms(self) -> int:
        return len(self.atom_types)

    def __len__(self) -> int:
        return self.n_atoms

    def __init__(
        self,
        atom_types: Iterable[int | str | Element] = tuple(),
        coords: np.ndarray = np.empty((0, 3), dtype=np.float64),
    ):
        coords = np.array(coords, dtype=np.float64)
        atom_types = [PERIODIC_TABLE[atom] for atom in atom_types]

        assert (
            len(coords.shape) == 2
            and coords.shape[1] == 3
            and len(atom_types) == coords.shape[0]
        )

        self.coords = coords
        self.atom_types = tuple(atom_types)

    @classmethod
    def from_xyz_file(
        cls, path: PathLike
    ) -> BaseGeometry | tuple[BaseGeometry, str]:
        """Create a Geometry from an XYZ file."""

        dt = np.dtype(
            [
                ("atom", "U5"),  # Unicode string up to 5 characters
                ("x", "f8"),  # 64-bit float
                ("y", "f8"),
                ("z", "f8"),
            ]
        )

        data = np.loadtxt(path, skiprows=2, dtype=dt, comments=None)

        atom_types = tuple([PERIODIC_TABLE[atom] for atom in data["atom"]])
        coords = np.column_stack((data["x"], data["y"], data["z"]))

        return cls(atom_types=atom_types, coords=coords)

    def xyz_str(self, comment: Optional[str] = None) -> str:
        """
        returns the xyz representation of this geometry as a string

        :param comment: comment for 2nd line of xyz file
        :return: xyz representation of this geometry
        """
        xyz = str(self.n_atoms) + "\n"
        if comment is not None:
            xyz += comment + "\n"
        else:
            xyz += "\n"

        for type, coords in zip(self.atom_types, self.coords):
            xyz += (
                f"{type.symbol:s} {coords[0]:.8f} {coords[1]:.8f} "
                f"{coords[2]:.8f}\n"
            )

        return xyz


def default_connectivity_cutoff(atom_types: tuple[Element, Element]) -> float:
    return sum(COVALENT_RADII[a] for a in atom_types) * 1.2

CONNECTIVITY_CUTOFF_FUNC: Callable[[tuple[Element, Element]], float] = (
    default_connectivity_cutoff
)


class _DefaultFuncDict(dict):
    default_func: Callable
    """
    A dictionary that calls a default function with keys as arguments,
    when a key is missing.
    """

    def __init__(
        self,
        *args,
        default_func: Callable[
            [tuple[Element, Element]], float
        ] = CONNECTIVITY_CUTOFF_FUNC,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.default_func = default_func

    def __missing__(self, key: tuple[Element, Element]) -> float:
        if len(key) != 2:
            raise KeyError(f"Key {key} must be a tuple of two elements.")
        
        ret = self.default_func(key)
        
        self[key] = ret

        return ret

    def array(self, atom_types: Sequence[Element]) -> np.ndarray:
        n_atoms = len(atom_types)
        array = np.zeros((n_atoms, n_atoms))
        for (atom1, atom_type1), (atom2, atom_type2) in \
                combinations(enumerate(atom_types), 2):
            value = self[atom_type1, atom_type2]
            array[atom1][atom2] = value
            array[atom2][atom1] = value
        return array


class BondsFromDistance:

    def __init__(self, connectivity_cutoff: Callable[[tuple[Element, Element]], float] = CONNECTIVITY_CUTOFF_FUNC):
        self.connectivity_cutoff = _DefaultFuncDict(default_func=connectivity_cutoff)

    def __call__(self, distance:float, atom_types: tuple[Element, Element]):
        if distance < 0:
            raise ValueError('distance can not be negative')
        else:
            return 1 if distance < self.connectivity_cutoff[atom_types] else 0

    def array(self, coords, atom_types: Sequence[Element]) -> np.ndarray:
        return np.where(pairwise_distances(coords) < self.connectivity_cutoff.array(atom_types), 1, 0)

