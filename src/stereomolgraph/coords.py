from __future__ import annotations

from collections import deque
from itertools import combinations
from typing import TYPE_CHECKING

import numpy as np

from stereomolgraph import (
    COVALENT_RADII,
    PERIODIC_TABLE,
    Element,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Mapping, Sequence
    from os import PathLike
    from typing import Literal, TypeVar

    NP_FLOAT = TypeVar(
        "NP_FLOAT", bound=np.dtype[np.floating], contravariant=True
    )
    N = TypeVar("N", bound=int)
    ONE = Literal[1]
    #TWO = Literal[2]
    THREE = Literal[3]
    FOUR = Literal[4]

def are_planar(points: np.ndarray[tuple[int, Literal[3]], np.dtype[np.float64]],
               threshold: float = 0.5
               ) -> np.bool_:
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
    if threshold <= 0:
        raise ValueError("threshold has to be bigger than 0")
    if len(points) < 4:
        return np.bool_(True)

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
                return np.bool_(False)
    return np.bool_(True)

def are_planar_volume(
    coords: np.ndarray[tuple[FOUR, THREE], NP_FLOAT], threshold: float = 3.0
) -> np.ndarray[tuple[ONE], np.dtype[np.bool_]]:
    if coords.shape[-2] > 4:
        raise NotImplementedError(
            "are_planar_volume is not implemented for more than 4 points"
        )
        #for comb in combinations(range(coords.shape[-2]), 4):
        #    if not are_planar_volume(coords[..., list(comb), :], threshold=threshold):
        #        return np.array([False], dtype=np.bool_)
        #    else:
        #        return np.array([True], dtype=np.bool_)

    assert coords.shape[-1] == 3
    assert coords.shape[-2] == 4

    v1 = coords[..., 0, :] - coords[..., 1, :]
    v2 = coords[..., 2, :] - coords[..., 3, :]
    v3 = coords[..., 3, :] - coords[..., 1, :]

    normal = np.cross(v2, v3, axis=-1)

    volumes = np.abs(np.sum(v1 * normal, axis=-1))

    return volumes < threshold


def handedness(
    coords: np.ndarray[tuple[FOUR, THREE], NP_FLOAT],
) -> np.ndarray[tuple[ONE], np.dtype[np.int8]]:
    """
    Calculates the orientation of the atom 4 from the plane defined
    by the first three atoms from their coordinates.


    """
    vec1 = coords[..., 0, :] - coords[..., 1, :]
    vec2 = coords[..., 2, :] - coords[..., 3, :]
    vec3 = coords[..., 3, :] - coords[..., 1, :]

    normal = np.cross(vec1, vec2, axis=-1)

    # Normalize normal vectors (in-place operations)
    norm = np.linalg.norm(normal, axis=-1, keepdims=True)
    np.divide(normal, norm, out=normal)  # in-place division

    # Compute dot product (in-place multiply and sum)
    np.multiply(normal, vec3, out=normal)
    dot_product = np.sum(normal, axis=-1)
    result = np.sign(dot_product).astype(np.int8)
    return result


def angle_from_coords(
    coords: np.ndarray[tuple[THREE, THREE], NP_FLOAT],
    out: None | np.ndarray[tuple[ONE], NP_FLOAT] = None,
) -> np.ndarray[tuple[int, ...], NP_FLOAT]:
    """ """
    assert np.issubdtype(coords.dtype, np.floating)
    assert coords.shape[-1] == 3
    assert coords.shape[-2] == 3

    BA = coords[..., 0, :] - coords[..., 1, :]
    BC = coords[..., 2, :] - coords[..., 1, :]

    # Compute dot product and norms
    dot_product = np.einsum("...i,...i->...", BA, BC)
    assert dot_product is not None
    norm_BA = np.linalg.norm(BA, axis=-1)
    norm_BC = np.linalg.norm(BC, axis=-1)

    if out is not None:
        assert out.shape == dot_product.shape, "Output array must match shape"
        assert out.dtype == coords.dtype, "Output array must match dtype"
    else:
        output_shape = coords.shape[:-2]
        if output_shape == ():
            output_shape = (1,)

        out = np.empty(output_shape, dtype=coords.dtype)  # type: ignore
        assert out is not None, "Output array must be provided or created"
        assert out.shape == output_shape, (
            "Output array must be created with correct shape"
        )

    # Compute cosine of angle with safe division
    denominator = np.multiply(norm_BA, norm_BC, out=out)  # Reuse out as temp
    np.divide(dot_product, denominator, where=(denominator != 0), out=out)

    # Clip and compute angle
    np.clip(out, -1.0, 1.0, out=out)
    np.arccos(out, out=out)
    np.rad2deg(out, out=out)
    assert out is not None
    return out  # .squeeze()[()] if out.ndim == 0 else out


def pairwise_distances(
    coords: np.ndarray[tuple[N, THREE], NP_FLOAT],
) -> np.ndarray[tuple[N, N], NP_FLOAT]:
    if coords.shape[-1] != 3:
        raise ValueError("Last dimension must be size 3 for 3D coordinates")

    # Compute differences using broadcasting
    diff = (
        coords[..., :, None, :] - coords[..., None, :, :]
    )  # shape (..., N, N, 3)

    # Square differences and sum along last dimension
    np.square(diff, out=diff)  # in-place squaring
    summed = np.sum(diff, axis=-1)  # shape (..., N, N)

    # Compute distances in-place
    distances = np.sqrt(summed, out=summed)
    return distances


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
        coords: np.ndarray[tuple[int, int], NP_FLOAT] = np.empty(
            (0, 3), dtype=np.float64
        ),
    ):
        self.coords = np.array(coords, dtype=np.float64)
        self.atom_types = tuple([PERIODIC_TABLE[atom] for atom in atom_types])

        assert (
            len(self.coords.shape) == 2
            and self.coords.shape[1] == 3
            and len(self.atom_types) == self.coords.shape[0]
        )

    @classmethod
    def from_xyz_file(cls, path: PathLike[str] | str) -> Geometry:
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

        atom_types = [PERIODIC_TABLE[atom] for atom in data["atom"]]
        coords = np.column_stack((data["x"], data["y"], data["z"]))

        return cls(atom_types=atom_types, coords=coords)

    def xyz_str(self, comment: None | str = None) -> str:
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


class _DefaultFuncDict(dict[tuple[Element, Element], float]):
    """
    A dictionary that calls a default function with keys as arguments,
    when a key is missing.
    """

    default_func: Callable[[tuple[Element, Element]], float]

    def __init__(
        self,
        *,
        default_func: Callable[
            [tuple[Element, Element]], float
        ] = CONNECTIVITY_CUTOFF_FUNC,
        **kwargs: Mapping[tuple[Element, Element], float],
    ):
        super().__init__(**kwargs)
        self.default_func = default_func

    def __missing__(self, key: tuple[Element, Element]) -> float:
        if len(key) != 2:
            raise KeyError(f"Key {key} must be a tuple of two elements.")
        elif (ret := self.get((key[1], key[0]), None)) is not None:
            pass
        else:
            ret = self.default_func(key)

        self[key] = ret

        return ret

    def array(self, atom_types: Sequence[Element]) -> np.ndarray:
        n_atoms = len(atom_types)
        array = np.zeros((n_atoms, n_atoms))
        for (atom1, atom_type1), (atom2, atom_type2) in combinations(
            enumerate(atom_types), 2
        ):
            value = self[atom_type1, atom_type2]
            array[atom1][atom2] = value
            array[atom2][atom1] = value
        return array


class BondsFromDistance:
    def __init__(
        self,
        connectivity_cutoff: Callable[
            [tuple[Element, Element]], float
        ] = CONNECTIVITY_CUTOFF_FUNC,
    ):
        self.connectivity_cutoff = _DefaultFuncDict(
            default_func=connectivity_cutoff
        )

    def __call__(
        self, distance: float, atom_types: tuple[Element, Element]
    ) -> Literal[0, 1]:
        if distance < 0:
            raise ValueError("distance can not be negative")
        else:
            return 1 if distance < self.connectivity_cutoff[atom_types] else 0

    def array(
        self,
        coords: np.ndarray[tuple[int, Literal[3]], np.dtype[np.float64]],
        atom_types: Sequence[Element],
    ) -> np.ndarray:
        return np.where(
            pairwise_distances(coords)
            < self.connectivity_cutoff.array(atom_types),
            1,
            0,
        )
