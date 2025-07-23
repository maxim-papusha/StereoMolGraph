from __future__ import annotations
import sys
import itertools
from abc import ABC, abstractmethod
from collections import Counter
from typing import TYPE_CHECKING, Any, Generic, Literal, TypeVar

import numpy as np

from stereomolgraph.coords import are_planar, handedness, angle_from_coords

if TYPE_CHECKING:
    from collections.abc import Generator, Set

    from typing_extensions import Self

    from stereomolgraph.graphs.mg import AtomId, Bond

if sys.version_info >= (3, 13):
    from typing import TypeVar
else:
    from typing_extensions import TypeVar


A = TypeVar("A", bound=tuple[int, ...],
            covariant=True,
            default=tuple[int, ...])
P = TypeVar("P", bound=None | Literal[1, 0, -1],
            covariant=True,
            default=None | Literal[1, 0, -1])


class Stereo(ABC, Generic[A, P]):
    """
    :class:`~typing.Generic` Class to represent the orientation of a group of
    atoms in space.
    This is used to represent local stereochemistry and simultanously the
    hybridization of atoms.
    """

    atoms: A
    """Atoms are a order dependent tuple of integers."""

    parity: P
    """parity is a number that defines the orientation of the atoms. If None,
    the relative orientation of the atoms is not defined.
    If 0 the orientation is defined and part of a achiral stereochemistry.
    If 1 or -1 the orientation is defined and part of a chiral stereochemistry.
    """

    PERMUTATION_GROUP: frozenset[A]
    """Defines all allowed permutations defined by the symmetry group under
    which the stereochemistry is invariant."""

    @abstractmethod
    def __init__(self, atoms: A, parity: P = None): ...

    @abstractmethod
    def __eq__(self, other: Any) -> bool: ...

    @abstractmethod
    def __hash__(self) -> int: ...

    @abstractmethod
    def invert(self) -> Self:
        """Inverts the stereo. If the stereo is achiral, it returns itself."""

    @abstractmethod
    def get_isomers(self) -> Set[Self]:
        """Returns all stereoisomers of the stereochemistry. Not just the
        inverted ones, but all possible stereoisomers."""


class AtomStereo(Stereo[A, P], ABC, Generic[A, P]):
    @property
    def central_atom(self) -> AtomId:
        return self.atoms[0]


class BondStereo(Stereo[A, P], ABC, Generic[A, P]):
    @property
    def bond(self) -> Bond:
        bond = frozenset(self.atoms[2:4])
        assert len(bond) == 2
        return bond


class _StereoMixin(Stereo[A, P], ABC, Generic[A, P]):
    __slots__ = ("atoms", "parity")

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.atoms}, {self.parity})"

    def __init__(self, atoms: A, parity: P = None):
        self.atoms = atoms
        self.parity = parity

    def _perm_atoms(self) -> Generator[tuple[int, ...], None, None]:
        if self.parity is None:
            return (
                tuple([self.atoms[i] for i in perm])
                for perm in itertools.permutations(range(len(self.atoms)))
            )
        else:
            return (
                tuple([self.atoms[i] for i in perm])
                for perm in self.PERMUTATION_GROUP
            )


class _ChiralStereoMixin(
    _StereoMixin[A, None | Literal[1, -1]], ABC, Generic[A]
):
    __slots__ = ("inversion",)
    inversion: A

    def invert(self) -> Self:
        if self.parity is None:
            return self
        return self.__class__(self.atoms, self.parity * -1)

    def _inverted_atoms(self) -> A:
        atoms = tuple([self.atoms[i] for i in self.inversion])
        assert len(atoms) == len(self.atoms) == len(self.inversion)
        return atoms  # type: ignore[return-value]

    def __eq__(self, other: Any) -> bool:
        if other.parity == 0:
            return False

        s_atoms, o_atoms = self.atoms, other.atoms
        set_s_atoms = set(s_atoms)
        set_o_atoms = set(o_atoms)

        if len(s_atoms) != len(o_atoms) or not set_s_atoms.issuperset(
            set_o_atoms
        ):
            return False

        if self.parity is None or other.parity is None:
            if set_s_atoms == set_o_atoms:
                return True
            return False

        elif self.parity == other.parity:
            if o_atoms == s_atoms or any(
                o_atoms == p for p in self._perm_atoms()
            ):
                return True
            return False

        elif self.parity * -1 == other.parity:
            if any(other._inverted_atoms() == p for p in self._perm_atoms()):
                return True
            return False

        raise RuntimeError("This should not happen")

    def __hash__(self) -> int:
        if self.parity is None:
            return hash(frozenset(Counter(self.atoms).items()))
        perm = frozenset(
            {
                tuple([self.atoms[i] for i in perm])
                for perm in self.PERMUTATION_GROUP
            }
        )

        inverted_perm = frozenset(
            {
                tuple([self._inverted_atoms()[i] for i in perm])
                for perm in self.PERMUTATION_GROUP
            }
        )

        if self.parity == 1:
            return hash((perm, inverted_perm))
        elif self.parity == -1:
            return hash((inverted_perm, perm))
        else:
            raise RuntimeError("Something is wrong with parity")


class _AchiralStereoMixin(_StereoMixin[A, None | Literal[0]], ABC, Generic[A]):
    __slots__: tuple[str, ...] = ()

    def invert(self) -> Self:
        return self

    def __eq__(self, other: Any) -> bool:
        if other.parity in (1, -1):
            return False

        s_atoms, o_atoms = self.atoms, other.atoms
        set_s_atoms = set(s_atoms)
        set_o_atoms = set(o_atoms)

        if len(s_atoms) != len(o_atoms) or not set_s_atoms.issuperset(
            set_o_atoms
        ):
            return False

        if self.parity is None or other.parity is None:
            if set_s_atoms == set_o_atoms:
                return True
            return False

        if self.parity == other.parity:
            if o_atoms == s_atoms or o_atoms in self._perm_atoms():
                return True
            return False

        raise RuntimeError("This should not happen!")

    def __hash__(self) -> int:
        if self.parity is None:
            return hash(frozenset(Counter(self.atoms).items()))
        elif self.parity == 0:
            perm = frozenset(
                {
                    tuple([self.atoms[i] for i in perm])
                    for perm in self.PERMUTATION_GROUP
                }
            )
            return hash(perm)
        raise RuntimeError("Something is wrong with parity")


class Tetrahedral(
    _ChiralStereoMixin[tuple[int, int, int, int, int]],
    AtomStereo[tuple[int, int, int, int, int], None | Literal[1, -1]],
):
    r"""Represents all possible configurations of atoms for a Tetrahedral
    Stereochemistry::

       parity = 1      parity = -1
           4                4
           |                |
           0                0
        /  ¦  \          /  ¦  \
       2   1   3        3   1   2

    Atoms of the tetrahedral stereochemistry are ordered in a way that when the
    first atom is rotated to the back, the other atoms in order are rotated in
    the direction defined by the stereo.

    :ivar atoms: Atoms of the stereochemistry
    :ivar parity: Stereochemistry
    :ivar PERMUTATION_GROUP: Permutations allowed by the stereochemistry
    """

    __slots__ = ()
    inversion = (0, 2, 1, 3, 4)
    PERMUTATION_GROUP = frozenset(
        {
            (0, 1, 2, 3, 4),
            (0, 3, 1, 2, 4),
            (0, 2, 3, 1, 4),
            (0, 1, 4, 2, 3),
            (0, 2, 1, 4, 3),
            (0, 4, 2, 1, 3),
            (0, 1, 3, 4, 2),
            (0, 4, 1, 3, 2),
            (0, 3, 4, 1, 2),
            (0, 2, 4, 3, 1),
            (0, 3, 2, 4, 1),
            (0, 4, 3, 2, 1),
        }
    )

    def __init__(
        self,
        atoms: tuple[int, int, int, int, int],
        parity: None | Literal[1, -1] = None,
    ):
        super().__init__(atoms=atoms, parity=parity)

    def get_isomers(self) -> set[Self]:
        return {
            self.__class__(atoms=self.atoms, parity=1),
            self.__class__(atoms=self.atoms, parity=-1),
        }

    @classmethod
    def from_coords(
        cls,
        atoms: tuple[int, int, int, int, int],
        coords: np.ndarray[
            tuple[Literal[5], Literal[3]], np.dtype[np.float64]
        ],
    ) -> Self:
        """
        Creates the representation of a Tetrahedral Stereochemistry
        from the coordinates of the atoms.

        :param atoms: Atoms of the stereochemistry
        :param coords: nAtomsx3 numpy array with cartesian coordinates
        """
        orientation = handedness(coords.take((1, 2, 3, 4), axis=0))
        int_orientation = int(orientation)
        assert int_orientation in (1, -1), (
            f"Orientation {orientation} is not valid for Tetrahedral "
            "stereochemistry.")
        return cls(atoms, int_orientation)


class SquarePlanar(
    _AchiralStereoMixin[tuple[int, int, int, int, int]],
    AtomStereo[tuple[int, int, int, int, int], None | Literal[0]],
):
    r""" Represents all possible configurations of atoms for a
    SquarePlanar Stereochemistry::

        1     4
         \   /
           0
         /   \
        2     3

    Atoms of the Square Planar stereochemistry are ordered in a way that


    :ivar atoms: Atoms of the stereochemistry
    :ivar parity: Stereochemistry
    """

    __slots__ = ()
    PERMUTATION_GROUP = frozenset(
        {
            (0, 1, 2, 3, 4),
            (0, 2, 3, 4, 1),
            (0, 3, 4, 1, 2),
            (0, 4, 1, 2, 3),
            (0, 4, 3, 2, 1),
            (0, 3, 2, 1, 4),
            (0, 2, 1, 4, 3),
            (0, 1, 4, 3, 2),
        }
    )

    @property
    def central_atom(self) -> AtomId:
        return self.atoms[0]

    def get_isomers(self) -> set[SquarePlanar]:
        return {
            SquarePlanar(atoms=atoms, parity=0)
            for perm in itertools.permutations(self.atoms[1:])
            if len(atoms := (self.atoms[0], *perm)) == 5
        }


class TrigonalBipyramidal(
    _ChiralStereoMixin[tuple[int, int, int, int, int, int]],
    AtomStereo[tuple[int, int, int, int, int, int], None | Literal[1, -1]],
):
    r"""Represents all possible configurations of atoms for a
    TrigonalBipyramidal Stereochemistry::

       parity = 1             parity = -1
        3   1                     1   3
         ◁  ¦                    ¦  ▷
            0  — 5           5 —  0
         ◀  ¦                    ¦  ▶
        4   2                     2   4

    Atoms of the trigonal bipyramidal stereochemistry are ordered in a way that
    when the first two atoms are the top and bottom of the bipyramid. The last
    three equatorial atoms are ordered in a way that when the first atom is
    rotated to the back, the other atoms in order are rotated in the direction
    defined by the stereo.

    :ivar atoms: Atoms of the stereochemistry
    :ivar parity: Stereochemistry
    """

    __slots__ = ()
    inversion = (0, 1, 2, 3, 5, 4)
    PERMUTATION_GROUP = frozenset(
        {
            (0, 1, 2, 3, 4, 5),
            (0, 1, 2, 5, 3, 4),
            (0, 1, 2, 4, 5, 3),
            (0, 2, 1, 3, 5, 4),
            (0, 2, 1, 5, 4, 3),
            (0, 2, 1, 4, 3, 5),
        }
    )

    def get_isomers(self) -> set[Self]:
        return {
            self.__class__(atoms=atoms, parity=p)
            for perm in itertools.permutations(self.atoms[1:])
            for p in (1, -1)
            if len(atoms := (self.atoms[0], *perm)) == 6
        }

    @classmethod
    def from_coords(
        cls: type[TrigonalBipyramidal],
        atoms: tuple[int, int, int, int, int, int],
        coords: np.ndarray[
            tuple[Literal[6], Literal[3]], np.dtype[np.float64]
        ],
    ) -> TrigonalBipyramidal:
        """
        calculates the distance of the atom 5 from the plane defined by the
        first three atoms in Angstrom. The sign of the distance is determined
        by the side of the plane that atom 5 is on.
        """

        # coords_dict
        _index_central_atom = 0
        indices = (1, 2, 3, 4, 5)

        if np.any(are_planar(coords[[1,2,3,4,5]])):
            raise ValueError("Four atoms are planar!")
        
        lst = np.array([[i, 0, j] for i, j
                        in itertools.combinations(indices, 2)], dtype=np.int8)

        # The atoms with the largest angle are the axial atoms
        angles = angle_from_coords(coords[lst])
        i, j = lst[angles.argmax()][[0, 2]] # axial atoms
        
        equatorial = [a for a in indices if a not in (i, j)]
        i_rotation = handedness(coords.take( [*equatorial, i], axis=0))
        j_rotation = -1 * handedness(coords.take([*equatorial, j], axis=0))

        assert int(i_rotation) == int(j_rotation)

        atoms_in_new_order = (i, j, *equatorial)
        orientation = int(i_rotation)
        tb_atoms = (atoms[0], *atoms_in_new_order)
        assert len(tb_atoms) == 6
        assert orientation in (1, -1)
        return TrigonalBipyramidal(tb_atoms, orientation)
        
class Octahedral(
    _ChiralStereoMixin[tuple[int, int, int, int, int, int, int]],
    AtomStereo[
        tuple[int, int, int, int, int, int, int], None | Literal[1, -1]
    ],
):
    """Represents all possible configurations of atoms for a Octahedral
    Stereochemistry::

        parity = 1             parity = -1
         3  1   6                3  2  6
          ◁ ¦ /                  ◁ ¦ /
            0                       0
          / ¦ ▶                  / ¦  ▶
         4  2  5                4   1  5
    """

    __slots__ = ()
    inversion = (0, 2, 1, 3, 4, 5, 6)
    PERMUTATION_GROUP = frozenset(
        {
            (0, 1, 2, 3, 4, 5, 6),
            (0, 1, 2, 6, 3, 4, 5),
            (0, 1, 2, 5, 6, 3, 4),
            (0, 1, 2, 4, 5, 6, 3),
            (0, 2, 1, 4, 3, 6, 5),
            (0, 2, 1, 5, 4, 3, 6),
            (0, 2, 1, 6, 5, 4, 3),
            (0, 2, 1, 3, 6, 5, 4),
            (0, 3, 5, 2, 4, 1, 6),
            (0, 3, 5, 6, 2, 4, 1),
            (0, 3, 5, 1, 6, 2, 4),
            (0, 3, 5, 4, 1, 6, 2),
            (0, 5, 3, 1, 4, 2, 6),
            (0, 5, 3, 6, 1, 4, 2),
            (0, 5, 3, 2, 6, 1, 4),
            (0, 5, 3, 4, 2, 6, 1),
            (0, 4, 6, 3, 2, 5, 1),
            (0, 4, 6, 1, 3, 2, 5),
            (0, 4, 6, 5, 1, 3, 2),
            (0, 4, 6, 2, 5, 1, 3),
            (0, 6, 4, 3, 1, 5, 2),
            (0, 6, 4, 2, 3, 1, 5),
            (0, 6, 4, 5, 2, 3, 1),
            (0, 6, 4, 1, 5, 2, 3),
        }
    )

    def get_isomers(self) -> set[Octahedral]:
        return {
            Octahedral(atoms=atoms, parity=p)
            for perm in itertools.permutations(self.atoms[1:])
            for p in (1, -1)
            if len((atoms := (self.atoms[0], *perm))) == 7
        }


class PlanarBond(
    _AchiralStereoMixin[tuple[int, int, int, int, int, int]],
    BondStereo[tuple[int, int, int, int, int, int], None | Literal[0]],
):
    r""" Represents all possible configurations of atoms for a
    Planar Structure and should be used for aromatic and double bonds::

        0        4
         \      /
          2 == 3
         /      \
        1        5

    All atoms of the double bond are in one plane. Atoms 2 and 3 are the center
    Atoms 0 and 1 are bonded to 2 and atoms 4 and 5 are bonded to 3.
    The stereochemistry is defined by the relative orientation
    of the atoms 0, 1, 4 and 5.

    :ivar atoms: Atoms of the stereochemistry
    :ivar parity: Stereochemistry
    :ivar PERMUTATION_GROUP: Permutations allowed by the stereochemistry
    """

    __slots__ = ()
    PERMUTATION_GROUP = frozenset(
        {
            (0, 1, 2, 3, 4, 5),
            (1, 0, 2, 3, 5, 4),
            (4, 5, 3, 2, 0, 1),
            (5, 4, 3, 2, 1, 0),
        }
    )

    def get_isomers(self) -> set[PlanarBond]:
        reordered_atoms = tuple(self.atoms[i] for i in (0, 1, 2, 3, 5, 4))
        assert len(reordered_atoms) == 6
        return {
            PlanarBond(self.atoms, 0),
            PlanarBond(reordered_atoms, 0),
        }

    def __init__(
        self,
        atoms: tuple[int, int, int, int, int, int],
        parity: None | Literal[0] = None,
    ):
        super().__init__(atoms, parity)

    @classmethod
    def from_coords(
        cls,
        atoms: tuple[int, int, int, int, int, int],
        coords: np.ndarray[
            tuple[Literal[6], Literal[3]], np.dtype[np.float64]
        ],
    ) -> PlanarBond:
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
        return cls(new_atoms, 0)


class AtropBond(
    _ChiralStereoMixin[tuple[int, int, int, int, int, int]],
    BondStereo[tuple[int, int, int, int, int, int], None | Literal[1, -1]],
):
    r"""
    Represents all possible configurations of atoms for a
    Atropostereoisomer bond::
    
        parity = 1          parity = -1
        1       5           1        5
         \     /            ◀      /
          2 - 3               2 - 3
        ◀      \            /      \
        0        4         0         4


    """

    __slots__ = ()
    inversion = (1, 0, 2, 3, 4, 5)
    PERMUTATION_GROUP = frozenset(
        {
            (0, 1, 2, 3, 4, 5),
            (1, 0, 2, 3, 5, 4),
            (4, 5, 3, 2, 1, 0),
            (5, 4, 3, 2, 0, 1),
        }
    )

    def get_isomers(self) -> set[AtropBond]:
        other_atoms = tuple(self.atoms[i] for i in (0, 1, 2, 3, 5, 4))
        assert len(other_atoms) == 6
        return {
            AtropBond(self.atoms, 1),
            AtropBond(other_atoms, -1),
        }
