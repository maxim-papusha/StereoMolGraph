from __future__ import annotations

import itertools
from abc import ABC, abstractmethod
from collections import Counter
from typing import TYPE_CHECKING, Protocol

import numpy as np

from stereomolgraph.cartesian import are_planar, handedness

if TYPE_CHECKING:
    import sys
    from collections.abc import Iterator
    from typing import Any, Literal, Optional

    from stereomolgraph.graph import AtomId

    # Self is included in typing from 3.11
    if sys.version_info >= (3, 11):
        from typing import Self
    else:
        from typing_extensions import Self


class Stereo(Protocol):
    """Base Class to represent the orientation of a group of atoms in space and
    their allowed permutations PERMUTATION_GROUP refers to the all allowed
    permutations of the atoms which are usually only rotations. Inversions are
    not chemically relevant and therefore not included in the permutations.

    :ivar atoms: Atoms
    :vartype atoms: tuple[int, ...]
    :ivar stereo: Stereochemistry
    :vartype stereo: Stereo
    """

    atoms: tuple[int, ...]
    parity: Optional[Literal[1, 0, -1]]
    PERMUTATION_GROUP: tuple[tuple[int, ...]]

    def __eq__(self, other: Any) -> bool: ...
    def __hash__(self) -> int: ...
    def get_isomers(self) -> tuple[Self]:
        """Returns all possible isomers of the stereochemistry"""


class _BaseStereo(ABC):
    __slots__ = ("atoms", "parity")

    atoms: tuple[int, ...]

    @property
    @abstractmethod
    def PERMUTATION_GROUP(self) -> tuple[tuple[int, ...], ...]: ...

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.atoms}, {self.parity})"

    def __init__(self, atoms: tuple[int, ...], parity: None | Literal[1, 0, -1] = None):
        self.atoms = atoms
        self.parity = parity

    def _perm_atoms(self) -> Iterator[tuple[int, ...]]:
        if self.parity is None:
            return (
                tuple([self.atoms[i] for i in perm])
                for perm in itertools.permutations(range(len(self.atoms)))
            )
        else:
            return (
                tuple([self.atoms[i] for i in perm]) for perm in self.PERMUTATION_GROUP
            )

    def get_isomers(self) -> tuple[Self]:
        """Returns all possible isomers of the stereochemistry"""
        return (self,)


class NoStereo(_BaseStereo):
    __slots__ = ()
    atoms: tuple[int, ...]
    parity: None

    def __init__(self, atoms: tuple[int, ...], parity: None = None):
        super().__init__(atoms=atoms, parity=parity)
        if parity is not None:
            raise ValueError("No stereo defined for this class")

    @property
    def PERMUTATION_GROUP(self) -> tuple[tuple[int, ...], ...]:
        return tuple(itertools.permutations(range(len(self.atoms))))


class _BaseChiralStereo(_BaseStereo, ABC):
    __slots__ = ()
    parity: Optional[Literal[1, -1]]

    @staticmethod
    @abstractmethod
    def _invert_atoms(atoms: tuple[int, ...]) -> tuple[int, ...]: ...

    @abstractmethod
    def get_isomers(self) -> tuple[Self, Self]: ...

    def invert(self) -> Self:
        if self.parity is None:
            return self
        return self.__class__(self.atoms, self.parity * -1)

    def _inverted_atoms(self) -> tuple[int, ...]:
        if self.parity is None:
            return self.atoms
        else:
            return self._invert_atoms(self.atoms)

    def __eq__(self, other: Any) -> bool:
        if other.parity == 0:
            return False

        s_atoms, o_atoms = self.atoms, other.atoms
        set_s_atoms = set(s_atoms)
        set_o_atoms = set(o_atoms)

        if len(s_atoms) != len(o_atoms) or not set_s_atoms.issuperset(set_o_atoms):
            return False

        if self.parity is None or other.parity is None:
            if set_s_atoms == set_o_atoms:
                return True
            return False

        elif self.parity == other.parity:
            if o_atoms == s_atoms or any(o_atoms == p for p in self._perm_atoms()):
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
            {tuple([self.atoms[i] for i in perm]) for perm in self.PERMUTATION_GROUP}
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


class _BaseAchiralStereo(_BaseStereo):
    __slots__ = ()
    parity: Optional[Literal[0]]

    def __eq__(self, other: Any) -> bool:
        if other.parity in (1, -1):
            return False

        s_atoms, o_atoms = self.atoms, other.atoms
        set_s_atoms = set(s_atoms)
        set_o_atoms = set(o_atoms)

        if len(s_atoms) != len(o_atoms) or not set_s_atoms.issuperset(set_o_atoms):
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


class AtomStereo(Stereo):
    @property
    def central_atom(self) -> AtomId:
        return self.atoms[0]


class BondStereo(Stereo):
    @property
    def bond(self) -> tuple[int, int]:
        return tuple(sorted(self.atoms[2:4]))


class Tetrahedral(_BaseChiralStereo, AtomStereo):
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
    atoms: tuple[int, int, int, int, int]
    parity: None | Literal[1, -1]

    def __init__(
        self,
        atoms: tuple[int, int, int, int, int],
        parity: None | Literal[1, -1] = None,
    ):
        if len(atoms) != 5:
            raise ValueError("Tetrahedral stereochemistry has 5 atoms")
        super().__init__(atoms=atoms, parity=parity)

    def get_isomers(self) -> tuple[Self, Self]:
        return {
            Tetrahedral(atoms=self.atoms, parity=1),
            Tetrahedral(atoms=self.atoms, parity=-1),
        }

    @property
    def PERMUTATION_GROUP(self) -> tuple[tuple[int, int, int, int, int], ...]:
        return (
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
        )

    @staticmethod
    def _invert_atoms(atoms) -> tuple[int, ...]:
        return tuple([atoms[i] for i in (0, 2, 1, 3, 4)])

    @classmethod
    def from_coords(
        cls,
        atoms: tuple[int, int, int, int, int],
        _,  # central atom
        p1: np.ndarray,
        p2: np.ndarray,
        p3: np.ndarray,
        p4: np.ndarray,
    ) -> Self:
        """
        Creates the representation of a Tetrahedral Stereochemistry
        from the coordinates of the atoms.
        :param atoms: Atoms of the stereochemistry
        :type atoms: tuple[int, int, int, int]
        :param p1: coordinates of atom 1
        :type p1: np.ndarray
        :param p2: coordinates of atom 2
        :type p2: np.ndarray
        :param p3: coordinates of atom 3
        :type p3: np.ndarray
        :param p4: coordinates of atom 4
        :type p4: np.ndarray
        ...
        :return: Tetrahedral
        :rtype: Tetrahedral
        """
        orientation = handedness(p1, p2, p3, p4)
        return cls(atoms, orientation)


class SquarePlanar(_BaseAchiralStereo, AtomStereo):
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
    atoms: tuple[int, int, int, int, int]
    parity: None | Literal[0]

    @property
    def central_atom(self) -> AtomId:
        return self.atoms[0]

    def get_isomers(self) -> tuple[Self, ...]:
        return {
            SquarePlanar(atoms=(self.atoms[0], *perm), parity=0)
            for perm in itertools.permutations(self.atoms[1:])
        }

    @property
    def PERMUTATION_GROUP(self) -> tuple[tuple[int, int, int, int, int], ...]:
        return (
            (0, 1, 2, 3, 4),
            (0, 2, 3, 4, 1),
            (0, 3, 4, 1, 2),
            (0, 4, 1, 2, 3),
            (0, 4, 3, 2, 1),
            (0, 3, 2, 1, 4),
            (0, 2, 1, 4, 3),
            (0, 1, 4, 3, 2),
        )


class TrigonalBipyramidal(_BaseChiralStereo, AtomStereo):
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
    atoms: tuple[int, int, int, int, int, int]
    parity: None | Literal[1, -1]

    def get_isomers(self) -> tuple[Self, ...]:
        return {
            TrigonalBipyramidal(atoms=(self.atoms[0], *perm), parity=p)
            for perm in itertools.permutations(self.atoms[1:])
            for p in (1, -1)
        }

    @property
    def PERMUTATION_GROUP(self) -> tuple[tuple[int, int, int, int, int, int], ...]:
        return (
            (0, 1, 2, 3, 4, 5),
            (0, 1, 2, 5, 3, 4),
            (0, 1, 2, 4, 5, 3),
            (0, 2, 1, 3, 5, 4),
            (0, 2, 1, 5, 4, 3),
            (0, 2, 1, 4, 3, 5),
        )

    @staticmethod
    def _invert_atoms(
        atoms: tuple[int, int, int, int, int, int],
    ) -> tuple[int, int, int, int, int, int]:
        return tuple([atoms[i] for i in (0, 1, 2, 3, 5, 4)])

    @classmethod
    def from_coords(
        cls: type[Self],
        atoms: tuple[int, int, int, int, int, int],
        _,
        p1: np.ndarray,
        p2: np.ndarray,
        p3: np.ndarray,
        p4: np.ndarray,
        p5: np.ndarray,
    ) -> TrigonalBipyramidal:
        """
        calculates the distance of the atom 5 from the plane defined by the
        first three atoms in Angstrom. The sign of the distance is determined
        by the side of the plane that atom 5 is on.
        """
        # For a trigonal bipyramidal geometry there are three atoms equatorial,
        # one on the top and one on the bottom.
        # first all combinations of three atoms are generated.
        # If Tetrahedral.from_coords(a, b, c, d)
        # == Tetrahedral.from_coords(a,b,c,e).invert()
        # Than the atoms a, b, c are equatorial, because they have an atom
        # above and below them.
        # If four atoms are in one plane the structure is a tetragonal pyramid.

        # coords_dict
        central_atom, *atoms = atoms
        cd: dict[int, Any] = {
            i: coords for i, coords in zip(atoms, [p1, p2, p3, p4, p5])
        }

        for comb in itertools.combinations(iterable=cd.keys(), r=3):
            i, j = [x for x in set(cd.keys()) - set(comb)]
            if are_planar(cd[comb[0]], cd[comb[1]], cd[comb[2]], cd[i]) or are_planar(
                cd[comb[0]], cd[comb[1]], cd[comb[2]], cd[j]
            ):
                raise ValueError("four atoms are planar")

            i_rotation = handedness(cd[comb[0]], cd[comb[1]], cd[comb[2]], cd[i])
            j_rotation = handedness(cd[comb[0]], cd[comb[1]], cd[comb[2]], cd[j]) * -1

            comb_is_equatorial = i_rotation == j_rotation

            if comb_is_equatorial is True:
                atoms_in_new_order = (i, j, *comb)
                orientation = i_rotation

                if orientation == 1:
                    return TrigonalBipyramidal((central_atom, *atoms_in_new_order), 1)
                elif orientation == -1:
                    return TrigonalBipyramidal((central_atom, *atoms_in_new_order), -1)
        else:
            raise ValueError("something went wrong")


class Octahedral(_BaseChiralStereo, AtomStereo):
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
    atoms: tuple[int, int, int, int, int, int]
    parity: None | Literal[1, -1]

    def get_isomers(self) -> tuple[Self, ...]:
        return {
            Octahedral(atoms=(self.atoms[0], *perm), parity=p)
            for perm in itertools.permutations(self.atoms[1:])
            for p in (1, -1)
        }

    @property
    def PERMUTATION_GROUP(self) -> tuple[tuple[int, int, int, int, int, int, int], ...]:
        return (
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
        )

    @staticmethod
    def _invert_atoms(
        atoms: tuple[int, int, int, int, int, int],
    ) -> tuple[int, int, int, int, int, int]:
        return tuple([atoms[i] for i in (0, 2, 1, 3, 4, 5, 6)])


class PlanarBond(_BaseAchiralStereo, BondStereo):
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
    atoms: tuple[int, int, int, int, int, int]
    parity: None | Literal[0]

    def get_isomers(self) -> tuple[Self, ...]:
        return (
            PlanarBond(self.atoms, 0),
            PlanarBond(tuple(self.atoms[i] for i in (0, 1, 2, 3, 5, 4)), 0),
        )

    @property
    def PERMUTATION_GROUP(
        self,
    ) -> tuple[tuple[int, int, int, int, int, int], ...]:
        return (
            (0, 1, 2, 3, 4, 5),
            (1, 0, 2, 3, 5, 4),
            (4, 5, 3, 2, 0, 1),
            (5, 4, 3, 2, 1, 0),
        )

    def __init__(
        self,
        atoms: tuple[int, int, int, int, int, int],
        parity: Optional[Literal[0]] = None,
    ):
        if len(atoms) != 6:
            raise ValueError("PlanarBond needs 6 atoms")
        super().__init__(atoms, parity)

    @classmethod
    def from_coords(
        cls,
        atoms: tuple[int, int, int, int, int, int],
        p0: np.ndarray,
        p1: np.ndarray,
        p2: np.ndarray,
        p3: np.ndarray,
        p4: np.ndarray,
        p5: np.ndarray,
    ) -> Self:
        a = (p0 - p1) / np.linalg.norm(p0 - p1)
        b = (p4 - p5) / np.linalg.norm(p4 - p5)
        result = int(np.sign(np.dot(a, b)))
        if result == 1:
            return cls(atoms, 0)
        elif result == -1:
            atoms = tuple(atoms[i] for i in (1, 0, 2, 3, 4, 5))
            return cls(atoms, 0)
        elif result == 0:
            raise ValueError("atoms are tetrahedral")
        else:
            raise ValueError("something went wrong")


class AtropBond(_BaseChiralStereo, BondStereo):
    r"""
    Stereochemistry:::
        parity = 1          parity = -1
        1       5           1        5
         \     /            ◀      /
          2 - 3               2 - 3
        ◀      \            /      \
        0        4         0         4


    """

    def get_isomers(self) -> tuple[Self, ...]:
        return (
            AtropBond(self.atoms, 0),
            AtropBond(tuple(self.atoms[i] for i in (0, 1, 2, 3, 5, 4)), 0),
        )

    @property
    def PERMUTATION_GROUP(self) -> tuple[tuple[int, int, int, int, int, int], ...]:
        return (
            (0, 1, 2, 3, 4, 5),
            (1, 0, 2, 3, 5, 4),
            (4, 5, 3, 2, 1, 0),
            (5, 4, 3, 2, 0, 1),
        )

    @staticmethod
    def _invert_atoms(
        atoms: tuple[int, int, int, int, int, int],
    ) -> tuple[int, int, int, int, int, int]:
        return tuple([atoms[i] for i in (1, 0, 2, 3, 4, 5, 6)])
