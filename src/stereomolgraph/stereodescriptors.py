from __future__ import annotations

import itertools
from abc import ABC, abstractmethod
from collections import Counter
from typing import TYPE_CHECKING, Any, Generic, Literal, TypeVar

import numpy as np

from stereomolgraph.cartesian import are_planar, handedness

if TYPE_CHECKING:
    from collections.abc import Generator, Set

    from stereomolgraph.graphs.mg import AtomId

A = TypeVar("A", bound=tuple[int, ...], covariant=True)
P = TypeVar("P", bound=None | Literal[1, 0, -1], covariant=True)


class Stereo(ABC, Generic[A, P]):
    """
    Base Class to represent the orientation of a group of atoms in space and
    their allowed permutations PERMUTATION_GROUP refers to the all allowed
    permutations of the atoms which are usually only rotations. Inversions are
    not chemically relevant and therefore not included in the permutations.

    :ivar atoms: Atoms
    :vartype atoms: tuple[int, ...]
    :ivar stereo: Stereochemistry
    :vartype stereo: Stereo
    """

    atoms: A
    parity: P
    PERMUTATION_GROUP: frozenset[A]

    def __eq__(self, other: Any) -> bool: ...
    def __hash__(self) -> int: ...
    def get_isomers(self: Stereo[A, P]) -> Set[Stereo[A, P]]:
        ...
        # """Returns all possible isomers of the stereochemistry"""


class AtomStereo(Stereo[A, P], Generic[A, P]):
    @property
    def central_atom(self) -> AtomId:
        return self.atoms[0]


class BondStereo(Stereo[A, P], Generic[A, P]):
    @property
    def bond(self) -> tuple[int, int]:
        return tuple(sorted(self.atoms[2:4]))  # # type: ignore[return-value]


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

    @abstractmethod
    def get_isomers(self: Stereo[A, P]) -> Set[Stereo[A, P]]:
        ...
        # """Returns all possible isomers of the stereochemistry"""


class _ChiralStereoMixin(
    _StereoMixin[A, None | Literal[1, -1]], ABC, Generic[A]
):
    __slots__ = "inversion"
    inversion: A

    @abstractmethod
    def get_isomers(
        self: _ChiralStereoMixin[A],
    ) -> Set[_ChiralStereoMixin[A]]: ...

    def invert(self: _ChiralStereoMixin[A]) -> _ChiralStereoMixin[A]:
        if self.parity is None:
            return self
        return self.__class__(self.atoms, self.parity * -1)

    def _inverted_atoms(self) -> A:
        if self.parity is None:
            return self.atoms
        else:
            return tuple([self.atoms[i] for i in self.inversion])  # type: ignore[return-value]

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
    __slots__ = ()

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

    def get_isomers(self) -> set[Tetrahedral]:
        return {
            Tetrahedral(atoms=self.atoms, parity=1),
            Tetrahedral(atoms=self.atoms, parity=-1),
        }

    @classmethod
    def from_coords(
        cls,
        atoms: tuple[int, int, int, int, int],
        _,  # central atom
        p1: np.ndarray,
        p2: np.ndarray,
        p3: np.ndarray,
        p4: np.ndarray,
    ) -> Tetrahedral:
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
            SquarePlanar(atoms=(self.atoms[0], *perm), parity=0)  # type: ignore[arg-type]
            for perm in itertools.permutations(self.atoms[1:])
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

    def get_isomers(self) -> set[TrigonalBipyramidal]:
        return {
            TrigonalBipyramidal(atoms=(self.atoms[0], *perm), parity=p)  # type: ignore[arg-type]
            for perm in itertools.permutations(self.atoms[1:])
            for p in (1, -1)
        }

    @classmethod
    def from_coords(
        cls: type[TrigonalBipyramidal],
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
        central_atom, *outer_atoms = atoms
        cd: dict[int, Any] = {
            i: coords for i, coords in zip(outer_atoms, [p1, p2, p3, p4, p5])
        }

        for comb in itertools.combinations(iterable=cd.keys(), r=3):
            i, j = [x for x in set(cd.keys()) - set(comb)]
            if are_planar(
                cd[comb[0]], cd[comb[1]], cd[comb[2]], cd[i]
            ) or are_planar(cd[comb[0]], cd[comb[1]], cd[comb[2]], cd[j]):
                raise ValueError("four atoms are planar")

            i_rotation = handedness(
                cd[comb[0]], cd[comb[1]], cd[comb[2]], cd[i]
            )
            j_rotation = (
                handedness(cd[comb[0]], cd[comb[1]], cd[comb[2]], cd[j]) * -1
            )

            comb_is_equatorial = i_rotation == j_rotation

            if comb_is_equatorial is True:
                atoms_in_new_order = (i, j, *comb)
                orientation = i_rotation

                if orientation == 1:
                    return TrigonalBipyramidal(
                        (central_atom, *atoms_in_new_order), 1
                    )
                elif orientation == -1:
                    return TrigonalBipyramidal(
                        (central_atom, *atoms_in_new_order), -1
                    )
        else:
            raise ValueError("something went wrong")


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
            Octahedral(atoms=(self.atoms[0], *perm), parity=p)  # type: ignore[arg-type]
            for perm in itertools.permutations(self.atoms[1:])
            for p in (1, -1)
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
        return {
            PlanarBond(self.atoms, 0),
            PlanarBond(tuple(self.atoms[i] for i in (0, 1, 2, 3, 5, 4)), 0),  # type: ignore[arg-type]
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
        p0: np.ndarray,
        p1: np.ndarray,
        p2: np.ndarray,
        p3: np.ndarray,
        p4: np.ndarray,
        p5: np.ndarray,
    ) -> PlanarBond:
        a = (p0 - p1) / np.linalg.norm(p0 - p1)
        b = (p4 - p5) / np.linalg.norm(p4 - p5)
        result = int(np.sign(np.dot(a, b)))
        if result == 1:
            return cls(atoms, 0)
        elif result == -1:
            atoms = tuple(atoms[i] for i in (1, 0, 2, 3, 4, 5))  # type: ignore[arg-type]
            return cls(atoms, 0)
        elif result == 0:
            raise ValueError("atoms are tetrahedral")
        else:
            raise ValueError("something went wrong")


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
        return {
            AtropBond(self.atoms, 1),
            AtropBond(tuple(self.atoms[i] for i in (0, 1, 2, 3, 5, 4)), -1),  # type: ignore[arg-type]
        }
