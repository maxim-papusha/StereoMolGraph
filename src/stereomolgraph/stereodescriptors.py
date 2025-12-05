from __future__ import annotations

import itertools
import sys
from collections import Counter
from typing import (
    TYPE_CHECKING,
    Any,
    Generic,
    Literal,
    Protocol,
    TypeVar,
    runtime_checkable,
)

if sys.version_info >= (3, 13):
    from typing import TypeVar
else:
    from typing_extensions import TypeVar

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Set

    from typing_extensions import Self

    from stereomolgraph.graphs.mg import AtomId, Bond


OInt = None | int
"Optional Integer"

A = TypeVar(
    "A", bound=tuple[OInt, ...],
    covariant=True,
    default=tuple[OInt, ...]
)
P = TypeVar(
    "P",
    covariant=True,
    bound=None | Literal[1, 0, -1],
    default=None | Literal[1, 0, -1],
)


@runtime_checkable
class Stereo(Protocol, Generic[A, P]):
    """
    Protocol to represent the orientation of a group of atoms in space.
    This is used to represent local stereochemistry and simultaneously the
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
    @property
    def PERMUTATION_GROUP(self,) -> Iterable[A]:
        """Defines all allowed permutations defined by the symmetry group under
        which the stereochemistry is invariant."""
        ...

    def __init__(self, atoms: A, parity: P = None): ...

    def __eq__(self, other: Any) -> bool: ...

    def __hash__(self) -> int: ...

    def invert(self) -> Self:
        """Inverts the stereo. If the stereo is achiral, it returns itself."""
        ...

    def get_isomers(self) -> Set[Self]:
        """Returns all stereoisomers of the stereochemistry. Not just the
        inverted ones, but all possible stereoisomers."""
        ...


@runtime_checkable
class AtomStereo(Stereo[A, P], Protocol, Generic[A, P]):
    @property
    def central_atom(self) -> AtomId: ...


@runtime_checkable
class BondStereo(Stereo[A, P], Protocol, Generic[A, P]):
    @property
    def bond(self) -> Bond: ...


class _StereoMixin(Generic[A, P]):
    PERMUTATION_GROUP: Iterable[A]
    inversion: None | A
    atoms: A
    parity: P

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.atoms}, {self.parity})"

    def __init__(self, atoms: A, parity: P = None):
        assert len(atoms) == len(self.PERMUTATION_GROUP[0])
        self.atoms = atoms
        self.parity = parity

    def _perm_atoms(self) -> Iterator[A]:
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

    def invert(self) -> Self:
        if self.parity is None:
            return self
        if self.parity == 0:
            return self
        new_parity = -self.parity
        assert new_parity in (1, -1)
        return self.__class__(self.atoms, new_parity) # type: ignore[return-value]

    def _inverted_atoms(self) -> A:
        if self.inversion is None:
            return self.atoms
        atoms = tuple([self.atoms[i] for i in self.inversion])
        assert len(atoms) == len(self.atoms) == len(self.inversion)
        return atoms  # type: ignore[return-value]

    def __eq__(self, other: Any) -> bool:
        if not hasattr(other, "atoms") or not hasattr(other, "parity"):
            return NotImplemented
        s_atoms, o_atoms = self.atoms, other.atoms
        set_s_atoms = set(s_atoms)
        set_o_atoms = set(o_atoms)

        if self.parity is None or other.parity is None:
            if set_s_atoms == set_o_atoms:
                return True
            return False

        if self.parity in (1, -1):
            if other.parity == 0:
                return False

            if len(s_atoms) != len(o_atoms) or not set_s_atoms.issuperset(
                set_o_atoms
            ):
                return False

            elif self.parity == other.parity:
                return o_atoms == s_atoms or any(
                    o_atoms == p for p in self._perm_atoms()
                )

            elif self.parity * -1 == other.parity:
                return any(
                    other._inverted_atoms() == p for p in self._perm_atoms()
                )

        if self.parity == 0:
            if other.parity in (1, -1):
                return False

            if len(s_atoms) != len(o_atoms) or not set_s_atoms.issuperset(
                set_o_atoms
            ):
                return False

            if other.parity is None:
                return set_s_atoms == set_o_atoms

            if self.parity == other.parity:
                return o_atoms == s_atoms or o_atoms in self._perm_atoms()
        raise RuntimeError(
            "This should not happen! "
            "Either the parity is not set correctly or the atoms are not "
            "ordered correctly."
        )

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
        # else parity in (1, -1):
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


class Tetrahedral(
    _StereoMixin[tuple[OInt, OInt, OInt, OInt, OInt], None | Literal[1, -1]],
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

    inversion = (0, 2, 1, 3, 4)
    PERMUTATION_GROUP = (
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

    def get_isomers(self) -> set[Self]:
        return {
            self.__class__(atoms=self.atoms, parity=1),
            self.__class__(atoms=self.atoms, parity=-1),
        }

    @property
    def central_atom(self) -> AtomId:
        return self.atoms[0]


class SquarePlanar(
    _StereoMixin[tuple[OInt, OInt, OInt, OInt, OInt], None | Literal[0]],
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

    inversion = None
    PERMUTATION_GROUP = (
            (0, 1, 2, 3, 4),
            (0, 2, 3, 4, 1),
            (0, 3, 4, 1, 2),
            (0, 4, 1, 2, 3),
            (0, 4, 3, 2, 1),
            (0, 3, 2, 1, 4),
            (0, 2, 1, 4, 3),
            (0, 1, 4, 3, 2),
    )

    def get_isomers(self) -> set[SquarePlanar]:
        return {
            SquarePlanar(atoms=atoms, parity=0)
            for perm in itertools.permutations(self.atoms[1:])
            if len(atoms := (self.atoms[0], *perm)) == 5
        }

    @property
    def central_atom(self) -> AtomId:
        return self.atoms[0]


class TrigonalBipyramidal(
    _StereoMixin[tuple[OInt, OInt, OInt, OInt, OInt, OInt], None | Literal[1, -1]],
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

    inversion = (0, 1, 2, 3, 5, 4)
    PERMUTATION_GROUP = (
            (0, 1, 2, 3, 4, 5),
            (0, 1, 2, 5, 3, 4),
            (0, 1, 2, 4, 5, 3),
            (0, 2, 1, 3, 5, 4),
            (0, 2, 1, 5, 4, 3),
            (0, 2, 1, 4, 3, 5),
    )

    def get_isomers(self) -> set[Self]:
        return {
            self.__class__(atoms=atoms, parity=p)
            for perm in itertools.permutations(self.atoms[1:])
            for p in (1, -1)
            if len(atoms := (self.atoms[0], *perm)) == 6
        }

    @property
    def central_atom(self) -> AtomId:
        return self.atoms[0]


class Octahedral(
    _StereoMixin[
        tuple[OInt, OInt, OInt, OInt, OInt, OInt, OInt], None | Literal[1, -1]
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

    inversion = (0, 2, 1, 3, 4, 5, 6)
    PERMUTATION_GROUP = (
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

    def get_isomers(self) -> set[Octahedral]:
        return {
            Octahedral(atoms=atoms, parity=p)
            for perm in itertools.permutations(self.atoms[1:])
            for p in (1, -1)
            if len((atoms := (self.atoms[0], *perm))) == 7
        }

    @property
    def central_atom(self) -> AtomId:
        return self.atoms[0]


class PlanarBond(
    _StereoMixin[tuple[OInt, OInt, OInt, OInt, OInt, OInt], None | Literal[0]],
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

    inversion = None
    PERMUTATION_GROUP = (
            (0, 1, 2, 3, 4, 5),
            (1, 0, 2, 3, 5, 4),
            (4, 5, 3, 2, 0, 1),
            (5, 4, 3, 2, 1, 0),
    )

    def get_isomers(self) -> set[PlanarBond]:
        reordered_atoms = tuple(self.atoms[i] for i in (0, 1, 2, 3, 5, 4))
        assert len(reordered_atoms) == 6
        return {
            PlanarBond(self.atoms, 0),
            PlanarBond(reordered_atoms, 0),
        }

    @property
    def bond(self) -> Bond:
        bond = frozenset(self.atoms[2:4])
        assert len(bond) == 2
        return bond


class AtropBond(
    _StereoMixin[tuple[OInt, OInt, OInt, OInt, OInt, OInt], None | Literal[1, -1]],
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

    inversion = (1, 0, 2, 3, 4, 5)
    PERMUTATION_GROUP = (
            (0, 1, 2, 3, 4, 5),
            (1, 0, 2, 3, 5, 4),
            (4, 5, 3, 2, 1, 0),
            (5, 4, 3, 2, 0, 1),
    )

    def get_isomers(self) -> set[AtropBond]:
        other_atoms = tuple(self.atoms[i] for i in (0, 1, 2, 3, 5, 4))
        assert len(other_atoms) == 6
        return {
            AtropBond(self.atoms, 1),
            AtropBond(other_atoms, -1),
        }

    @property
    def bond(self) -> Bond:
        bond = frozenset(self.atoms[2:4])
        assert len(bond) == 2
        return bond
