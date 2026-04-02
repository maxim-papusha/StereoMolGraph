from __future__ import annotations

from typing import Literal

from typing_extensions import Self

from stereomolgraph import Bond
from stereomolgraph.stereodescriptors import OInt, _StereoMixin

BondAutoCls = int  # bond automorphism class
MappingId = int


class NonRotatableBond33(
    _StereoMixin[
        tuple[OInt, OInt, OInt, int, int, OInt, OInt, OInt], None | Literal[1, -1]
    ],
):
    r"""
    Represents a bond that cannot freely rotate::
            parity = 1
             0    5
             |    |
        1  ▷ 3 - 4 ◁ 6
            ◀     ▶
           2        7

            parity = -1
             0    5
             |    |
        1  ▷ 3 - 4 ◁ 7
            ◀     ▶
           2        6
    """

    parity = 0
    inversion = (0, 2, 1, 3, 4, 5, 7, 6)
    _bond: Bond
    PERMUTATION_GROUP = (
        (0, 1, 2, 3, 4, 5, 6, 7),
        (5, 7, 6, 4, 3, 0, 2, 1),
        (1, 2, 0, 3, 4, 6, 7, 5),
        (6, 5, 7, 4, 3, 1, 0, 2),
        (2, 0, 1, 3, 4, 7, 5, 6),
        (7, 6, 5, 4, 3, 2, 1, 0),
    )

    def get_isomers(self) -> set[Self]:
        return {self}

    @property
    def bond(self) -> Bond:
        bond = frozenset(self.atoms[3:5])
        assert len(bond) == 2
        return bond


class NonRotatableBond23(
    _StereoMixin[tuple[OInt, OInt, int, int, OInt, OInt, OInt], None | Literal[1, -1]],
):
    r"""
    Represents a bond that cannot freely rotate::
            parity = 1
           0     4
            \    |
             2 - 3 ◁ 5
            /     ▶
           1        6

            parity = -1
           0     4
            \    |
             2 - 3 ◁ 6
            /     ▶
           1        5
    """

    parity = 0
    inversion = (0, 1, 2, 3, 4, 6, 5)
    _bond: Bond
    PERMUTATION_GROUP = ((0, 1, 2, 3, 4, 5, 6),)

    def get_isomers(self) -> set[Self]:
        return {self}

    @property
    def bond(self) -> Bond:
        bond = frozenset(self.atoms[2:4])
        assert len(bond) == 2
        return bond


class NonRotatableBond13(
    _StereoMixin[tuple[OInt, int, int, OInt, OInt, OInt], None | Literal[1, -1]],
):
    r"""
    Represents a bond that cannot freely rotate::
            parity = 1
          0     3
           \    |
            1 - 2 ◁ 4
                  ▶
                   5

            parity = -1
          0     3
           \    |
            1 - 2 ◁ 5
                  ▶
                   4
    """

    parity = 0
    inversion = (0, 1, 2, 3, 5, 4)
    _bond: Bond
    PERMUTATION_GROUP = ((0, 1, 2, 3, 4, 5),)

    def get_isomers(self) -> set[Self]:
        return {self}

    @property
    def bond(self) -> Bond:
        bond = frozenset(self.atoms[1:3])
        assert len(bond) == 2
        return bond
