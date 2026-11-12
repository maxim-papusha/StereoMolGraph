"""Regression tests for exact stereodescriptor equality and hash consistency."""

import pytest

from stereomolgraph.stereodescriptors import (
    AtropBond,
    Octahedral,
    PlanarBond,
    SquarePlanar,
    Tetrahedral,
    TrigonalBipyramidal,
)


@pytest.mark.parametrize(
    ("stereo_type", "atoms", "assigned_parity"),
    [
        (Tetrahedral, (6, 0, 1, 2, 3), 1),
        (Tetrahedral, (6, 0, 1, 2, 3), -1),
        (SquarePlanar, (6, 0, 1, 2, 3), 0),
        (TrigonalBipyramidal, (6, 0, 1, 2, 3, 4), 1),
        (Octahedral, (6, 0, 1, 2, 3, 4, 5), -1),
        (PlanarBond, (0, 1, 2, 3, 4, 5), 0),
        (AtropBond, (0, 1, 2, 3, 4, 5), 1),
    ],
)
def test_undefined_is_not_equal_to_assigned(stereo_type, atoms, assigned_parity):
    undefined = stereo_type(atoms, None)
    assigned = stereo_type(atoms, assigned_parity)

    assert undefined != assigned
    assert assigned != undefined
    assert len({undefined, assigned}) == 2


def test_undefined_equality_ignores_atom_order_but_preserves_multiplicity():
    first = Tetrahedral((6, 0, 1, 2, None), None)
    reordered = Tetrahedral((6, None, 2, 0, 1), None)
    different_multiplicity = Tetrahedral((6, 0, 0, 2, None), None)

    assert first == reordered
    assert hash(first) == hash(reordered)
    assert first != different_multiplicity
    assert len({first, reordered, different_multiplicity}) == 2


def test_undefined_stereo_requires_same_geometry_type():
    atoms = (6, 0, 1, 2, 3)
    tetrahedral = Tetrahedral(atoms, None)
    square_planar = SquarePlanar(atoms, None)

    assert tetrahedral != square_planar
    assert square_planar != tetrahedral
    assert len({tetrahedral, square_planar}) == 2


def test_undefined_stereo_does_not_bridge_distinct_assigned_configurations():
    atoms = (6, 0, 1, 2, 3)
    unknown = Tetrahedral(atoms, None)
    first = Tetrahedral(atoms, 1)
    second = Tetrahedral(atoms, -1)

    assert unknown != first
    assert unknown != second
    assert first != second
    assert len({unknown, first, second}) == 3
