from copy import deepcopy
from itertools import permutations

import numpy as np
import pytest

from stereomolgraph.coords import Geometry, are_planar
from stereomolgraph.stereodescriptors import (
    AtropBond,
    Octahedral,
    PlanarBond,
    RigidBond12,
    RigidBond13,
    RigidBond23,
    RigidBond33,
    SquarePlanar,
    Tetrahedral,
    TrigonalBipyramidal,
)
from stereomolgraph.xyz2graph import atom_stereo_from_coords


class TestPlanar:
    def test_are_planar_true(self):
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.1],
            ],
            dtype=np.float64,
        )
        assert are_planar(coords, threshold=0.5)

    def test_are_planar_false(self):
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 1.0],
            ],
            dtype=np.float64,
        )
        assert not are_planar(coords, threshold=0.1)


class TestTetrahedral:
    def test_from_coords(self, enantiomer_geos):
        coords1 = enantiomer_geos[0].coords
        coords2 = enantiomer_geos[1].coords
        atoms = (3, 0, 1, 2, 4)

        stereo1 = atom_stereo_from_coords(atoms, coords1.take(atoms, axis=0))
        stereo2 = atom_stereo_from_coords(atoms, coords2.take(atoms, axis=0))
        assert stereo1 is not None and stereo2 is not None
        assert stereo1.parity == -1
        assert stereo2.parity == 1

    def test_from_permuted_coords(self, enantiomer_geos):
        coords = enantiomer_geos[0].coords
        different_perms = {
            atom_stereo_from_coords(atoms := (3, *perm), coords.take(atoms, axis=0))
            for perm in permutations((0, 1, 2, 4))
        }
        assert len(different_perms) == 1

    def test_equality(self):
        stereo1 = Tetrahedral((6, 0, 1, 2, 3), 1)
        stereo2 = Tetrahedral((6, 1, 2, 0, 3), 1)
        stereo3 = Tetrahedral((6, 0, 2, 1, 3), -1)
        stereo4 = Tetrahedral((6, 0, 2, 1, 3), 1)
        assert stereo1 == stereo2 == stereo3 != stereo4
        assert hash(stereo1) == hash(stereo2) == hash(stereo3) != hash(stereo4)

    def test_equality_with_none(self):
        stereo1 = Tetrahedral((6, 0, 1, 2, 3), None)
        stereo2 = Tetrahedral((6, 1, 2, 3, 0), None)
        assert stereo1 == stereo2

    def test_permutations(self):
        stereo1 = Tetrahedral((6, 0, 1, 2, 3), 1)
        stereo2 = Tetrahedral((6, 1, 2, 0, 3), 1)
        assert set(stereo1._perm_atoms()) == set(stereo2._perm_atoms())

    def test_is_immutable(self):
        stereo = Tetrahedral((6, 0, 1, 2, 3), 1)

        with pytest.raises(AttributeError):
            stereo.parity = -1

        with pytest.raises(AttributeError):
            stereo.atoms = (6, 1, 0, 2, 3)


class TestTrigonalBipyramidal:
    def test_from_coords(self, data_path):
        pcl5 = Geometry.from_xyz_file(data_path / "PCl5.xyz")
        result = atom_stereo_from_coords(
            (0, 1, 2, 3, 4, 5), pcl5.coords.take((0, 1, 2, 3, 4, 5), axis=0)
        )
        assert result is not None
        assert result.parity == 1
        assert set(result.atoms) == {3, 4, 0, 1, 2, 5}

    def test_equality(self):
        stereo1 = TrigonalBipyramidal((6, 0, 1, 2, 3, 4), 1)
        stereo2 = TrigonalBipyramidal((6, 1, 0, 2, 3, 4), -1)
        stereo3 = TrigonalBipyramidal((6, 1, 0, 3, 4, 2), -1)
        stereo4 = TrigonalBipyramidal((6, 1, 0, 3, 4, 2), 1)
        stereo5 = deepcopy(stereo4)
        assert stereo1 == stereo2 == stereo3 != stereo4
        assert (
            hash(stereo1)
            == hash(stereo2)
            == hash(stereo3)
            != hash(stereo4)
            == hash(stereo5)
        )

    def test_equality_with_none(self):
        stereo1 = TrigonalBipyramidal((0, 1, 2, 3, 4, 5), None)
        stereo2 = TrigonalBipyramidal((0, 2, 1, 3, 4, 5), None)
        different_center = TrigonalBipyramidal((1, 0, 2, 3, 4, 5), None)
        assert stereo1 == stereo2 != different_center


class TestPlanarBond:
    def test_equality(self):
        stereo1 = PlanarBond((5, 4, 3, 2, 1, 0), 0)
        stereo2 = PlanarBond((4, 5, 3, 2, 0, 1), 0)
        stereo3 = PlanarBond((1, 0, 2, 3, 5, 4), 0)
        stereo4 = PlanarBond((4, 5, 3, 2, 1, 0), 0)
        stereo5 = deepcopy(stereo4)
        assert stereo1 == stereo2 == stereo3 != stereo4 == stereo5
        assert (
            hash(stereo1)
            == hash(stereo2)
            == hash(stereo3)
            != hash(stereo4)
            == hash(stereo5)
        )

    def test_equality_with_none(self):
        stereo3 = PlanarBond((1, 0, 2, 3, 5, 4), None)
        stereo4 = PlanarBond((5, 4, 3, 2, 1, 0), None)
        assert stereo3 == stereo4


@pytest.mark.parametrize(
    "descriptor,parities,ligand_count",
    [
        (Tetrahedral, (1, -1), 4),
        (SquarePlanar, (0,), 4),
        (TrigonalBipyramidal, (1, -1), 5),
        (Octahedral, (1, -1), 6),
    ],
)
def test_undefined_atom_stereo_preserves_center(descriptor, parities, ligand_count):
    atoms = tuple(range(ligand_count + 1))
    unknown = descriptor(atoms, None)
    expected = {(0, *neighbors) for neighbors in permutations(atoms[1:])}
    assert set(unknown._perm_atoms()) == expected
    reordered = descriptor((0, *reversed(atoms[1:])), None)
    assert unknown == reordered
    assert hash(unknown) == hash(reordered)
    assert unknown.canonical_form() == reordered.canonical_form()
    assert unknown >= reordered and not unknown > reordered

    for reordered_atoms in expected:
        for parity in parities:
            assigned = descriptor(reordered_atoms, parity)
            assert unknown > assigned and unknown >= assigned
            assert assigned < unknown and assigned <= unknown
            assert not assigned >= unknown
            assert unknown != assigned
    changed_center = (1, 0, *atoms[2:])
    for parity in (None, *parities):
        different = descriptor(changed_center, parity)
        assert unknown != different
        assert not unknown >= different
        assert not different >= unknown


@pytest.mark.parametrize(
    "descriptor,parities,left_count,right_count",
    [
        (PlanarBond, (0,), 2, 2),
        (AtropBond, (1, -1), 2, 2),
        (RigidBond33, (1, -1), 3, 3),
        (RigidBond23, (1, -1), 2, 3),
        (RigidBond13, (1, -1), 1, 3),
        (RigidBond12, (0,), 1, 2),
    ],
)
def test_undefined_bond_stereo_preserves_connectivity(
    descriptor, parities, left_count, right_count
):
    left = tuple(range(left_count))
    a, b = left_count, left_count + 1
    right = tuple(range(b + 1, b + 1 + right_count))
    atoms = (*left, a, b, *right)
    unknown = descriptor(atoms, None)
    expected = {
        (*first, a, b, *second)
        for first in permutations(left)
        for second in permutations(right)
    }
    if left_count == right_count:
        expected.update(
            (*second, b, a, *first)
            for first in permutations(left)
            for second in permutations(right)
        )
    assert set(unknown._perm_atoms()) == expected
    for reordered_atoms in expected:
        equivalent = descriptor(reordered_atoms, None)
        assert unknown == equivalent
        assert hash(unknown) == hash(equivalent)
        assert unknown.canonical_form() == equivalent.canonical_form()
        assert unknown >= equivalent and not unknown > equivalent
        for parity in parities:
            assigned = descriptor(reordered_atoms, parity)
            assert unknown > assigned and assigned < unknown
            assert unknown >= assigned and assigned <= unknown
            assert unknown != assigned and not assigned >= unknown

    # Exchanging endpoints alone, moving a ligand to the other endpoint, or
    # exchanging a center with a ligand changes connectivity.
    for i, j in ((a, b), (0, len(atoms) - 1), (0, a)):
        changed = list(atoms)
        changed[i], changed[j] = changed[j], changed[i]
        for parity in (None, *parities):
            different = descriptor(tuple(changed), parity)
            assert unknown != different
            assert not unknown >= different
            assert not different >= unknown


@pytest.mark.parametrize(
    "descriptor,atoms,expected_count",
    [(PlanarBond, (0, 1, 2, 3, 4, 5), 2), (SquarePlanar, (0, 1, 2, 3, 4), 3)],
)
def test_undefined_includes_distinct_planar_isomers(descriptor, atoms, expected_count):
    unknown = descriptor(atoms, None)
    isomers = unknown.get_isomers()
    assert len(isomers) == expected_count
    for isomer in isomers:
        assert isomer.parity == 0
        assert unknown > isomer
        for other in isomers:
            assert (isomer >= other) == (isomer == other)
            assert not isomer > other


@pytest.mark.parametrize(
    "descriptor,atoms,equivalent_atoms,different_atoms,parity",
    [
        (Tetrahedral, (0, 1, 1, 2, None), (0, None, 1, 2, 1), (0, 1, 2, 2, None), 1),
        (
            PlanarBond,
            (None, 1, 2, 3, None, 4),
            (1, None, 2, 3, 4, None),
            (None, None, 2, 3, 1, 4),
            0,
        ),
        (
            RigidBond33,
            (None, None, 1, 3, 4, 5, 5, 6),
            (1, None, None, 3, 4, 6, 5, 5),
            (None, 1, 1, 3, 4, 5, 5, 6),
            1,
        ),
    ],
)
def test_undefined_preserves_optional_atom_multiplicities(
    descriptor, atoms, equivalent_atoms, different_atoms, parity
):
    unknown = descriptor(atoms, None)
    equivalent = descriptor(equivalent_atoms, None)
    assert unknown == equivalent
    assert hash(unknown) == hash(equivalent)
    assert unknown.canonical_form() == equivalent.canonical_form()
    assert unknown > descriptor(equivalent_atoms, parity)
    assert unknown != descriptor(different_atoms, None)
    assert not unknown >= descriptor(different_atoms, parity)
    assigned = descriptor(atoms, parity)
    permuted = descriptor(
        tuple(atoms[i] for i in descriptor.PERMUTATION_GROUP[-1]), parity
    )
    assert assigned == permuted
    assert assigned.canonical_form() == permuted.canonical_form()
