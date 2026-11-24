from copy import deepcopy
from itertools import permutations

import numpy as np
import pytest

from stereomolgraph.coords import Geometry, are_planar
from stereomolgraph.stereodescriptors import (
    AtropBond,
    PlanarBond,
    Tetrahedral,
    TrigonalBipyramidal,
)
from stereomolgraph.xyz2graph import atom_stereo_from_coords


@pytest.mark.parametrize(
    "stereo, expected",
    [
        (Tetrahedral((9, 4, 3, 2, 1), 1), (1, 9, 1, 2, 3, 4)),
        (Tetrahedral((9, 4, 3, 2, 1), -1), (1, 9, 1, 2, 4, 3)),
        (Tetrahedral((9, 4, 3, 2, 1), None), (2, 1, 2, 3, 4, 9)),
        (Tetrahedral((9, 3, 4, 2, 1), -1), (1, 9, 1, 2, 3, 4)),
    ],
)
def test_canonical_form_uses_smallest_order_and_normalizes_parity(stereo, expected):
    assert stereo.canonical_form() == expected


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
        stereo2 = TrigonalBipyramidal((1, 0, 2, 3, 4, 5), None)
        assert stereo1 == stereo2


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


class TestAtropBond:
    @staticmethod
    def helicity(atoms, coords):
        """Sign of the triple product (x0 - x1, x4 - x5, x3 - x2): a pseudoscalar,
        unchanged by rotations and flipped by a reflection."""
        x = {a: np.asarray(coords[a], dtype=float) for a in atoms}
        a0, a1, a2, a3, a4, a5 = atoms
        return np.sign(np.dot(x[a0] - x[a1], np.cross(x[a4] - x[a5], x[a3] - x[a2])))

    # A twisted axis 2-3 along z: 0, 1 on 2 in the xz plane, 4, 5 on 3 in the
    # yz plane (an allene, or a biaryl held at 90 degrees).
    coords = {
        0: (1.0, 0.0, -1.5),
        1: (-1.0, 0.0, -1.5),
        2: (0.0, 0.0, -0.7),
        3: (0.0, 0.0, 0.7),
        4: (0.0, 1.0, 1.5),
        5: (0.0, -1.0, 1.5),
    }

    def test_permutation_group_keeps_the_helicity(self):
        atoms = (0, 1, 2, 3, 4, 5)
        ref = self.helicity(atoms, self.coords)
        for perm in AtropBond.PERMUTATION_GROUP:
            permuted = tuple(atoms[i] for i in perm)
            assert self.helicity(permuted, self.coords) == ref, perm
        inverted = tuple(atoms[i] for i in AtropBond.inversion)
        assert self.helicity(inverted, self.coords) == -ref

    def test_equality(self):
        stereo = AtropBond((0, 1, 2, 3, 4, 5), 1)
        c2_image = AtropBond((5, 4, 3, 2, 1, 0), 1)
        mirror_image = AtropBond((4, 5, 3, 2, 1, 0), 1)
        assert stereo == c2_image
        assert hash(stereo) == hash(c2_image)
        assert stereo != mirror_image
        assert stereo == mirror_image.invert()

    def test_twisted_c2_symmetric_graph_differs_from_its_mirror_image(self):
        from stereomolgraph import StereoMolGraph

        def twisted_butene(parity):
            g = StereoMolGraph()
            for i, s in enumerate("CCCCHH"):
                g.add_atom(i, s)
            for b in [(0, 1), (1, 2), (2, 3), (1, 4), (2, 5)]:
                g.add_bond(*b)
            g.set_bond_stereo(AtropBond((0, 4, 1, 2, 3, 5), parity))
            return g

        assert twisted_butene(1) != twisted_butene(-1)

    def test_get_isomers_gives_both_helicities(self):
        stereo = AtropBond((0, 1, 2, 3, 4, 5), None)
        isomers = stereo.get_isomers()
        assert isomers == {
            AtropBond((0, 1, 2, 3, 4, 5), 1),
            AtropBond((0, 1, 2, 3, 4, 5), -1),
        }
        assert len(isomers) == 2
