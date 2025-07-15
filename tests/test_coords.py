import tempfile
import textwrap
from pathlib import Path

import numpy as np
import pytest

from stereomolgraph.coords import Geometry, are_planar, angle_from_coords, handedness, pairwise_distances


class TestGeometry:
    @pytest.mark.parametrize(
        "xyz_content, comment, coords",
        [
            (  # with comment
                textwrap.dedent("""\
                    2
                        this is a comment comment with whitespace
                    C                 -3.7    0.02    0.2
                    C                 -3.1   -1.1   -0.2
                """),
                "this is a comment comment with whitespace",
                [[-3.7, 0.02, 0.2], [-3.1, -1.1, -0.2]],
            ),
            (  # with empty lines at end
                textwrap.dedent("""\
                    2

                    C                 -3.7    0.02    0.2
                    C                  0.0    0.0     0.0

                """),
                "",
                [[-3.7, 0.02, 0.2], [0.0, 0.0, 0.0]],
            ),
        ],
        ids=["with comment", "empty lines at end"],
    )
    def test_from_xyz_file(self, xyz_content, comment, coords):
        with tempfile.NamedTemporaryFile("w+", delete=False) as tmp:
            tmp.write(xyz_content)
            tmp.flush()
            fake_path = Path(tmp.name)

            geo = Geometry.from_xyz_file(fake_path)

        np.testing.assert_equal(geo.coords, coords)




class TestGeometryFunctions:
    """Test suite for geometry functions"""
    
    # Test data
    PLANAR_POINTS = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0.5, 0.5, 0]  # Perfectly planar
    ], dtype=np.float64) * 100
    
    NON_PLANAR_POINTS = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0.5, 0.5, 0.1]  # Slightly out of plane
    ], dtype=np.float64) * 100
    
    RIGHT_ANGLE_POINTS = np.array([
        [1, 0, 0],  # A
        [0, 0, 0],  # B (vertex)
        [0, 1, 0]   # C
    ], dtype=np.float64) * 100
    
    STRAIGHT_ANGLE_POINTS = np.array([
        [1, 0, 0],
        [0, 0, 0],
        [-1, 0, 0]
    ], dtype=np.float64) * 100
    
    CHIRAL_TETRAHEDRON = np.array([
        [0, 0, 0],  # Central atom
        [1, 0, 0],  # Right
        [0, 1, 0],  # Front
        [0, 0, 1]   # Up (R configuration)
    ], dtype=np.float64) * 100
    
    CHIRAL_TETRAHEDRON_S = np.array([
        [0, 0, 0],  # Central atom
        [1, 0, 0],  # Right
        [0, 1, 0],  # Front
        [0, 0, -1]  # Down (S configuration)
    ], dtype=np.float64) * 100
    
    def test_are_planar_non_vectorized(self):
        """Test planar check with single point set"""
        assert are_planar(self.PLANAR_POINTS)
        assert not are_planar(self.NON_PLANAR_POINTS)
        assert are_planar(self.NON_PLANAR_POINTS/100, threshold=1)  # Should pass with larger threshold
        
    @pytest.mark.skip("Not implemented yet")
    def test_are_planar_vectorized(self):
        """Test planar check with multiple point sets"""
        multiple_sets = np.stack([self.PLANAR_POINTS, self.NON_PLANAR_POINTS])
        result = are_planar(multiple_sets)
        assert result.shape == (2,)
        assert result[0] is np.bool_(True)
        assert result[1] is np.bool_(False)
        
    @pytest.mark.skip("Not implemented yet")
    def test_are_planar_edge_cases(self):
        with pytest.raises(AssertionError):
            assert are_planar(self.PLANAR_POINTS[:3])
        
        # Test with colinear points (degenerate case)
        colinear = np.array([
            [0, 0, 0],
            [1, 0, 0],
            [2, 0, 0],
            [3, 0, 0]
        ])
        assert are_planar(colinear)  # Technically planar
        
    def test_angle_from_coords_non_vectorized(self):
        """Test angle calculation with single point set"""
        angle = angle_from_coords(self.RIGHT_ANGLE_POINTS)
        assert np.isclose(angle, 90.0)
        
        angle = angle_from_coords(self.STRAIGHT_ANGLE_POINTS)
        assert np.isclose(angle, 180.0)
        
        # Test degenerate case (zero-length vectors) not included
        degenerate = np.array([
            [0, 0, 0],
            [0, 0, 0],  # All points same
            [0, 0, 0]
        ], dtype=np.float64)
        angle = angle_from_coords(degenerate)
        #assert np.isnan(angle)
        
    def test_angle_from_coords_vectorized(self):
        """Test angle calculation with multiple point sets"""
        multiple_angles = np.stack([
            self.RIGHT_ANGLE_POINTS,
            self.STRAIGHT_ANGLE_POINTS,
            np.array([[1,1,0], [0,0,0], [-1,1,0]])  # 90° angle
        ])
        
        angles = angle_from_coords(multiple_angles)
        assert angles.shape == (3,)
        assert np.allclose(angles, [90.0, 180.0, 90.0])
        
    def test_handedness_non_vectorized(self):
        """Test chirality calculation with single point set"""
        # Test R configuration
        result = handedness(self.CHIRAL_TETRAHEDRON)
        assert result == -1
        
        # Test S configuration
        result = handedness(self.CHIRAL_TETRAHEDRON_S)
        assert result == 1
        
        # Test planar case (should be 0)
        planar = np.vstack([self.PLANAR_POINTS[:3], [0.5, 0.5, 0]])
        result = handedness(planar)
        assert result == 0
        
    def test_handedness_vectorized(self):
        """Test chirality calculation with multiple point sets"""
        multiple_configs = np.stack([
            self.CHIRAL_TETRAHEDRON,
            self.CHIRAL_TETRAHEDRON_S,
            np.vstack([self.PLANAR_POINTS[:3], [0.5, 0.5, 0]])  # Planar
        ])
        
        results = handedness(multiple_configs)
        assert results.shape == (3,)
        assert np.array_equal(results, [-1, 1, 0])
        
    def test_pairwise_distances_non_vectorized(self):
        """Test distance matrix calculation"""
        points = np.array([
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0]
        ], dtype=np.float64)
        
        dist_matrix = pairwise_distances(points)
        expected = np.array([
            [0, 1, 1],
            [1, 0, np.sqrt(2)],
            [1, np.sqrt(2), 0]
        ])
        
        assert dist_matrix.shape == (3, 3)
        assert np.allclose(dist_matrix, expected)
        
    def test_pairwise_distances_vectorized(self):
        """Test batched distance matrix calculation"""
        batch_points = np.stack([
            self.PLANAR_POINTS,
            self.PLANAR_POINTS,
            self.PLANAR_POINTS,
            self.PLANAR_POINTS,
            self.PLANAR_POINTS
        ])
        result = pairwise_distances(self.PLANAR_POINTS)
        vec_results = pairwise_distances(batch_points)

        assert vec_results.shape == (5, 4, 4)
        for r in vec_results:
            assert r.shape == (4, 4)
            assert np.allclose(r, result)

    @pytest.mark.skip(reason="NaN handling not implemented yet")
    def test_all_functions_nan_handling(self):
        """Test how functions handle NaN inputs"""
        nan_points = np.full((4, 3), np.nan)
        
        with pytest.raises(ValueError):
            are_planar(nan_points)
            
        with pytest.raises(ValueError):
            handedness(nan_points[:4])
            
        with pytest.raises(ValueError):
            angle_from_coords(nan_points[:3])
            
        with pytest.raises(ValueError):
            pairwise_distances(nan_points)