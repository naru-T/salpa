"""
Test for linear_alignment function
"""

import unittest
import numpy as np
from shapely.geometry import Point
import geopandas as gpd
from geopandas import GeoDataFrame
from salpy.linear_alignment import linear_alignment


class TestLinearAlignment(unittest.TestCase):
    """Test cases for linear_alignment function"""

    def setUp(self):
        """Set up test data"""
        # Create a simple GeoDataFrame with points that should form a line
        # but with some noise
        x = np.array([0, 1, 2, 3, 4])
        y = np.array([0, 1, 2, 3, 4])

        # Add some noise
        noise = np.random.normal(0, 0.1, 5)
        y_noisy = y + noise

        # Create shot numbers (all points belong to the same shot)
        shot_numbers = ['SHOT123456789_1'] * 5

        # Create geometry
        geometry = [Point(x[i], y_noisy[i]) for i in range(5)]

        # Create GeoDataFrame
        self.gdf = GeoDataFrame(
            {'shot_number': shot_numbers},
            geometry=geometry,
            crs="EPSG:4326"
        )

        # CRS code for testing
        self.crs_code = 4326

    def test_linear_alignment_basic(self):
        """Test that linear_alignment returns a GeoDataFrame"""
        result = linear_alignment(self.gdf, self.crs_code)
        self.assertIsInstance(result, GeoDataFrame)
        self.assertEqual(len(result), len(self.gdf))

    def test_linear_alignment_columns(self):
        """Test that linear_alignment adds the expected columns"""
        result = linear_alignment(self.gdf, self.crs_code)
        self.assertIn('orig_x', result.columns)
        self.assertIn('orig_y', result.columns)
        self.assertIn('coords_x', result.columns)
        self.assertIn('coords_y', result.columns)

    def test_linear_alignment_projection(self):
        """Test that points are projected onto a line"""
        result = linear_alignment(self.gdf, self.crs_code)

        # Get the projected coordinates
        coords_x = result['coords_x'].values
        coords_y = result['coords_y'].values

        # Calculate the slope and intercept of the line
        # For our test data, the line should be y = x
        slope = (coords_y[-1] - coords_y[0]) / (coords_x[-1] - coords_x[0])

        # The slope should be close to 1
        self.assertAlmostEqual(slope, 1.0, places=1)

        # Check that the points are more aligned than before
        # by comparing the variance of the distances to the line
        orig_x = result['orig_x'].values
        orig_y = result['orig_y'].values

        # Calculate distances from original points to the line y = x
        orig_distances = np.abs(orig_y - orig_x) / np.sqrt(2)

        # Calculate distances from projected points to the line y = x
        proj_distances = np.abs(coords_y - coords_x) / np.sqrt(2)

        # The variance of the projected distances should be smaller
        self.assertLess(np.var(proj_distances), np.var(orig_distances))


if __name__ == '__main__':
    unittest.main()
