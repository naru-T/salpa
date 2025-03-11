"""
Salpy: Satellite LiDAR Point Adjustment for Python
==================================================

A Python package for optimizing and adjusting satellite LiDAR point positions
using genetic algorithms.

Main Functions
-------------
linear_alignment : Aligns satellite LiDAR footprints to a best-fit line
position_adjustment : Adjusts point locations with specified offsets
minimize_loss : Optimizes point positions using genetic algorithms
perform_distance : Calculates distances between point data and reference surfaces
extract_values : Extracts values from raster data at point locations
positional_correction : Corrects positions based on reference data
"""

from .extract_values import extract_values
from .get_loss import get_loss
from .linear_alignment import linear_alignment
from .minimize_loss import minimize_loss
from .perform_distance import perform_distance
from .position_adjustment import position_adjustment
from .positional_correction import positional_correction

__version__ = '0.0.1'
__author__ = 'Narumasa Tsutsumida'
__email__ = 'rsnaru.jp@gmail.com'
