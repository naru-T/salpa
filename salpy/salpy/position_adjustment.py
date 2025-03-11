"""
Position Adjustment of Satellite LiDAR Footprints

This module provides functionality to adjust the point locations of satellite LiDAR footprints
with specified offsets.
"""

import numpy as np
from geopandas import GeoDataFrame
from shapely.geometry import Point


def position_adjustment(
    gdf: GeoDataFrame,
    add_x: float,
    add_y: float,
    crs_code: int
) -> GeoDataFrame:
    """
    Adjusts the point locations of satellite LiDAR footprints with specified offsets.

    Parameters
    ----------
    gdf : GeoDataFrame
        Input satellite LiDAR footprints as a GeoDataFrame
    add_x : float
        X-axis offset
    add_y : float
        Y-axis offset
    crs_code : int
        Coordinate reference system code

    Returns
    -------
    GeoDataFrame
        Adjusted satellite LiDAR footprints as a GeoDataFrame with updated geometry

    Examples
    --------
    >>> import geopandas as gpd
    >>> gdf = gpd.read_file('path_to_your_shapefile.shp')
    >>> crs_code = 3857
    >>> adjusted_gdf = position_adjustment(gdf, add_x=10, add_y=5, crs_code=crs_code)
    """
    # Ensure input is a GeoDataFrame with valid geometry
    if not isinstance(gdf, GeoDataFrame):
        raise TypeError("Input gdf must be a GeoDataFrame")

    if gdf.geometry.is_empty.any():
        raise ValueError("Input gdf contains empty geometries")

    # Validate CRS code
    if not isinstance(crs_code, int):
        raise TypeError("crs_code must be an integer EPSG code")

    # Transform to specified CRS
    gdf = gdf.to_crs(epsg=crs_code)

    # Get original coordinates
    coords = np.array([[geom.x, geom.y] for geom in gdf.geometry.centroid])

    # Store original coordinates
    gdf_copy = gdf.copy()
    gdf_copy['orig_x'] = coords[:, 0]
    gdf_copy['orig_y'] = coords[:, 1]

    # Create new coordinates with offsets
    new_coords_x = coords[:, 0] + add_x
    new_coords_y = coords[:, 1] + add_y

    # Create new geometry with adjusted coordinates
    new_geometry = [Point(x, y) for x, y in zip(new_coords_x, new_coords_y)]

    # Create new GeoDataFrame with adjusted coordinates
    adjusted_gdf = GeoDataFrame(
        gdf_copy.drop(columns='geometry'),
        geometry=new_geometry,
        crs=f"EPSG:{crs_code}"
    )

    return adjusted_gdf
