"""
Linear Alignment of Points

This module provides functionality to align satellite LiDAR footprints to a best-fit line
using an analytical total least squares approach.
"""

import numpy as np
import pandas as pd
from geopandas import GeoDataFrame
from shapely.geometry import Point


def linear_alignment(gdf: GeoDataFrame, crs_code: int) -> GeoDataFrame:
    """
    Aligns satellite LiDAR footprints to a best-fit line using an analytical total least squares approach.

    The method computes the principal component of the coordinates for each shot group to minimize
    the perpendicular distances between the points and the line, and then projects the points onto this line.

    Parameters
    ----------
    gdf : GeoDataFrame
        A GeoDataFrame containing spatial features. It must include a column named 'shot_number'.
    crs_code : int
        An integer representing the coordinate reference system (CRS) code to which the GeoDataFrame will be transformed.

    Returns
    -------
    GeoDataFrame
        A GeoDataFrame with adjusted coordinates that best fit a line, including original coordinates
        stored in 'orig_x' and 'orig_y'.

    Examples
    --------
    >>> import geopandas as gpd
    >>> gdf = gpd.read_file('path_to_your_shapefile.shp')
    >>> crs_code = 3857
    >>> adjusted_gdf = linear_alignment(gdf, crs_code)
    """
    # Validate input
    if not isinstance(gdf, GeoDataFrame):
        raise TypeError("gdf must be a GeoDataFrame")

    if not isinstance(crs_code, int):
        raise TypeError("crs_code must be an integer EPSG code")

    if 'shot_number' not in gdf.columns:
        raise ValueError("The GeoDataFrame does not contain a 'shot_number' column, which is required for the function to work correctly.")

    # Transform to specified CRS
    gdf = gdf.to_crs(epsg=crs_code)

    # Extract unique shot groups (using the first 12 characters of shot_number)
    shot_ids = gdf['shot_number'].astype(str).str[:12].unique()

    adjusted_gdf_list = []

    for shot in shot_ids:
        # Subset the GeoDataFrame for the current shot group
        gdf_subset = gdf[gdf['shot_number'].astype(str).str[:12] == shot]

        # Extract coordinates
        coords = np.array([[geom.x, geom.y] for geom in gdf_subset.geometry])

        # If there's only one point, use the point as is
        if len(coords) == 1:
            gdf_subset['orig_x'] = coords[0, 0]
            gdf_subset['orig_y'] = coords[0, 1]
            gdf_subset['coords_x'] = coords[0, 0]
            gdf_subset['coords_y'] = coords[0, 1]
        else:
            # Compute the centroid
            centroid = np.mean(coords, axis=0)

            # Center the coordinates
            centered = coords - centroid

            # Compute the covariance matrix and perform eigen decomposition
            cov_mat = np.cov(centered, rowvar=False)
            eigenvalues, eigenvectors = np.linalg.eigh(cov_mat)

            # Get the principal eigenvector (direction of maximum variance)
            # The eigenvectors are sorted by eigenvalues in ascending order, so we take the last one
            direction = eigenvectors[:, -1]

            # Project each point onto the line passing through the centroid in the direction of 'direction'
            projections = np.zeros_like(coords)
            for i, pt in enumerate(coords):
                # Calculate the projection of (pt - centroid) onto the direction vector
                proj_scalar = np.dot(pt - centroid, direction)
                # Calculate the projected point
                projections[i] = centroid + proj_scalar * direction

            # Store original and projected coordinates
            gdf_subset['orig_x'] = coords[:, 0]
            gdf_subset['orig_y'] = coords[:, 1]
            gdf_subset['coords_x'] = projections[:, 0]
            gdf_subset['coords_y'] = projections[:, 1]

        # Create an adjusted GeoDataFrame with the new projected coordinates
        geometry = [Point(x, y) for x, y in zip(gdf_subset['coords_x'], gdf_subset['coords_y'])]
        adjusted_gdf = GeoDataFrame(
            gdf_subset.drop(columns='geometry'),
            geometry=geometry,
            crs=f"EPSG:{crs_code}"
        )

        # Append to the list
        adjusted_gdf_list.append(adjusted_gdf)

    # Merge all adjusted GeoDataFrames into one
    final_adjusted_gdf = pd.concat(adjusted_gdf_list, ignore_index=True)

    return final_adjusted_gdf
