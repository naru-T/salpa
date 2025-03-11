"""
Calculate Distances Between Point Data and Reference Surfaces

This module provides functionality to calculate various distance metrics between
point data and reference surfaces.
"""

import dtw
import numpy as np
from geopandas import GeoDataFrame
from scipy.spatial.distance import cityblock, correlation, euclidean
from scipy.stats import hausdorff_distance

from .extract_values import extract_values


def perform_distance(
    gdf: GeoDataFrame,
    raster,
    raster_transform,
    method: str = "euclidean",
    buf: float = 0,
    target_variable: str = "mean",
    crs_code: int = None,
    value_column: str = None
) -> float:
    """
    Calculates various distance metrics between point data and reference surfaces.

    Parameters
    ----------
    gdf : GeoDataFrame
        Input GeoDataFrame containing point geometries
    raster : numpy.ndarray
        Input raster data as a numpy array
    raster_transform : affine.Affine
        Affine transformation that maps pixel coordinates to coordinates in the CRS
    method : str, optional
        Method for calculating the distance, by default "euclidean".
        Options are "dtw", "euclidean", "manhattan", "correlation", "area", "hausdorff".
    buf : float, optional
        Buffer size around each point, by default 0
    target_variable : str, optional
        Statistic to compute for buffered areas, by default "mean"
    crs_code : int, optional
        Coordinate reference system code, by default None
    value_column : str, optional
        Name of the column containing measurement values, by default None

    Returns
    -------
    float
        Calculated distance value

    Examples
    --------
    >>> import geopandas as gpd
    >>> import rasterio
    >>> gdf = gpd.read_file('points.shp')
    >>> with rasterio.open('dem.tif') as src:
    ...     raster = src.read(1)
    ...     transform = src.transform
    >>> distance = perform_distance(
    ...     gdf=gdf,
    ...     raster=raster,
    ...     raster_transform=transform,
    ...     method="euclidean",
    ...     buf=12.5,
    ...     target_variable="mean",
    ...     crs_code=3857,
    ...     value_column="elevation"
    ... )
    """
    # Ensure input is a GeoDataFrame with valid geometry
    if not isinstance(gdf, GeoDataFrame):
        raise TypeError("Input gdf must be a GeoDataFrame")

    # Check if value_column exists
    if value_column is not None and value_column not in gdf.columns:
        raise ValueError(f"Column '{value_column}' not found in gdf")

    # Extract values from raster
    extracted_gdf = extract_values(
        gdf,
        raster,
        raster_transform,
        buf=buf,
        stat=target_variable,
        crs_code=crs_code
    )

    # If no value_column is provided, return a default distance
    if value_column is None:
        return 0.0

    # Get point values and extracted raster values
    point_values = extracted_gdf[value_column].values
    raster_values = extracted_gdf['extracted_value'].values

    # Remove NaN values
    valid_indices = ~(np.isnan(point_values) | np.isnan(raster_values))
    point_values = point_values[valid_indices]
    raster_values = raster_values[valid_indices]

    # If no valid data, return a high distance value
    if len(point_values) == 0 or len(raster_values) == 0:
        return 1000.0

    # Calculate distance based on the specified method
    if method == "dtw":
        alignment = dtw.dtw(point_values, raster_values, keep_internals=True)
        distance_value = alignment.distance
    elif method == "euclidean":
        distance_value = euclidean(point_values, raster_values)
    elif method == "manhattan":
        distance_value = cityblock(point_values, raster_values)
    elif method == "correlation":
        # Higher correlation means lower distance, so we use 1 - correlation
        corr = correlation(point_values, raster_values)
        distance_value = 1 - corr if not np.isnan(corr) else 1000.0
    elif method == "area":
        # Calculate the area between the curves (simple approximation)
        distance_value = np.sum(np.abs(point_values - raster_values))
    elif method == "hausdorff":
        distance_value = hausdorff_distance(
            point_values.reshape(-1, 1),
            raster_values.reshape(-1, 1)
        )
    else:
        raise ValueError(f"Unsupported distance method: {method}")

    return distance_value
