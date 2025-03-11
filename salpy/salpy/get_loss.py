"""
Calculate Loss Between LiDAR and Reference Data

This module provides functionality to calculate the loss (difference) between
LiDAR point measurements and reference raster data.
"""

import dtw
import numpy as np
from geopandas import GeoDataFrame
from scipy.spatial.distance import cityblock, correlation, euclidean, directed_hausdorff

from .extract_values import extract_values
from .position_adjustment import position_adjustment


# Custom Hausdorff distance function since scipy.stats.hausdorff_distance is not available
def hausdorff_distance(u, v):
    """
    Compute the Hausdorff distance between two arrays.

    Parameters
    ----------
    u : array_like
        Input array.
    v : array_like
        Input array.

    Returns
    -------
    float
        The Hausdorff distance between arrays.
    """
    # Use scipy's directed_hausdorff function
    forward, _, _ = directed_hausdorff(u, v)
    backward, _, _ = directed_hausdorff(v, u)
    return max(forward, backward)


def get_loss(
    lidar_footprints: GeoDataFrame,
    add_x: float,
    add_y: float,
    buf: float,
    input_rast,
    raster_transform,
    minimizing_method: str = "euclidean",
    target_variable: str = "mean",
    crs_code: int = None,
    lidar_value: str = None
) -> dict:
    """
    Calculates the loss (difference) between LiDAR point measurements and reference raster data.

    Parameters
    ----------
    lidar_footprints : GeoDataFrame
        Input satellite LiDAR footprints as a GeoDataFrame
    add_x : float
        X-axis offset
    add_y : float
        Y-axis offset
    buf : float
        Buffer size around each point
    input_rast : numpy.ndarray
        Input raster data as a numpy array
    raster_transform : affine.Affine
        Affine transformation that maps pixel coordinates to coordinates in the CRS
    minimizing_method : str, optional
        Method for calculating the loss, by default "euclidean".
        Options are "dtw", "euclidean", "manhattan", "correlation", "area", "hausdorff".
    target_variable : str, optional
        Statistic to compute for buffered areas, by default "mean"
    crs_code : int, optional
        Coordinate reference system code, by default None
    lidar_value : str, optional
        Name of the column containing LiDAR measurement values, by default None

    Returns
    -------
    dict
        Dictionary containing the loss value and adjusted footprints

    Examples
    --------
    >>> import geopandas as gpd
    >>> import rasterio
    >>> gdf = gpd.read_file('lidar_points.shp')
    >>> with rasterio.open('dem.tif') as src:
    ...     raster = src.read(1)
    ...     transform = src.transform
    >>> result = get_loss(
    ...     lidar_footprints=gdf,
    ...     add_x=10,
    ...     add_y=5,
    ...     buf=12.5,
    ...     input_rast=raster,
    ...     raster_transform=transform,
    ...     minimizing_method="euclidean",
    ...     target_variable="mean",
    ...     crs_code=3857,
    ...     lidar_value="elevation"
    ... )
    """
    # Ensure input is a GeoDataFrame with valid geometry
    if not isinstance(lidar_footprints, GeoDataFrame):
        raise TypeError("Input lidar_footprints must be a GeoDataFrame")

    # Check if lidar_value column exists
    if lidar_value is not None and lidar_value not in lidar_footprints.columns:
        raise ValueError(f"Column '{lidar_value}' not found in lidar_footprints")

    # Adjust positions
    adjusted_footprints = position_adjustment(
        lidar_footprints, add_x, add_y, crs_code if crs_code else lidar_footprints.crs.to_epsg()
    )

    # Extract values from raster
    extracted_footprints = extract_values(
        adjusted_footprints,
        input_rast,
        raster_transform,
        buf=buf,
        stat=target_variable,
        crs_code=crs_code
    )

    # If no lidar_value is provided, return a default loss
    if lidar_value is None:
        return {
            "loss_value": 0,
            "adjusted_footprints": extracted_footprints
        }

    # Get LiDAR values and extracted raster values
    lidar_values = extracted_footprints[lidar_value].values
    raster_values = extracted_footprints['extracted_value'].values

    # Remove NaN values
    valid_indices = ~(np.isnan(lidar_values) | np.isnan(raster_values))
    lidar_values = lidar_values[valid_indices]
    raster_values = raster_values[valid_indices]

    # If no valid data, return a high loss value
    if len(lidar_values) == 0 or len(raster_values) == 0:
        return {
            "loss_value": 1000,
            "adjusted_footprints": extracted_footprints
        }

    # Calculate loss based on the specified method
    if minimizing_method == "dtw":
        alignment = dtw.dtw(lidar_values, raster_values, keep_internals=True)
        loss_value = alignment.distance
    elif minimizing_method == "euclidean":
        loss_value = euclidean(lidar_values, raster_values)
    elif minimizing_method == "manhattan":
        loss_value = cityblock(lidar_values, raster_values)
    elif minimizing_method == "correlation":
        # Higher correlation means lower loss, so we use 1 - correlation
        corr = correlation(lidar_values, raster_values)
        loss_value = 1 - corr if not np.isnan(corr) else 1000
    elif minimizing_method == "area":
        # Calculate the area between the curves (simple approximation)
        loss_value = np.sum(np.abs(lidar_values - raster_values))
    elif minimizing_method == "hausdorff":
        loss_value = hausdorff_distance(
            lidar_values.reshape(-1, 1),
            raster_values.reshape(-1, 1)
        )
    else:
        raise ValueError(f"Unsupported minimizing method: {minimizing_method}")

    return {
        "loss_value": loss_value,
        "adjusted_footprints": extracted_footprints
    }
