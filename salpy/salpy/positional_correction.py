"""
Positional Correction of Satellite LiDAR Footprints

This module provides functionality to correct the positions of satellite LiDAR footprints
based on reference data using optimization techniques.
"""

import numpy as np
from geopandas import GeoDataFrame

from .minimize_loss import minimize_loss
from .position_adjustment import position_adjustment


def positional_correction(
    lidar_footprints: GeoDataFrame,
    reference_raster,
    raster_transform,
    minimizing_method: str = "euclidean",
    target_variable: str = "mean",
    buf: float = 12.5,
    crs_code: int = None,
    lidar_value: str = None,
    lower_bounds: list = None,
    upper_bounds: list = None,
    pop_size: int = 50,
    max_iter: int = 100,
    parallel: bool = False
) -> GeoDataFrame:
    """
    Corrects the positions of satellite LiDAR footprints based on reference data
    using optimization techniques.

    Parameters
    ----------
    lidar_footprints : GeoDataFrame
        Input satellite LiDAR footprints as a GeoDataFrame
    reference_raster : numpy.ndarray
        Reference raster data as a numpy array
    raster_transform : affine.Affine
        Affine transformation that maps pixel coordinates to coordinates in the CRS
    minimizing_method : str, optional
        Method for calculating the loss, by default "euclidean"
    target_variable : str, optional
        Statistic to compute for buffered areas, by default "mean"
    buf : float, optional
        Buffer size around each point, by default 12.5
    crs_code : int, optional
        Coordinate reference system code, by default None
    lidar_value : str, optional
        Name of the column containing LiDAR measurement values, by default None
    lower_bounds : list, optional
        Lower bounds for optimization parameters [x, y], by default [-30, -30]
    upper_bounds : list, optional
        Upper bounds for optimization parameters [x, y], by default [30, 30]
    pop_size : int, optional
        Population size for genetic algorithm, by default 50
    max_iter : int, optional
        Maximum number of iterations for genetic algorithm, by default 100
    parallel : bool, optional
        Whether to use parallel processing, by default False

    Returns
    -------
    GeoDataFrame
        Corrected satellite LiDAR footprints as a GeoDataFrame with updated geometry

    Examples
    --------
    >>> import geopandas as gpd
    >>> import rasterio
    >>> gdf = gpd.read_file('lidar_points.shp')
    >>> with rasterio.open('dem.tif') as src:
    ...     raster = src.read(1)
    ...     transform = src.transform
    >>> corrected_gdf = positional_correction(
    ...     lidar_footprints=gdf,
    ...     reference_raster=raster,
    ...     raster_transform=transform,
    ...     minimizing_method="euclidean",
    ...     target_variable="mean",
    ...     buf=12.5,
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

    # Run optimization to find the best offset
    optim_result = minimize_loss(
        lidar_footprints=lidar_footprints,
        input_rast=reference_raster,
        raster_transform=raster_transform,
        minimizing_method=minimizing_method,
        target_variable=target_variable,
        buf=buf,
        crs_code=crs_code,
        lidar_value=lidar_value,
        lower_bounds=lower_bounds,
        upper_bounds=upper_bounds,
        pop_size=pop_size,
        max_iter=max_iter,
        parallel=parallel
    )

    # If optimization failed, return the original data
    if optim_result is None:
        print("Optimization failed. Returning original data.")
        return lidar_footprints

    # Apply the optimal offset to the original data
    best_x = optim_result["best_x"]
    best_y = optim_result["best_y"]
    best_value = optim_result["best_value"]

    print(f"Optimal offset found: x={best_x}, y={best_y}, loss={best_value}")

    # Apply the position adjustment with the optimal offset
    corrected_footprints = position_adjustment(
        lidar_footprints,
        add_x=best_x,
        add_y=best_y,
        crs_code=crs_code if crs_code else lidar_footprints.crs.to_epsg()
    )

    return corrected_footprints
