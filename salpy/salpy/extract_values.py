"""
Extract Values from Raster at Point Locations

This module provides functionality to extract values from a raster at point locations,
optionally using a buffer around each point.
"""

import numpy as np
import rasterio
from geopandas import GeoDataFrame
from rasterio.mask import mask
from shapely.geometry import mapping


def extract_values(
    gdf: GeoDataFrame,
    raster,
    raster_transform,
    buf: float = 0,
    stat: str = "mean",
    crs_code: int = None
) -> GeoDataFrame:
    """
    Extracts values from a raster at point locations, optionally using a buffer.

    Parameters
    ----------
    gdf : GeoDataFrame
        Input GeoDataFrame containing point geometries
    raster : numpy.ndarray
        Input raster data as a numpy array
    raster_transform : affine.Affine
        Affine transformation that maps pixel coordinates to coordinates in the CRS
    buf : float, optional
        Buffer size around each point, by default 0
    stat : str, optional
        Statistic to compute for buffered areas ('mean', 'median', 'min', 'max', 'sum'),
        by default "mean"
    crs_code : int, optional
        Coordinate reference system code, by default None

    Returns
    -------
    GeoDataFrame
        Input GeoDataFrame with an additional column containing the extracted values

    Examples
    --------
    >>> import geopandas as gpd
    >>> import rasterio
    >>> gdf = gpd.read_file('points.shp')
    >>> with rasterio.open('dem.tif') as src:
    ...     raster = src.read(1)
    ...     transform = src.transform
    >>> result = extract_values(gdf, raster, transform, buf=10, stat='mean')
    """
    # Ensure input is a GeoDataFrame with valid geometry
    if not isinstance(gdf, GeoDataFrame):
        raise TypeError("Input gdf must be a GeoDataFrame")

    # Transform to specified CRS if provided
    if crs_code is not None:
        gdf = gdf.to_crs(epsg=crs_code)

    # Create a copy of the input GeoDataFrame
    result_gdf = gdf.copy()

    # Create a column for the extracted values
    result_gdf['extracted_value'] = np.nan

    # Create a memory raster dataset
    with rasterio.io.MemoryFile() as memfile:
        with memfile.open(
            driver='GTiff',
            height=raster.shape[0],
            width=raster.shape[1],
            count=1,
            dtype=raster.dtype,
            transform=raster_transform,
            crs=f"EPSG:{gdf.crs.to_epsg()}" if crs_code is None else f"EPSG:{crs_code}"
        ) as dataset:
            dataset.write(raster, 1)

            # Extract values for each point
            for idx, row in result_gdf.iterrows():
                geom = row.geometry

                # Apply buffer if specified
                if buf > 0:
                    geom = geom.buffer(buf)

                # Extract values using mask
                try:
                    out_image, out_transform = mask(dataset, [mapping(geom)], crop=True)

                    # Filter out nodata values
                    valid_data = out_image[0][out_image[0] != dataset.nodata]

                    if len(valid_data) > 0:
                        # Calculate the requested statistic
                        if stat == 'mean':
                            result_gdf.at[idx, 'extracted_value'] = np.mean(valid_data)
                        elif stat == 'median':
                            result_gdf.at[idx, 'extracted_value'] = np.median(valid_data)
                        elif stat == 'min':
                            result_gdf.at[idx, 'extracted_value'] = np.min(valid_data)
                        elif stat == 'max':
                            result_gdf.at[idx, 'extracted_value'] = np.max(valid_data)
                        elif stat == 'sum':
                            result_gdf.at[idx, 'extracted_value'] = np.sum(valid_data)
                        else:
                            raise ValueError(f"Unsupported statistic: {stat}")
                except Exception as e:
                    print(f"Error extracting values for feature {idx}: {e}")

    return result_gdf
