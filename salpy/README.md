# Salpy: Satellite LiDAR Point Adjustment for Python

Salpy is a Python package for optimizing and adjusting satellite LiDAR point positions using genetic algorithms. It provides tools for linear adjustment and location correction of spatial point data, particularly focused on satellite-based LiDAR observations.

## Features

- Linear alignment of point data to best-fit lines
- Position adjustment with specified offsets
- Loss minimization using genetic algorithms
- Distance calculations between point data and reference surfaces
- Value extraction from raster data at point locations
- Positional correction based on reference data

## Installation

```bash
pip install salpy
```

## Usage

```python
import geopandas as gpd
from salpy import linear_alignment, position_adjustment, minimize_loss

# Example: Linear alignment
gdf = gpd.read_file('path_to_your_shapefile.shp')
crs_code = 3857
adjusted_gdf = linear_alignment(gdf, crs_code)

# Example: Position adjustment
adjusted_gdf = position_adjustment(gdf, add_x=10, add_y=5, crs_code=3857)

# Example: Minimize loss (optimization)
import rasterio
with rasterio.open('reference_dem.tif') as src:
    raster = src.read(1)
    transform = src.transform

result = minimize_loss(
    lidar_footprints=gdf,
    input_rast=raster,
    raster_transform=transform,
    minimizing_method="euclidean",
    target_variable="mean",
    buf=12.5,
    crs_code=3857,
    lidar_value="elevation",
    lower_bounds=[-30, -30],
    upper_bounds=[30, 30],
    pop_size=50,
    max_iter=100,
    parallel=False
)
```

## Dependencies

- numpy
- pandas
- geopandas
- rasterio
- scipy
- scikit-learn
- pyproj
- shapely
- dtw-python
- deap (for genetic algorithms)

## License

GPL-3.0

## Credits

This package is a Python port of the R package [salpa](https://github.com/naru-T/salpa) by Narumasa Tsutsumida.
