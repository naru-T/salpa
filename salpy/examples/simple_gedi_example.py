"""
Simple example using real GEDI LiDAR data with salpy package
"""

import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from salpy import linear_alignment, position_adjustment

# Load GEDI data
gedi_file = "/Users/nt/Dropbox/Workplace/Github/salpa/inst/extdata/gedi_l2a_shot_sample.gpkg"
gedi_data = gpd.read_file(gedi_file)

# Use elev_lowestmode as our elevation column
elevation_column = "elev_lowestmode"

# Convert elevation column to float
gedi_data[elevation_column] = gedi_data[elevation_column].astype(float)

# Print basic information
print(f"GEDI data loaded: {len(gedi_data)} points")
print(f"CRS: {gedi_data.crs}")
print(f"Elevation range: {gedi_data[elevation_column].min():.2f} - {gedi_data[elevation_column].max():.2f} meters")

# Example 1: Linear alignment
aligned_data = linear_alignment(gedi_data, gedi_data.crs.to_epsg())
print(f"Original points: {len(gedi_data)}")
print(f"Aligned points: {len(aligned_data)}")

# Example 2: Position adjustment
adjusted_data = position_adjustment(gedi_data, add_x=10, add_y=5, crs_code=gedi_data.crs.to_epsg())
print(f"Adjusted {len(adjusted_data)} points with offset (10, 5)")

# Visualize results
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Plot original points
gedi_data.plot(ax=axs[0], column=elevation_column, cmap="viridis", legend=True, markersize=20)
axs[0].set_title("Original GEDI Points")

# Plot aligned points
aligned_data.plot(ax=axs[1], column=elevation_column, cmap="viridis", legend=True, markersize=20)
axs[1].set_title("Aligned GEDI Points")

# Plot adjusted points
adjusted_data.plot(ax=axs[2], column=elevation_column, cmap="viridis", legend=True, markersize=20)
axs[2].set_title("Adjusted GEDI Points")

plt.tight_layout()
output_file = "gedi_simple_example.png"
plt.savefig(output_file)
print(f"Visualization saved as {output_file}")