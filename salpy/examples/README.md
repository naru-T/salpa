# Salpy Examples

This directory contains example scripts demonstrating how to use the Salpy package.

## Running Examples

To run an example, navigate to the examples directory and run:

```bash
python basic_usage.py
```

## Example Files

- `basic_usage.py`: Demonstrates basic usage of all main functions in the Salpy package

## Creating Your Own Examples

When creating your own examples, consider the following:

1. Import the necessary functions from the Salpy package
2. Create sample data or load real data
3. Apply the Salpy functions to the data
4. Visualize or analyze the results

Example:

```python
import geopandas as gpd
from salpy import linear_alignment

# Load data
gdf = gpd.read_file('path_to_your_shapefile.shp')

# Apply linear alignment
aligned_gdf = linear_alignment(gdf, crs_code=3857)

# Visualize results
import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
gdf.plot(ax=ax1, color='red')
ax1.set_title('Original Points')
aligned_gdf.plot(ax=ax2, color='blue')
ax2.set_title('Aligned Points')
plt.tight_layout()
plt.show()
