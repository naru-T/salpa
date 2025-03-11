"""
Basic usage examples for the salpy package
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
import rasterio
from rasterio.transform import Affine
import matplotlib.pyplot as plt

from salpy import (
    linear_alignment,
    position_adjustment,
    extract_values,
    perform_distance,
    minimize_loss,
    positional_correction
)

# Create a sample GeoDataFrame with LiDAR points
def create_sample_data():
    """Create sample LiDAR point data"""
    # Create a simple grid of points
    x = np.linspace(0, 100, 10)
    y = np.linspace(0, 100, 10)
    xx, yy = np.meshgrid(x, y)

    # Flatten the grid
    x_flat = xx.flatten()
    y_flat = yy.flatten()

    # Add some random elevation values
    elevation = 100 + 0.5 * x_flat + 0.3 * y_flat + np.random.normal(0, 5, len(x_flat))

    # Create shot numbers (group points into shots)
    shot_numbers = []
    for i in range(len(x_flat)):
        shot_group = i // 10  # 10 points per shot
        shot_numbers.append(f"SHOT{shot_group:012d}_{i % 10}")

    # Create geometry
    geometry = [Point(x, y) for x, y in zip(x_flat, y_flat)]

    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(
        {
            'shot_number': shot_numbers,
            'elevation': elevation
        },
        geometry=geometry,
        crs="EPSG:4326"
    )

    return gdf

# Create a sample raster
def create_sample_raster():
    """Create a sample raster DEM"""
    # Create a simple elevation model
    x = np.linspace(0, 100, 100)
    y = np.linspace(0, 100, 100)
    xx, yy = np.meshgrid(x, y)

    # Create a surface with a slope and some random noise
    dem = 100 + 0.5 * xx + 0.3 * yy + np.random.normal(0, 2, (100, 100))

    # Create an affine transform
    transform = Affine(1.0, 0.0, 0.0, 0.0, 1.0, 0.0)

    return dem, transform

def main():
    """Run examples of salpy functionality"""
    # Create sample data
    print("Creating sample data...")
    lidar_points = create_sample_data()
    dem, transform = create_sample_raster()

    # Example 1: Linear alignment
    print("\nExample 1: Linear alignment")
    aligned_points = linear_alignment(lidar_points, 4326)
    print(f"Original points: {len(lidar_points)}")
    print(f"Aligned points: {len(aligned_points)}")
    print("Columns added:", set(aligned_points.columns) - set(lidar_points.columns))

    # Example 2: Position adjustment
    print("\nExample 2: Position adjustment")
    adjusted_points = position_adjustment(lidar_points, add_x=5, add_y=10, crs_code=4326)
    print(f"Adjusted {len(adjusted_points)} points with offset (5, 10)")

    # Example 3: Extract values from raster
    print("\nExample 3: Extract values from raster")
    extracted_points = extract_values(
        lidar_points,
        dem,
        transform,
        buf=2,
        stat="mean",
        crs_code=4326
    )
    print("Extracted values statistics:")
    print(extracted_points['extracted_value'].describe())

    # Example 4: Calculate distance
    print("\nExample 4: Calculate distance")
    distance = perform_distance(
        lidar_points,
        dem,
        transform,
        method="euclidean",
        buf=2,
        target_variable="mean",
        crs_code=4326,
        value_column="elevation"
    )
    print(f"Distance between LiDAR and DEM: {distance:.2f}")

    # Example 5: Minimize loss (optimization)
    print("\nExample 5: Minimize loss (optimization)")
    result = minimize_loss(
        lidar_footprints=lidar_points,
        input_rast=dem,
        raster_transform=transform,
        minimizing_method="euclidean",
        target_variable="mean",
        buf=2,
        crs_code=4326,
        lidar_value="elevation",
        lower_bounds=[-20, -20],
        upper_bounds=[20, 20],
        pop_size=20,
        max_iter=10,  # Small number for quick example
        parallel=False
    )
    if result:
        print(f"Optimal offset: x={result['best_x']:.2f}, y={result['best_y']:.2f}")
        print(f"Optimal loss value: {result['best_value']:.2f}")

    # Example 6: Positional correction
    print("\nExample 6: Positional correction")
    corrected_points = positional_correction(
        lidar_footprints=lidar_points,
        reference_raster=dem,
        raster_transform=transform,
        minimizing_method="euclidean",
        target_variable="mean",
        buf=2,
        crs_code=4326,
        lidar_value="elevation",
        lower_bounds=[-20, -20],
        upper_bounds=[20, 20],
        pop_size=20,
        max_iter=10,  # Small number for quick example
        parallel=False
    )
    print(f"Corrected {len(corrected_points)} points")

    # Visualize results
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))

    # Plot original points
    lidar_points.plot(ax=axs[0, 0], column='elevation', cmap='viridis',
                     legend=True, markersize=20)
    axs[0, 0].set_title('Original LiDAR Points')

    # Plot aligned points
    aligned_points.plot(ax=axs[0, 1], column='elevation', cmap='viridis',
                       legend=True, markersize=20)
    axs[0, 1].set_title('Aligned LiDAR Points')

    # Plot adjusted points
    adjusted_points.plot(ax=axs[1, 0], column='elevation', cmap='viridis',
                        legend=True, markersize=20)
    axs[1, 0].set_title('Adjusted LiDAR Points')

    # Plot corrected points
    corrected_points.plot(ax=axs[1, 1], column='elevation', cmap='viridis',
                         legend=True, markersize=20)
    axs[1, 1].set_title('Corrected LiDAR Points')

    plt.tight_layout()
    plt.savefig('salpy_examples.png')
    print("\nVisualization saved as 'salpy_examples.png'")

if __name__ == "__main__":
    main()