# salpa: Satellite LiDAR Point Adjustment <img src="inst/hex/salpa_hex.png" align="right" height="150" />

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/929837682.svg)](https://doi.org/10.5281/zenodo.15041709)
<!-- badges: end -->

## Overview

`salpa` is an R package designed for optimizing and adjusting satellite LiDAR point positions. The package provides tools for improving the spatial accuracy of satellite-based LiDAR observations through two main adjustment methods:

1. **Linear Alignment**: Aligns points to their best-fit line
2. **Positional Correction**: Fine-tunes point locations with specified offsets to minimize elevation differences using optimization algorithms

## Installation

You can install the development version of salpa from GitHub:

```r
# install.packages("devtools")
devtools::install_github("naru-T/salpa")
```

## Features

- **Linear Alignment**: Corrects systematic linear offsets in LiDAR data
- **Positional Correction**: Optimizes point positions using advanced optimization methods
- **Multiple Optimization Methods**: Supports both Genetic Algorithm (GA) and Particle Swarm Optimization (PSO)
- **Multiple Distance Metrics**: Supports various distance calculations (Euclidean, DTW, Manhattan, etc.)
- **Flexible Buffer Options**: Customizable buffer sizes for point adjustments
- **Parallel Processing**: Support for parallel computation in optimization

## Dependencies

### Required Packages
- sf: For spatial data handling
- GA: For genetic algorithm optimization
- pso: For particle swarm optimization
- exactextractr: For raster value extraction
- stats: For statistical calculations
- dplyr: For data manipulation

### Suggested Packages
- terra: For raster data handling
- tmap & tmaptools: For mapping
- elevatr: For elevation data access
- ggplot2: For visualization
- tidyr: For data reshaping
- testthat: For unit testing

## Basic Usage

```r
library(salpa)
library(sf)
library(terra)

# Load sample data
gpkg_path <- system.file("extdata", "gedi_l2a_shot_sample.gpkg", package = "salpa")
lidar_footprints <- st_read(gpkg_path)

# Step 1: Linear Alignment
aligned_points <- linear_alignment(lidar_footprints, crs_code = 3857)

# Step 2: Positional Correction with Genetic Algorithm (default)
ga_corrected_positions <- positional_correction(
  lidar_footprints = aligned_points,
  input_rast = dem_rast,           # Your reference DEM
  minimizing_method = "euclidean", # Distance metric for comparison
  target_variable = "mean",
  lidar_value = "elev_lowestmode", # Column containing LiDAR elevation values
  buf = 12.5,                      # Buffer size in meters
  crs_code = 3857,
  optimization_method = "ga",      # Use Genetic Algorithm (default)
  parallel = TRUE
)

# Step 2 (Alternative): Positional Correction with Particle Swarm Optimization
pso_corrected_positions <- positional_correction(
  lidar_footprints = aligned_points,
  input_rast = dem_rast,
  minimizing_method = "euclidean",
  target_variable = "mean",
  lidar_value = "elev_lowestmode",
  buf = 12.5,
  crs_code = 3857,
  optimization_method = "pso",     # Use Particle Swarm Optimization (faster alternative)
  parallel = TRUE
)
```

## Optimization Methods

The package supports multiple optimization methods:

1. **Genetic Algorithm (GA)**: The original implementation, effective for complex optimization problems with noisy spatial data
2. **Particle Swarm Optimization (PSO)**: An alternative that may offer performance advantages for certain datasets

Performance varies by dataset characteristics, with each method having strengths for different problem types:
- **GA**: Often performs better on discrete problems or when the search space has many local optima
- **PSO**: May perform better on continuous problems with smoother gradients

Both methods support parallel processing for improved performance and achieve similar accuracy levels. We recommend benchmarking both methods on your specific data to determine the optimal approach.

## Distance Metrics

The package supports multiple distance metrics for comparing sequences:

1. **Euclidean** (default): Standard point-to-point distance
2. **DTW**: Dynamic Time Warping for sequence alignment
3. **Manhattan**: Sum of absolute differences
4. **Correlation**: Based on Pearson correlation
5. **Area**: Area between curves
6. **Hausdorff**: Maximum minimum distance

Each metric handles NA values appropriately and can be selected using the `minimizing_method` parameter.

## Documentation

For more detailed information and examples, please refer to:
- Package vignette: `vignette("salpa")`
- Function documentation: `?linear_alignment`, `?positional_correction`, `?perform_distance`

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## Author

Narumasa Tsutsumida (ORCID: 0000-0002-6333-0301)

## License

This project is licensed under the GPL-3 License - see the LICENSE file for details.

## Citation

If you use this package in your research, please cite it as:

```r
Tsutsumida N. (2025) salpa: Satellite LiDAR Point Adjustment, R package version 0.0.1.3, https://github.com/naru-T/salpa
```

## Contact

For questions and feedback:
- GitHub issues: [https://github.com/naru-T/salpa/issues](https://github.com/naru-T/salpa/issues)

## Notes

The photo used in the hex sticker is from [wikipedia](https://en.wikipedia.org/wiki/Salp#/media/File:Sea_Salp_Chain.jpg).
