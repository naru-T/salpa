# salpa: Satellite LiDAR Point Adjustment <img src="inst/hex/salpa_hex.png" align="right" height="150" />

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Overview

`salpa` is an R package designed for optimizing and adjusting satellite LiDAR point positions using genetic algorithms. The package provides tools for improving the spatial accuracy of satellite-based LiDAR observations through two main adjustment methods:

1. **Linear Alignment**: Aligns points to their best-fit line using genetic algorithm optimization
2. **Positional Correction**: Fine-tunes point locations with specified offsets to minimize elevation differences

## Installation

You can install the development version of salpa from GitHub:

```r
# install.packages("devtools")
devtools::install_github("naru-T/salpa")
```

## Features

- **Linear Alignment**: Corrects systematic linear offsets in LiDAR data
- **Positional Correction**: Optimizes point positions using genetic algorithms
- **Flexible Buffer Options**: Customizable buffer sizes for point adjustments
- **Parallel Processing**: Support for parallel computation in optimization

## Dependencies

### Required Packages
- sf: For spatial data handling
- GA: For genetic algorithm optimization
- exactextractr: For raster value extraction
- dtw: For dynamic time warping calculations
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

# Step 2: Positional Correction
corrected_positions <- positional_correction(
  lidar_footprints = aligned_points,
  input_rast = dem_rast,           # Your reference DEM
  minimizing_method = "dtw",
  target_variable = "mean",
  lidar_value = "elev_lowestmode", # Column containing LiDAR elevation values
  buf = 12.5,                      # Buffer size in meters
  crs_code = 3857,
  parallel = TRUE
)
```

## Documentation

For more detailed information and examples, please refer to:
- Package vignette: `vignette("salpa")`
- Function documentation: `?linear_alignment`, `?positional_correction`

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## Author

Narumasa Tsutsumida (ORCID: 0000-0002-6333-0301)

## License

This project is licensed under the GPL-3 License - see the LICENSE file for details.

## Citation

If you use this package in your research, please cite it as:

```r
citation("salpa")
```

## Contact

For questions and feedback:
- GitHub issues: [https://github.com/naru-T/salpa/issues](https://github.com/naru-T/salpa/issues)
- Email: rsnaru.jp@gmail.com 

## Notes

The photo used in the hex sticker is from [wikipedia](https://en.wikipedia.org/wiki/Salp#/media/File:Sea_Salp_Chain.jpg).