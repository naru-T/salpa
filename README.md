# salpa: Satellite LiDAR Point Alignment

This R package provides tools for aligning satellite LiDAR data with geographic reference data through optimized positional correction.

## Features

- Multiple optimization algorithms:
  - Genetic Algorithm (GA)
  - Particle Swarm Optimization (PSO)
  - Whale Optimization Algorithm (WOA)
  - Limited-memory BFGS with Box constraints (L-BFGS-B)
- Various error metrics:
  - Euclidean distance
  - Manhattan distance
  - Dynamic Time Warping (DTW)
  - Correlation-based
  - Area-based
  - Hausdorff distance
- Efficient batch processing
- Support for parallel computation

## Installation

```r
# Install from GitHub
devtools::install_github("naru-T/salpa")
```

## Basic Usage

```r
library(salpa)
library(sf)
library(terra)

# Load lidar data
lidar_points <- st_read("lidar_points.gpkg")

# Load reference DEM
dem <- rast("reference_dem.tif")

# Apply positional correction
result <- positional_correction(
  lidar_footprints = lidar_points,
  input_rast = dem,
  minimizing_method = "euclidean",
  target_variable = "mean",
  lidar_value = "elevation",
  optimization_method = "ga" 
)

# Access the corrected footprints
corrected_points <- result$corrected_footprints

# Get optimization results
optim_results <- result$optim_result
```

## Optimization Methods

The package implements four different optimization methods, each with distinct characteristics:

### Genetic Algorithm (GA)
- Mimics natural selection processes
- Good for complex search spaces with many local optima
- Parameter `optimization_method = "ga"`
- Additional parameters through `pop_size` and `max_iter`

```r
# Example with GA-specific parameters
ga_result <- positional_correction(
  lidar_footprints = lidar_points,
  input_rast = dem,
  optimization_method = "ga",
  pop_size = 50,      # Population size
  max_iter = 100,     # Maximum iterations
  parallel = TRUE     # Enable parallel processing
)
```

### Particle Swarm Optimization (PSO)
- Based on social behavior of bird flocking or fish schooling
- Efficient for continuous optimization problems
- Generally faster convergence than GA for many problems
- Parameter `optimization_method = "pso"`
- Additional parameters can be passed through `pso_params` list

```r
# Example with PSO-specific parameters
pso_result <- positional_correction(
  lidar_footprints = lidar_points,
  input_rast = dem,
  optimization_method = "pso",
  pop_size = 30,      # Swarm size
  max_iter = 100,     # Maximum iterations
  parallel = TRUE,    # Enable parallel processing
  pso_params = list(
    trace = FALSE,    # Whether to print progress
    abstol = 1e-4,    # Absolute tolerance
    reltol = 1e-4,    # Relative tolerance
    REPORT = 10       # Report every 10 iterations
  )
)
```

### Whale Optimization Algorithm (WOA)
- Mimics hunting behavior of humpback whales
- Balances exploration and exploitation phases
- Effective for multimodal optimization problems
- Parameter `optimization_method = "woa"`
- Additional parameters through `woa_params` list

```r
# Example with WOA-specific parameters
woa_result <- positional_correction(
  lidar_footprints = lidar_points,
  input_rast = dem,
  optimization_method = "woa",
  pop_size = 30,      # Number of whales
  max_iter = 100,     # Maximum iterations
  woa_params = list(
    verbose = TRUE,   # Whether to print progress
    batch_size = 10   # Process in batches (more memory efficient)
  )
)
```

### Limited-memory BFGS with Box constraints (L-BFGS-B)
- Gradient-based optimization method with memory efficiency
- Typically fastest for well-behaved, smooth functions
- Most efficient for small search spaces
- Parameter `optimization_method = "lbfgsb"`
- Additional parameters through `lbfgsb_params` list

```r
# Example with L-BFGS-B specific parameters
lbfgsb_result <- positional_correction(
  lidar_footprints = lidar_points,
  input_rast = dem,
  optimization_method = "lbfgsb",
  lbfgsb_params = list(
    factr = 1e7,     # Convergence factor
    pgtol = 1e-5,    # Projected gradient tolerance
    trace = 1        # Tracing level (0-6)
  )
)
```

## Performance Considerations

Each optimization method has strengths for different scenarios:

- **L-BFGS-B**: Fastest convergence for well-behaved functions with smooth gradients. For small search spaces (≤60 meters), L-BFGS-B can be 2-5x faster than other methods.

- **GA**: Most robust for complex, non-convex problems with many local optima. Requires more function evaluations but parallel processing can significantly improve performance.

- **PSO**: Good balance between exploration speed and convergence accuracy. Performance can be significantly improved by:
  - Using smaller swarm sizes (10-30) for faster convergence
  - Setting appropriate tolerance parameters to avoid unnecessary iterations
  - Using parallel processing for large datasets
  - Employing early stopping criteria for well-behaved functions

- **WOA**: Excellent at escaping local minima with its spiral search pattern. The built-in batch processing capability makes it memory efficient for large datasets.

### Method Selection Guidelines

- For small adjustments (a few dozen meters) or relatively flat terrain: **L-BFGS-B**
- For complex terrain with many local optima: **GA** or **WOA**
- For best balance of speed and accuracy: **PSO** with tuned parameters
- For large datasets that don't fit in memory: **WOA** with batch processing

## Distance Metrics

The package implements multiple distance metrics for comparing elevation profiles:

- **Euclidean Distance** (`minimizing_method = "euclidean"`): Standard point-to-point distance
- **Dynamic Time Warping** (`minimizing_method = "dtw"`): Handles temporal shifts in sequences
- **Manhattan Distance** (`minimizing_method = "manhattan"`): Sum of absolute differences
- **Correlation Distance** (`minimizing_method = "correlation"`): Based on correlation coefficients
- **Area Distance** (`minimizing_method = "area"`): Area between profiles
- **Hausdorff Distance** (`minimizing_method = "hausdorff"`): Maximum of all minimum distances

## Advanced Usage

### Parallel Processing

Enable parallel processing to speed up computation:

```r
result <- positional_correction(
  lidar_footprints = lidar_points,
  input_rast = dem,
  minimizing_method = "euclidean",
  target_variable = "mean",
  lidar_value = "elevation",
  optimization_method = "pso",
  parallel = TRUE  # Enable parallel processing
)
```

### Custom Search Bounds

Specify custom search bounds for the optimization:

```r
result <- positional_correction(
  lidar_footprints = lidar_points,
  input_rast = dem,
  minimizing_method = "euclidean",
  target_variable = "mean",
  lidar_value = "elevation",
  optimization_method = "pso",
  lower_bounds = c(-50, -50),  # Lower bounds for x and y
  upper_bounds = c(50, 50)     # Upper bounds for x and y
)

```

## Inquiry
For questions and feedback:
- GitHub issues: [https://github.com/naru-T/salpa/issues](https://github.com/naru-T/salpa/issues)

## License

MIT
