# salpa: Satellite LiDAR Point Alignment (v0.0.2.2)

R package for aligning satellite LiDAR footprints (e.g. GEDI) to a high-resolution DEM by optimizing a 2D positional offset.

Requires **R ≥ 4.1** (uses the native pipe `|>`). Population-based optimizers (GA / PSO / WOA) use a **custom C++ implementation** with fork-safe parallel fitness evaluation (`parallel::mclapply` + `terra::wrap`). There is **no runtime dependency** on the CRAN `GA` or `pso` packages.

## Features

- Optimization algorithms:
  - Genetic Algorithm (`"ga"`)
  - Particle Swarm Optimization (`"pso"`)
  - Whale Optimization Algorithm (`"woa"`)
  - Limited-memory BFGS with box constraints (`"lbfgsb"`)
- Distance / loss metrics (`minimizing_method`):
  - `"euclidean"`, `"manhattan"`, `"dtw"`, `"correlation"`, `"area"`, `"hausdorff"`
- Optional `geographic_coordinates = TRUE` to treat search bounds as metres and convert them to degrees for WGS84-style DEMs
- Parallel fitness evaluation via `parallel = TRUE` or a numeric worker count (e.g. `parallel = 8`)

## Installation

```r
# Recommended: install the released source tarball
install.packages(
  "https://github.com/naru-T/salpa/releases/download/0.0.2.2/salpa_0.0.2.2.tar.gz",
  repos = NULL,
  type = "source"
)

# Or from GitHub (needs compilation; LinkingTo Rcpp)
# remotes::install_github("naru-T/salpa")
```

## Basic usage

`lidar_value` is **required** (column of LiDAR elevations). DEM and footprints must share the same CRS as `crs_code`.

```r
library(salpa)
library(sf)
library(terra)

lidar_points <- st_read("lidar_points.gpkg")
dem <- rast("reference_dem.tif")

result <- positional_correction(
  lidar_footprints = lidar_points,
  input_rast = dem,
  minimizing_method = "euclidean",
  target_variable = "mean",
  lidar_value = "elevation",
  buf = 12.5,
  crs_code = 6677,                 # must match DEM / footprints
  optimization_method = "ga",
  lower_bounds = c(-25, -25),      # metres (default)
  upper_bounds = c(25, 25),
  pop_size = 50,
  max_iter = 100,
  parallel = 8
)

corrected_points <- result$corrected_footprints
optim_results <- result$optim_result   # best_x, best_y, best_value
```

### Geographic (lon/lat) DEMs

If the DEM is in geographic coordinates and you want bounds in **metres**:

```r
result <- positional_correction(
  lidar_footprints = lidar_points,
  input_rast = dem,
  lidar_value = "elevation",
  crs_code = 4326,
  geographic_coordinates = TRUE,   # convert ±25 m bounds to degrees
  lower_bounds = c(-25, -25),
  upper_bounds = c(25, 25),
  optimization_method = "ga",
  parallel = TRUE
)
```

## Optimization methods

Defaults for `positional_correction()`: `pop_size = 25`, `max_iter = 50`,
`lower_bounds = c(-25, -25)`, `upper_bounds = c(25, 25)`, `parallel = FALSE`.

### Genetic Algorithm (`optimization_method = "ga"`)

Custom GA (BLX-α crossover, tournament selection) with C++ population updates.

Optional knobs via `ga_params`:

| Key | Meaning |
|---|---|
| `pcrossover` | Crossover probability (default `0.8`) |
| `pmutation` | Mutation probability (default `0.1`) |
| `run` | Early-stop patience in generations |
| `elitism` | Number of elites kept |
| `tournament_k` | Tournament size |

```r
ga_result <- positional_correction(
  lidar_footprints = lidar_points,
  input_rast = dem,
  lidar_value = "elevation",
  crs_code = 6677,
  optimization_method = "ga",
  pop_size = 50,
  max_iter = 100,
  parallel = TRUE,
  ga_params = list(pcrossover = 0.8, pmutation = 0.1, run = 25)
)
```

### Particle Swarm Optimization (`optimization_method = "pso"`)

Custom PSO (not the CRAN `pso` package).

Optional knobs via `pso_params`:

| Key | Meaning |
|---|---|
| `w` | Inertia weight |
| `c1` | Cognitive coefficient |
| `c2` | Social coefficient |
| `v_max` | Max velocity |
| `abstol` | Improvement tolerance for stagnation (default `1e-3`) |
| `stagnation_limit` | Iterations without improvement before stop |
| `trace` | `TRUE` to print progress |

```r
pso_result <- positional_correction(
  lidar_footprints = lidar_points,
  input_rast = dem,
  lidar_value = "elevation",
  crs_code = 6677,
  optimization_method = "pso",
  pop_size = 30,
  max_iter = 100,
  parallel = TRUE,
  pso_params = list(
    w = 0.7,
    c1 = 1.5,
    c2 = 1.5,
    abstol = 1e-4,
    trace = FALSE
  )
)
```

### Whale Optimization Algorithm (`optimization_method = "woa"`)

Custom WOA with C++ position updates.

Optional knobs via `woa_params`:

| Key | Meaning |
|---|---|
| `verbose` | Print progress |
| `early_stop` | Enable early stopping (default `TRUE`) |
| `convergence_threshold` | Loss improvement threshold |
| `stagnation_limit` | Iterations without improvement |
| `memoization` | Cache objective evaluations (default `TRUE`) |

```r
woa_result <- positional_correction(
  lidar_footprints = lidar_points,
  input_rast = dem,
  lidar_value = "elevation",
  crs_code = 6677,
  optimization_method = "woa",
  pop_size = 30,
  max_iter = 100,
  parallel = TRUE,
  woa_params = list(verbose = TRUE, early_stop = TRUE)
)
```

### L-BFGS-B (`optimization_method = "lbfgsb"`)

Gradient-based local search (`stats::optim` method `"L-BFGS-B"`). Fast on smooth losses; can settle in local minima on flat multimodal surfaces (e.g. area).

Optional knobs via `lbfgsb_params` (passed into the optimizer control), including:

| Key | Meaning |
|---|---|
| `multi_start` | `TRUE` to run several starts (default `FALSE`) |
| `n_starts` | Number of starts when `multi_start = TRUE` |
| `factr` | Convergence factor for `optim` |
| `pgtol` | Projected gradient tolerance |
| `trace` | `optim` tracing level |

```r
lbfgsb_result <- positional_correction(
  lidar_footprints = lidar_points,
  input_rast = dem,
  lidar_value = "elevation",
  crs_code = 6677,
  optimization_method = "lbfgsb",
  max_iter = 100,
  lbfgsb_params = list(
    factr = 1e7,
    pgtol = 1e-5,
    multi_start = FALSE
  )
)
```

### Method selection guidelines

- Small shifts / relatively smooth loss: **L-BFGS-B**
- Complex terrain / many local optima: **GA** or **WOA**
- Speed / accuracy balance: **PSO**
- Parallelism helps GA / PSO / WOA (`parallel = TRUE` or a core count); L-BFGS-B is typically single-start

## Distance metrics

| `minimizing_method` | Role |
|---|---|
| `"euclidean"` | Mean squared / Euclidean profile difference |
| `"manhattan"` | Sum of absolute differences |
| `"dtw"` | Dynamic time warping |
| `"correlation"` | Correlation-based distance |
| `"area"` | Area between elevation profiles |
| `"hausdorff"` | Hausdorff distance between profiles |

## Custom search bounds

```r
result <- positional_correction(
  lidar_footprints = lidar_points,
  input_rast = dem,
  minimizing_method = "euclidean",
  target_variable = "mean",
  lidar_value = "elevation",
  crs_code = 6677,
  optimization_method = "pso",
  lower_bounds = c(-50, -50),
  upper_bounds = c(50, 50),
  parallel = 8
)
```

## Parallel evaluation

- `parallel = FALSE` — serial (default)
- `parallel = TRUE` — use `detectCores() - 1` workers (Unix fork / `mclapply`)
- `parallel = 8` — use 8 workers

On Windows, evaluation falls back to serial (no fork).

## Inquiry

Questions and feedback: [GitHub issues](https://github.com/naru-T/salpa/issues)

## License

GPL-3. See [LICENSE.md](LICENSE.md).
