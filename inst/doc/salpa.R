## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# # install.packages("remotes")
# remotes::install_github("naru-T/salpa")

## -----------------------------------------------------------------------------
library(salpa)
library(sf)
library(terra)
library(exactextractr)
library(raster)
library(ggplot2)
library(dplyr)

## -----------------------------------------------------------------------------
# Load GEDI L2A sample data
gpkg_path <- system.file("extdata", "gedi_l2a_shot_sample.gpkg", package = "salpa")
lidar_footprints <- st_read(gpkg_path, quiet = TRUE)

lidar_footprints <- lidar_footprints %>% mutate(gedi_dem = as.numeric(elev_lowestmode) - as.numeric(geoid))
# Look at the data
#head(lidar_footprints)
plot(lidar_footprints["gedi_dem"], pch = 20, main = "GEDI L2A Lidar Footprints")

## -----------------------------------------------------------------------------
# Define area of interest
bbox_4326 <- st_bbox(c(xmin = 139.40119, ymin = 36.86163,
                      xmax = 139.52543, ymax = 36.94841),
                    crs = 4326)
bbox_sf <- st_as_sf(st_as_sfc(bbox_4326))


# Get reference DEM
dem <- elevatr::get_elev_raster(bbox_sf, z = 12, clip = "bbox")
dem_3857 <- raster::projectRaster(dem, crs = 3857)
dem_rast <- terra::rast(dem_3857)

# Transform lidar footprints to Web Mercator
lidar_footprints_3857 <- st_transform(lidar_footprints, 3857)

# Plot DEM
plot(dem_rast, main = "Reference DEM")
plot(st_geometry(lidar_footprints_3857), add = TRUE, pch = 20, col = "red")

## -----------------------------------------------------------------------------
# Apply positional correction using GA
ga_result <- positional_correction(
  lidar_footprints = lidar_footprints_3857,
  input_rast = dem_rast,
  minimizing_method = "euclidean",
  target_variable = "mean",
  lidar_value = "gedi_dem",
  buf = 12.5,
  lower_bounds = c(-1000, -1000),
  upper_bounds = c(1000, 1000),
  crs_code = 3857,
  optimization_method = "ga",
  max_iter = 200,  
  pop_size = 20
)

# Print optimization result
print(ga_result$optim_result)

## -----------------------------------------------------------------------------
# Apply positional correction using PSO
pso_result <- positional_correction(
  lidar_footprints = lidar_footprints_3857,
  input_rast = dem_rast,
  minimizing_method = "euclidean",
  target_variable = "mean",
  lidar_value = "gedi_dem",
  buf = 12.5,
  lower_bounds = c(-1000, -1000),
  upper_bounds = c(1000, 1000),
  crs_code = 3857,
  optimization_method = "pso",
  max_iter = 200,  
  pop_size = 20,
  pso_params = list(
    trace = TRUE,     # Whether to print progress
    abstol = 1e-4,     # Absolute tolerance
    reltol = 1e-4,     # Relative tolerance
    REPORT = 10         # Report every 5 iterations
  )
)

# Print optimization result
print(pso_result$optim_result)

## -----------------------------------------------------------------------------
# Apply positional correction using WOA
woa_result <-  positional_correction(
    lidar_footprints = lidar_footprints_3857,
    input_rast = dem_rast,
    minimizing_method = "euclidean",
    target_variable = "mean",
    lidar_value = "gedi_dem",
    buf = 12.5,
    lower_bounds = c(-1000, -1000),
    upper_bounds = c(1000, 1000),
    crs_code = 3857,
    optimization_method = "woa",   # Whale Optimization Algorithm
    max_iter = 200,
    pop_size = 20,
    woa_params = list(
      verbose = TRUE,   # Whether to print progress
      batch_size = 10    # Process all positions in one batch
    )
  )

# Print optimization result if available
print(woa_result$optim_result)

## -----------------------------------------------------------------------------
# Apply positional correction using L-BFGS-B with error handling
lbfgsb_result <-  positional_correction(
    lidar_footprints = lidar_footprints_3857,
    input_rast = dem_rast,
    minimizing_method = "euclidean",
    target_variable = "mean",
    lidar_value = "gedi_dem",
    buf = 12.5,
    lower_bounds = c(-1000, -1000),
    upper_bounds = c(1000, 1000),
    crs_code = 3857,
    optimization_method = "lbfgsb",  # L-BFGS-B
    max_iter = 200,
    lbfgsb_params = list(
      factr = 1e7,        # Control precision
      pgtol = 1e-5,       # Gradient tolerance
      trace = 0           # No tracing
    )
  )

# Print optimization result if available
print(lbfgsb_result$optim_result)


## -----------------------------------------------------------------------------
# Function to safely get values from optimization results
posit_result_get <- function(result, default = 0) {
  if (!is.null(result) && !is.null(result$optim_result)) {
    return(result$optim_result)
  } else {
    # Return default values if result is NULL
    warning("Optimization result was NULL, using default values instead")
    return(list(best_x = default, best_y = default, best_value = NA))
  }
}

# Create a data frame of adjustment values
adjustments <- data.frame(
  Method = c("GA", "PSO", "WOA", "L-BFGS-B"),
  X_Adjustment = c(
    posit_result_get(ga_result)$best_x,
    posit_result_get(pso_result)$best_x,
    posit_result_get(woa_result)$best_x,
    posit_result_get(lbfgsb_result)$best_x
  ),
  Y_Adjustment = c(
    posit_result_get(ga_result)$best_y,
    posit_result_get(pso_result)$best_y,
    posit_result_get(woa_result)$best_y,
    posit_result_get(lbfgsb_result)$best_y
  ),
  Fitness_Value = c(
    posit_result_get(ga_result)$best_value,
    posit_result_get(pso_result)$best_value,
    posit_result_get(woa_result)$best_value,
    posit_result_get(lbfgsb_result)$best_value
  )
)

# Print the adjustment values
print(adjustments)

# Plot X and Y adjustments
ggplot(adjustments, aes(x = X_Adjustment, y = Y_Adjustment, color = Method)) +
  geom_point(size = 4) +
  geom_text(aes(label = Method), vjust = -1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  labs(title = "Comparison of X/Y Positional Adjustments",
       x = "X Adjustment (meters)",
       y = "Y Adjustment (meters)")

# Plot fitness values (excluding NA values)
adjustments_filtered <- adjustments[!is.na(adjustments$Fitness_Value),]
if(nrow(adjustments_filtered) > 0) {
  ggplot(adjustments_filtered, aes(x = Method, y = Fitness_Value, fill = Method)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Comparison of Fitness Values",
         y = "Fitness Value (lower is better)",
         x = "Method")
} else {
  cat("No valid fitness values to plot")
}

## -----------------------------------------------------------------------------
# Create buffered footprints for extraction
buffered_original <- st_buffer(lidar_footprints_3857, 12.5)

# Use safe buffering function to handle possible NULL results
make_buffer <- function(result, crs_code = 3857, buffer_size = 12.5) {
  if (!is.null(result)) {
    best_x <- result$optim_result$best_x
    best_y <- result$optim_result$best_y
    adjusted_footprints <- position_adjustment(lidar_footprints_3857, best_x, best_y, crs_code)
    return(st_buffer(adjusted_footprints, buffer_size))
  } else {
    # Return the original footprints if result is NULL
    warning("Optimization result was NULL, using original footprints instead")
    return(buffered_original)
  }
}

# Safely create buffered footprints
buffered_ga <- make_buffer(ga_result)
buffered_pso <- make_buffer(pso_result)
buffered_woa <- make_buffer(woa_result)
buffered_lbfgsb <- make_buffer(lbfgsb_result)

# Extract values from each buffered dataset
original_extracted <- exact_extract(dem_rast, buffered_original, "mean", progress = FALSE)
ga_extracted <- exact_extract(dem_rast, buffered_ga, "mean", progress = FALSE)
pso_extracted <- exact_extract(dem_rast, buffered_pso, "mean", progress = FALSE)
woa_extracted <- exact_extract(dem_rast, buffered_woa, "mean", progress = FALSE)
lbfgsb_extracted <- exact_extract(dem_rast, buffered_lbfgsb, "mean", progress = FALSE)

# Get lidar values
lidar_values <- as.numeric(lidar_footprints_3857$gedi_dem)

# Print summary stats for each dataset
cat("Original values summary:\n")
print(summary(original_extracted))
cat("\nGA corrected values summary:\n")
print(summary(ga_extracted))
cat("\nPSO corrected values summary:\n")
print(summary(pso_extracted))
cat("\nWOA corrected values summary:\n")
print(summary(woa_extracted))
cat("\nL-BFGS-B corrected values summary:\n")
print(summary(lbfgsb_extracted))
cat("\nLidar values summary:\n")
print(summary(lidar_values))

## -----------------------------------------------------------------------------
# Create a safe function to calculate MAE
calc_mae <- function(lidar_values, extracted_values, method_name) {
  if (length(extracted_values) == 0 || all(is.na(extracted_values))) {
    warning("No valid values for method: ", method_name)
    return(NA)
  }
  
  # Ensure both vectors are the same length
  min_length <- min(length(lidar_values), length(extracted_values))
  
  if (min_length == 0) {
    warning("Zero-length vectors for method: ", method_name)
    return(NA)
  }
  
  # Calculate MAE with proper handling of NA values
  mae <- mean(abs(lidar_values[1:min_length] - extracted_values[1:min_length]), na.rm = TRUE)
  
  if (is.nan(mae) || is.na(mae)) {
    warning("MAE calculation resulted in NaN for method: ", method_name)
    return(NA)
  }
  
  return(mae)
}

# Calculate MAE for each method
original_mae <- calc_mae(lidar_values, original_extracted, "Original")
ga_mae <- calc_mae(lidar_values, ga_extracted, "GA")
pso_mae <- calc_mae(lidar_values, pso_extracted, "PSO")
woa_mae <- calc_mae(lidar_values, woa_extracted, "WOA")
lbfgsb_mae <- calc_mae(lidar_values, lbfgsb_extracted, "L-BFGS-B")


# Create a data frame for plotting
mae_df <- data.frame(
  Method = c("Original", "GA", "PSO", "WOA", "L-BFGS-B"),
  MAE = c(original_mae, ga_mae, pso_mae, woa_mae, lbfgsb_mae)
)

# Print MAE values
print(mae_df)

# Filter out NA values for plotting
mae_df_filtered <- mae_df[!is.na(mae_df$MAE),]

# Plot MAE comparison if there are valid values
if (nrow(mae_df_filtered) > 0) {
  ggplot(mae_df_filtered, aes(x = Method, y = MAE, fill = Method)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Mean Absolute Error Comparison",
         y = "Mean Absolute Error (m)",
         x = "Method")
} else {
  cat("No valid MAE values to plot")
}

## -----------------------------------------------------------------------------
# Function to run optimization method with timing
time_optimization <- function(method, params) {
  start_time <- Sys.time()
  
  result <- tryCatch({
    # Create a new list with the optimization_method parameter
    all_params <- params
    all_params$optimization_method <- method
    
    # Call the function with do.call to properly pass all parameters
    do.call(positional_correction, all_params)
  }, error = function(e) {
    warning(paste(method, "optimization failed:", conditionMessage(e)))
    return(NULL)
  })
  
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  return(list(
    result = result,
    time = execution_time
  ))
}

# Common parameters for all methods
common_params <- list(
  lidar_footprints = lidar_footprints_3857,
  input_rast = dem_rast,
  minimizing_method = "euclidean",
  target_variable = "mean",
  lidar_value = "gedi_dem",
  buf = 12.5,
  lower_bounds = c(-1000, -1000),
  upper_bounds = c(1000, 1000),
  crs_code = 3857,
  max_iter = 100,
  pop_size = 20
)

# Time each method
cat("Running timing tests (this may take a while)...\n")

ga_timing <- time_optimization("ga", common_params)
pso_timing <- time_optimization("pso", common_params)
woa_timing <- time_optimization("woa", common_params)
lbfgsb_timing <-  time_optimization("lbfgsb", common_params)

# Create timing comparison data frame
timing_df <- data.frame(
  Method = c("GA", "PSO", "WOA", "L-BFGS-B"),
  Time_Seconds = c(
    ga_timing$time,
    pso_timing$time,
    woa_timing$time,
    lbfgsb_timing$time
  )
)

# Print timing results
print(timing_df)

# Plot timing comparison (excluding NA values)
timing_df_filtered <- timing_df[!is.na(timing_df$Time_Seconds),]

if(nrow(timing_df_filtered) > 0) {
  ggplot(timing_df_filtered, aes(x = Method, y = Time_Seconds, fill = Method)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Execution Time Comparison",
         y = "Time (seconds)",
         x = "Method") +
    geom_text(aes(label = sprintf("%.1f", Time_Seconds)), vjust = -0.5)
} else {
  cat("No valid timing values to plot")
}

## -----------------------------------------------------------------------------
# Function to create trajectory lines between original and corrected positions
create_trajectories <- function(original_sf, corrected_sf, method_name) {
  # Extract centroids
  orig_points <- st_centroid(original_sf)
  corr_points <- st_centroid(corrected_sf)
  
  # Create a data frame with line start and end coordinates
  n_points <- nrow(orig_points)
  
  # Take a sample if there are too many points
  max_sample <- 50
  if (n_points > max_sample) {
    sample_indices <- sample(1:n_points, max_sample)
    orig_points <- orig_points[sample_indices,]
    corr_points <- corr_points[sample_indices,]
    n_points <- max_sample
  }
  
  orig_coords <- st_coordinates(orig_points)
  corr_coords <- st_coordinates(corr_points)
  
  # Create line segments for each point pair
  lines_list <- list()
  for (i in 1:n_points) {
    pt1 <- c(orig_coords[i, 1], orig_coords[i, 2])
    pt2 <- c(corr_coords[i, 1], corr_coords[i, 2])
    lines_list[[i]] <- st_linestring(rbind(pt1, pt2))
  }
  
  # Convert to sf object
  lines_sf <- st_sf(
    ID = 1:length(lines_list),
    Method = rep(method_name, length(lines_list)),
    geometry = st_sfc(lines_list, crs = st_crs(original_sf))
  )
  
  return(lines_sf)
}

# Access the position_adjustment correctly based on package structure
get_adjusted_footprints <- function(result) {
  if (is.null(result)) return(NULL)
  
  # Check which slot contains the adjusted data
  if (!is.null(result$position_adjustment)) {
    return(result$position_adjustment)
  } else if (!is.null(result$corrected_footprints)) {
    return(result$corrected_footprints)
  } else {
    return(NULL)
  }
}

# Create trajectories if results are available
trajectories_list <- list()
trajectory_count <- 1

# Add GA trajectories
ga_adjusted <- get_adjusted_footprints(ga_result)
if (!is.null(ga_adjusted)) {
  ga_trajectories <- create_trajectories(lidar_footprints_3857, ga_adjusted, "GA")
  trajectories_list[[trajectory_count]] <- ga_trajectories
  trajectory_count <- trajectory_count + 1
}

# Add PSO trajectories
pso_adjusted <- get_adjusted_footprints(pso_result)
if (!is.null(pso_adjusted)) {
  pso_trajectories <- create_trajectories(lidar_footprints_3857, pso_adjusted, "PSO")
  trajectories_list[[trajectory_count]] <- pso_trajectories
  trajectory_count <- trajectory_count + 1
}

# Add WOA trajectories if available
woa_adjusted <- get_adjusted_footprints(woa_result)
if (!is.null(woa_adjusted)) {
  woa_trajectories <- create_trajectories(lidar_footprints_3857, woa_adjusted, "WOA")
  trajectories_list[[trajectory_count]] <- woa_trajectories
  trajectory_count <- trajectory_count + 1
}

# Add L-BFGS-B trajectories if available
lbfgsb_adjusted <- get_adjusted_footprints(lbfgsb_result)
if (!is.null(lbfgsb_adjusted)) {
  lbfgsb_trajectories <- create_trajectories(lidar_footprints_3857, lbfgsb_adjusted, "L-BFGS-B")
  trajectories_list[[trajectory_count]] <- lbfgsb_trajectories
  trajectory_count <- trajectory_count + 1
}

# Print number of available methods for plotting
cat("Number of methods with available trajectory data:", length(trajectories_list), "\n")

# Check available methods for plotting
if (length(trajectories_list) > 0) {
  all_trajectories <- do.call(rbind, trajectories_list)
  
  # Create a plot of the trajectories
  plot1 <- ggplot() +
    geom_sf(data = all_trajectories, aes(color = Method), size = 0.5) +
    geom_sf(data = st_centroid(lidar_footprints_3857), size = 1, color = "black") +
    theme_minimal() +
    labs(title = "Footprint Displacement Trajectories",
         subtitle = "Black points: Original positions, Colored lines: Displacement vectors") +
    scale_color_brewer(palette = "Set1") +
    coord_sf(datum = NA)
  
  print(plot1)
  
  # Create a faceted plot by method
  if (length(trajectories_list) > 1) {
    plot2 <- ggplot() +
      geom_sf(data = all_trajectories, aes(color = Method), size = 0.5) +
      geom_sf(data = st_centroid(lidar_footprints_3857), size = 1, color = "black") +
      facet_wrap(~ Method) +
      theme_minimal() +
      labs(title = "Footprint Displacement Trajectories by Method",
           subtitle = "Black points: Original positions, Colored lines: Displacement vectors") +
      scale_color_brewer(palette = "Set1") +
      coord_sf(datum = NA)
    
    print(plot2)
  }
} else {
  cat("No valid trajectory data to plot\n")
}

## -----------------------------------------------------------------------------
# Function to calculate displacement distances
calc_displacements <- function(original_sf, corrected_sf, method_name) {
  # Extract centroids
  orig_points <- st_centroid(original_sf)
  corr_points <- st_centroid(corrected_sf)
  
  # Calculate distances
  orig_coords <- st_coordinates(orig_points)
  corr_coords <- st_coordinates(corr_points)
  
  dx <- corr_coords[, 1] - orig_coords[, 1]
  dy <- corr_coords[, 2] - orig_coords[, 2]
  distance <- sqrt(dx^2 + dy^2)
  
  # Return as data frame
  data.frame(
    Distance = distance,
    Method = method_name
  )
}

# Create displacement data frames
displacement_list <- list()

if (!is.null(ga_result$corrected_footprints)) {
  ga_displacements <- calc_displacements(
    lidar_footprints_3857,
    ga_result$corrected_footprints,
    "GA"
  )
  displacement_list[[1]] <- ga_displacements
}

if (!is.null(pso_result$corrected_footprints)) {
  pso_displacements <- calc_displacements(
    lidar_footprints_3857,
    pso_result$corrected_footprints,
    "PSO"
  )
  displacement_list[[2]] <- pso_displacements
}

if (!is.null(woa_result) && !is.null(woa_result$corrected_footprints)) {
  woa_displacements <- calc_displacements(
    lidar_footprints_3857,
    woa_result$corrected_footprints,
    "WOA"
  )
  displacement_list[[3]] <- woa_displacements
}

if (!is.null(lbfgsb_result) && !is.null(lbfgsb_result$corrected_footprints)) {
  lbfgsb_displacements <- calc_displacements(
    lidar_footprints_3857,
    lbfgsb_result$corrected_footprints,
    "L-BFGS-B"
  )
  displacement_list[[4]] <- lbfgsb_displacements
}

# Combine all displacements
if (length(displacement_list) > 0) {
  all_displacements <- do.call(rbind, displacement_list)
  
  # Create a density plot of displacement distances
  ggplot(all_displacements, aes(x = Distance, fill = Method)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    labs(title = "Distribution of Displacement Distances",
         x = "Distance (meters)",
         y = "Density") +
    scale_fill_brewer(palette = "Set1")
  
  # Create a boxplot of displacement distances
  ggplot(all_displacements, aes(x = Method, y = Distance, fill = Method)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = "Distribution of Displacement Distances",
         x = "Method",
         y = "Distance (meters)") +
    scale_fill_brewer(palette = "Set1")
} else {
  cat("No valid displacement data to plot")
}

## ----eval=FALSE---------------------------------------------------------------
# # Using different distance metrics
# dtw_result <- positional_correction(
#   lidar_footprints = lidar_footprints_3857,
#   input_rast = dem_rast,
#   minimizing_method = "dtw",  # Can be "euclidean", "dtw", "manhattan", "correlation", "area", or "hausdorff"
#   target_variable = "mean",
#   lidar_value = "gedi_dem",
#   optimization_method = "ga"
# )

