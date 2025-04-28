#' Positional correction of satellite LiDAR footprints
#'
#' This function performs positional correction on satellite LiDAR footprints using optimization algorithms.
#'
#' @param lidar_footprints sf object containing satellite LiDAR footprints.
#' @param input_rast Raster object for reference data.
#' @param minimizing_method Character string for the minimizing method. Default is "euclidean".
#'        Options include "dtw", "euclidean", "manhattan", "correlation", "area", "hausdorff".
#' @param target_variable Character string for the target variable.
#' @param buf Numeric value for buffer size (default: 12.5).
#' @param crs_code Numeric value for CRS code.
#' @param lidar_value Name of the column containing LiDAR measurement values.
#' @param optimization_method Character string for optimization method. Options are "ga" (genetic algorithm) or "pso" (particle swarm optimization). Default is "ga".
#' @param lower_bounds Numeric vector for lower bounds in the optimization (default: c(-30, -30)).
#' @param upper_bounds Numeric vector for upper bounds in the optimization (default: c(30, 30)).
#' @param pop_size Numeric value for population size (GA) or swarm size (PSO) (default: 50).
#' @param max_iter Numeric value for maximum iterations (default: 100).
#' @param parallel Logical value for parallel processing (default: FALSE).
#' @param pso_params List of additional parameters for PSO (default: NULL).
#'        List can include: w (inertia), c1 (cognitive learning factor), c2 (social learning factor).
#' @param geographic_coordinates Logical value indicating if input coordinates are in geographic coordinates (default: FALSE).
#'        If TRUE, the bounds will be converted from meters to degrees.
#' @return List containing optimization results and corrected positions.
#' @export
positional_correction <- function(lidar_footprints, input_rast, minimizing_method = "euclidean",
                                    target_variable, buf = 12.5, crs_code, lidar_value,
                                    optimization_method = "ga",
                                    lower_bounds = c(-25, -25),
                                    upper_bounds = c(25, 25),
                                    pop_size = 50,
                                    max_iter = 100,
                                    parallel = FALSE,
                                    pso_params = NULL,
                                    geographic_coordinates = FALSE) {
    # Validate buffer size
    if (buf <= 0) {
        stop("Buffer size must be a positive number")
    }
    
    # Validate optimization method
    if (!(optimization_method %in% c("ga", "pso"))) {
        stop("optimization_method must be either 'ga' or 'pso'")
    }

    # Convert bounds to degrees if using geographic coordinates
    if (geographic_coordinates) {
        lower_bounds <- lower_bounds / 111320  # Convert meters to degrees (approximate)
        upper_bounds <- upper_bounds / 111320  # Convert meters to degrees (approximate)
    }

    optim_result <- tryCatch({
        results <- minimize_loss(lidar_footprints, input_rast, minimizing_method,
                                 target_variable, buf, crs_code, lidar_value,
                                 optimization_method, lower_bounds, upper_bounds, 
                                 pop_size, max_iter, parallel, pso_params)
        
        # If best_value is 1000 (worst case), set offsets to 0
        if (!is.null(results) && results$best_value == 1000) {
            data.frame(best_x = 0, best_y = 0, best_value = results$best_value)
        } else {
            data.frame(best_x = results$best_x, best_y = results$best_y, best_value = results$best_value)
        }
    },
    error = function(e) {
        warning("Loss minimization failed: ", e$message)
        data.frame(best_x = 0, best_y = 0, best_value = NA)
    })

    # Use default offsets (0, 0) if optimization failed
    if (is.na(optim_result$best_x) || is.na(optim_result$best_y)) {
        optim_result$best_x <- 0
        optim_result$best_y <- 0
    }

    # Compute the adjusted footprints using the position_adjustment function
    adjusted_footprints <- position_adjustment(lidar_footprints,
                                               optim_result$best_x,
                                               optim_result$best_y,
                                               crs_code)

    # Return the results with the adjusted sf object under the key 'position_adjustment'
    return(list(optim_result = optim_result, position_adjustment = adjusted_footprints))
}