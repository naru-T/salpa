#' Positional correction of satellite LiDAR footprints
#'
#' This function performs positional correction on satellite LiDAR footprints using genetic algorithm optimization.
#'
#' @param lidar_footprints sf object containing satellite LiDAR footprints
#' @param input_rast raster object for reference data
#' @param minimizing_method character string for the minimizing method. Default is "euclidean".
#'        Options are "dtw", "euclidean", "manhattan", "correlation", "area", "hausdorff".
#' @param target_variable character string for the target variable
#' @param buf numeric value for buffer size (default: 12.5)
#' @param crs_code numeric value for CRS code
#' @param lidar_value name of the column containing LiDAR measurement values
#' @param lower_bounds numeric vector for lower bounds in GA (default: c(-30, -30))
#' @param upper_bounds numeric vector for upper bounds in GA (default: c(30, 30))
#' @param pop_size numeric value for population size in GA (default: 50)
#' @param max_iter numeric value for maximum iterations in GA (default: 100)
#' @param parallel logical value for parallel processing
#' @return list containing optimization results and corrected positions
#' @export
positional_correction <- function(lidar_footprints, input_rast, minimizing_method = "dtw", target_variable, buf = 12.5, crs_code, lidar_value, lower_bounds = c(-30, -30), upper_bounds = c(30, 30), pop_size = 50, max_iter = 100, parallel = FALSE) {
    # Validate buffer size
    if (buf <= 0) {
        stop("Buffer size must be a positive number")
    }
    
    optim_result <- tryCatch({
        results <- minimize_loss(lidar_footprints, input_rast, minimizing_method, target_variable, buf, crs_code, lidar_value, lower_bounds, upper_bounds, pop_size, max_iter, parallel)
        data.frame(best_x = results$best_x, best_y = results$best_y, best_value = results$best_value)
    },
    error = function(e) {
        data.frame(best_x = NA, best_y = NA, best_value = NA)
    })

    position_adjustment <- position_adjustment(lidar_footprints, optim_result$best_x, optim_result$best_y, crs_code)
    return(list(optim_result = optim_result, position_adjustment = position_adjustment))
}