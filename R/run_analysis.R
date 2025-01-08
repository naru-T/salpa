#' Run the analysis for all orbits
#'
#' @param gedi_l2a sf object of gedi l2a
#' @param input_raster raster object for extraction
#' @param buf numeric value for buffer size
#' @param minimizing_method character string for the minimizing method
#' @param target_variable character string for the target variable
#' @param crs_code numeric value for CRS code
#' @param gedi_value numeric value for GEDI data
#' @param lower_bounds numeric vector for lower bounds in GA
#' @param upper_bounds numeric vector for upper bounds in GA
#' @param pop_size numeric value for population size in GA
#' @param max_iter numeric value for maximum iterations in GA
#' @param parallel logical value for parallel processing
#' @return list with results for all orbits
#' @export
run_analysis <- function(gedi_l2a, input_raster, minimizing_method = "dtw", target_variable, buf = 12.5, crs_code, gedi_value, lower_bounds = c(-30, -30), upper_bounds = c(30, 30), pop_size = 50, max_iter = 100, parallel = parallel) {
    optim_result <- tryCatch({
        results <- minimize_loss(gedi_l2a, input_raster, minimizing_method, target_variable, buf, crs_code, gedi_value, lower_bounds, upper_bounds, pop_size, max_iter, parallel)
        data.frame(best_x = results$best_x, best_y = results$best_y, best_value = results$best_value)
    },
    error = function(e) {
        data.frame(best_x = NA, best_y = NA, best_value = NA)
    })

    return(optim_result)
}