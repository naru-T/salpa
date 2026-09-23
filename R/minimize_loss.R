#' Minimize loss
#'
#' @param lidar_footprints sf object containing satellite LiDAR footprints
#' @param input_rast referenced raster object for extraction
#' @param minimizing_method character string for the minimizing method. Default is "euclidean".
#'        Options are "dtw", "euclidean", "manhattan", "correlation", "area", "hausdorff".
#' @param target_variable character string for the target variable (min, max, mean, median, etc.). See exact_extract
#' @param buf numeric value for buf size. Default is 12.5
#' @param crs_code numeric value for CRS code
#' @param lidar_value name of the column containing LiDAR measurement values
#' @param optimization_method character string for optimization method. Options are "ga" (genetic algorithm), "pso" (particle swarm optimization), "woa" (whale optimization algorithm), or "lbfgsb" (L-BFGS-B). Default is "ga".
#' @param lower_bounds numeric vector for lower bounds in optimization. Default is c(-30, -30)
#' @param upper_bounds numeric vector for upper bounds in optimization. Default is c(30, 30)
#' @param pop_size numeric value for population size (GA) or swarm size (PSO/WOA). Default is 50
#' @param max_iter numeric value for maximum iterations. Default is 100
#' @param pso_params list of additional parameters for PSO. Default is NULL.
#'        List can include: trace (logical), maxit (numeric), maxf (numeric), abstol (numeric), reltol (numeric), REPORT (numeric), trace.stats (logical)
#' @param woa_params list of additional parameters for WOA. Default is NULL.
#'        List can include: trace (logical), batch_size (numeric for batch processing)
#' @param lbfgsb_params list of additional parameters for L-BFGS-B. Default is NULL.
#'        List can include: factr (numeric), pgtol (numeric), trace (numeric), REPORT (integer)
#' @param parallel logical value indicating whether to use parallel processing. Default is FALSE.
#' @importFrom stats na.omit
#' @return list with `best_x`, `best_y`, and `best_value` elements
#' @examples
#' \donttest{
#'   lidar_sf <- sf::st_as_sf(
#'     data.frame(x = c(500, 510, 520), y = c(300, 310, 320), elev = c(100, 101, 99)),
#'     coords = c("x", "y"), crs = 32654
#'   )
#'   dem <- terra::rast(
#'     nrows = 10, ncols = 10,
#'     xmin = 450, xmax = 600, ymin = 250, ymax = 400,
#'     crs = "EPSG:32654", vals = runif(100, 98, 103)
#'   )
#'   result <- minimize_loss(lidar_sf, dem, minimizing_method = "euclidean",
#'                           target_variable = "mean", buf = 10, crs_code = 32654,
#'                           lidar_value = "elev", optimization_method = "lbfgsb",
#'                           max_iter = 10)
#' }
#' @export
minimize_loss <- function(lidar_footprints, input_rast, minimizing_method, target_variable, buf, crs_code, lidar_value, 
                          optimization_method = "ga", lower_bounds = c(-30, -30), upper_bounds = c(30, 30), 
                          pop_size = 50, max_iter = 100, pso_params = NULL, woa_params = NULL,
                          lbfgsb_params = NULL, parallel = FALSE) {

    # Check if there's data to process
    if (dim(lidar_footprints)[1] == 0) {
        message("Skipping because of no high quality data found")
        return(NULL)
    }

    # Create scalar objective function (used by lbfgsb and as fallback)
    objective_function <- create_objective_function(
        lidar_footprints = lidar_footprints,
        input_rast = input_rast,
        buf = buf,
        minimizing_method = minimizing_method,
        target_variable = target_variable,
        crs_code = crs_code,
        lidar_value = lidar_value
    )

    # GA/PSO/WOA use eval_population_fitness (mclapply parallel) which is faster
    # than batch_fn on real data (benchmark v0.0.2.1.0: GA 0.6x, PSO 0.32x, WOA 0.27x slower with batch_fn).
    # L-BFGS-B uses gradient descent so batch_fn is not applicable.
    batch_fn <- NULL

    # Run optimization using factory pattern
    tryCatch({
        optimization_factory(
            optimization_method = optimization_method,
            objective_function = objective_function,
            lower_bounds = lower_bounds,
            upper_bounds = upper_bounds,
            pop_size = pop_size,
            max_iter = max_iter,
            pso_params = pso_params,
            woa_params = woa_params,
            lbfgsb_params = lbfgsb_params,
            parallel = parallel,
            batch_fn = batch_fn
        )
    }, error = function(e) {
        handle_optimization_failure(e)
    })
}