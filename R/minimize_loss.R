#' Minimize loss
#'
#' @param gedi_l2a sf object of gedi_l2a
#' @param input_raster referenced raster object for extraction
#' @param minimizing_method character string for the minimizing method ("dtw" is only supported)
#' @param target_variable character string for the target variable (min, max, mean, median, etc.). See exact_extract
#' @param buf numeric value for buf size. Default is 12.5
#' @param crs_code numeric value for CRS code
#' @param gedi_value character string for the GEDI value
#' @param lower_bounds numeric vector for lower bounds in GA. Default is c(-30, -30)
#' @param upper_bounds numeric vector for upper bounds in GA. Default is c(30, 30)
#' @param pop_size numeric value for population size in GA. Default is 50
#' @param max_iter numeric value for maximum iterations in GA. Default is 100
#' @param parallel logical value for parallel processing. Default is FALSE
#' @importFrom GA ga
#' @return list with optimization results
#' @export
minimize_loss <- function(gedi_l2a, input_raster, minimizing_method, target_variable, buf, crs_code, gedi_value, lower_bounds, upper_bounds, pop_size, max_iter, parallel) {

    if (dim(gedi_l2a)[1] == 0) {
        print(paste("Skipping because of no high quality data found"))
        return(NULL)
    }

    objective_function <- function(x, y) {
        res <- get_loss(gedi_l2a = gedi_l2a, add_x = x, add_y = y, buf = buf, input_raster = input_raster, minimizing_method = minimizing_method, target_variable = target_variable, crs_code = crs_code, gedi_value = gedi_value)
        out <- res$loss_value
        if (is.na(out)) {
            print(paste("NA values produced at x: ", x, " y: ", y))
            out <- 1000
        }
        return(out)
    }


        fitness <- function(x) {
            # Calculate the fitness value with error handling
            tryCatch({
                fitness_value <- -objective_function(x[1], x[2])
                if(is.na(fitness_value) || !is.numeric(fitness_value)) {
                return(-Inf)  # Return -Inf for invalid solutions
                }
                return(fitness_value)
            }, error = function(e) {
                return(-Inf)  # Return -Inf for errors
            })
            }
        out <- ga(type = "real-valued",
                  fitness = fitness,
                  lower = lower_bounds,
                  upper = upper_bounds,
                  popSize = pop_size,
                  maxiter = max_iter,
                  parallel = parallel,
                  seed = 1118)

      optim_result <-  data.frame(
            best_x = (out@solution)[1],
            best_y = (out@solution)[2],
            best_value = -(out@fitnessValue))

    #print(paste("DONE... x:", optim_result$best_x, "y:", optim_result$best_y, "value:", optim_result$best_value))

    return(optim_result)
}
