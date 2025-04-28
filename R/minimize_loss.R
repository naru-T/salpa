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
#' @param optimization_method character string for optimization method. Options are "ga" (genetic algorithm) or "pso" (particle swarm optimization). Default is "ga".
#' @param lower_bounds numeric vector for lower bounds in optimization. Default is c(-30, -30)
#' @param upper_bounds numeric vector for upper bounds in optimization. Default is c(30, 30)
#' @param pop_size numeric value for population size (GA) or swarm size (PSO). Default is 50
#' @param max_iter numeric value for maximum iterations. Default is 100
#' @param parallel logical value for parallel processing. Default is FALSE
#' @param pso_params list of additional parameters for PSO. Default is NULL.
#'        List can include: trace (logical), maxit (numeric), maxf (numeric), abstol (numeric), reltol (numeric), REPORT (numeric), trace.stats (logical)
#' @importFrom GA ga
#' @importFrom stats na.omit
#' @return list with optimization results
#' @export
minimize_loss <- function(lidar_footprints, input_rast, minimizing_method, target_variable, buf, crs_code, lidar_value, 
                          optimization_method = "ga", lower_bounds = c(-30, -30), upper_bounds = c(30, 30), 
                          pop_size = 50, max_iter = 100, parallel = FALSE, pso_params = NULL) {

    if (dim(lidar_footprints)[1] == 0) {
        print(paste("Skipping because of no high quality data found"))
        return(NULL)
    }

    objective_function <- function(x, y) {
        res <- get_loss(lidar_footprints = lidar_footprints, add_x = x, add_y = y, buf = buf, input_rast = input_rast, minimizing_method = minimizing_method, target_variable = target_variable, crs_code = crs_code, lidar_value = lidar_value)
        out <- res$loss_value
        if (is.na(out)) {
            print(paste("NA values produced at x: ", x, " y: ", y))
            out <- 1000
        }
        return(out)
    }

    if (optimization_method == "ga") {
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
                  
        optim_result <- data.frame(
            best_x = (out@solution)[1],
            best_y = (out@solution)[2],
            best_value = -(out@fitnessValue))
            
    } else if (optimization_method == "pso") {
        # Load the pso package if not already loaded
        if (!requireNamespace("pso", quietly = TRUE)) {
            stop("Package 'pso' is needed for PSO optimization. Please install it with install.packages('pso')")
        }
        
        # If parallel is TRUE, set up parallel processing for PSO if possible
        if (parallel) {
            if (!requireNamespace("parallel", quietly = TRUE)) {
                warning("Package 'parallel' is needed for parallel PSO. Installing now...")
                install.packages("parallel")
            }
            
            # Set up parallel evaluation for PSO
            num_cores <- parallel::detectCores() - 1
            if (num_cores < 1) num_cores <- 1
            
            cl <- parallel::makeCluster(num_cores)
            on.exit(parallel::stopCluster(cl), add = TRUE)
            
            # Export necessary objects to cluster
            parallel::clusterExport(cl, varlist = c("objective_function", 
                                                   "lidar_footprints", 
                                                   "buf", 
                                                   "input_rast", 
                                                   "minimizing_method", 
                                                   "target_variable", 
                                                   "crs_code", 
                                                   "lidar_value"), 
                                  envir = environment())
            
            # Create parallel version of objective function
            pso_objective <- function(x) {
                tryCatch({
                    value <- objective_function(x[1], x[2])
                    if(is.na(value) || !is.numeric(value)) {
                        return(1000)  # Return high value for invalid solutions
                    }
                    return(value)
                }, error = function(e) {
                    return(1000)  # Return high value for errors
                })
            }
            
            # Define parallel evaluation function
            parallel_eval <- function(x) {
                parallel::parLapply(cl, as.list(1:nrow(x)), function(i) {
                    pso_objective(x[i,])
                })
            }
            
        } else {
            # Standard PSO objective function (PSO minimizes, so no need to negate)
            pso_objective <- function(x) {
                tryCatch({
                    value <- objective_function(x[1], x[2])
                    if(is.na(value) || !is.numeric(value)) {
                        return(1000)  # Return high value for invalid solutions
                    }
                    return(value)
                }, error = function(e) {
                    return(1000)  # Return high value for errors
                })
            }
        }
        
        # Default PSO control parameters - tuned for better performance
        control_params <- list(
            maxit = max_iter,     # Maximum iterations
            s = pop_size,         # Swarm size
            trace = FALSE,        # Don't print progress
            abstol = 1e-4,        # Absolute tolerance
            reltol = 1e-4,        # Relative tolerance
            maxf = max_iter * pop_size * 2  # Maximum function evaluations
        )
        
        # Update with user-provided parameters if any
        if (!is.null(pso_params)) {
            for (param_name in names(pso_params)) {
                control_params[[param_name]] <- pso_params[[param_name]]
            }
        }
        
        # Use optim with method="L-BFGS-B" for small search spaces as a comparison
        # This can be much faster for well-behaved functions
        if (all(abs(upper_bounds - lower_bounds) <= 60)) {
            # Try standard optimization first (often faster for simple problems)
            optim_result_std <- try({
                optim_std <- stats::optim(
                    par = c(0, 0),
                    fn = pso_objective,
                    method = "L-BFGS-B",
                    lower = lower_bounds,
                    upper = upper_bounds,
                    control = list(maxit = max_iter)
                )
                
                data.frame(
                    best_x = optim_std$par[1],
                    best_y = optim_std$par[2],
                    best_value = optim_std$value
                )
            }, silent = TRUE)
            
            if (!inherits(optim_result_std, "try-error")) {
                # Standard optimization worked, use its result
                return(optim_result_std)
            }
            # If standard optimization failed, fall back to PSO
        }
        
        # Run PSO optimization
        if (parallel) {
            # Use vectorized version for parallel execution
            pso_result <- pso::psoptim(
                par = c(0, 0),              # Starting position
                fn = pso_objective,         # Objective function
                lower = lower_bounds,       # Lower bounds
                upper = upper_bounds,       # Upper bounds
                control = control_params,   # Control parameters
                fnVectorized = TRUE,
                fnPar = parallel_eval
            )
        } else {
            # Standard serial execution
            pso_result <- pso::psoptim(
                par = c(0, 0),              # Starting position
                fn = pso_objective,         # Objective function
                lower = lower_bounds,       # Lower bounds
                upper = upper_bounds,       # Upper bounds
                control = control_params    # Control parameters
            )
        }
        
        optim_result <- data.frame(
            best_x = pso_result$par[1],
            best_y = pso_result$par[2],
            best_value = pso_result$value)
            
    } else {
        stop("Invalid optimization_method. Use 'ga' or 'pso'.")
    }

    #print(paste("DONE... x:", optim_result$best_x, "y:", optim_result$best_y, "value:", optim_result$best_value))

    return(optim_result)
}
