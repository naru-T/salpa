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
#' @param parallel logical value for parallel processing. Default is FALSE
#' @param pso_params list of additional parameters for PSO. Default is NULL.
#'        List can include: trace (logical), maxit (numeric), maxf (numeric), abstol (numeric), reltol (numeric), REPORT (numeric), trace.stats (logical)
#' @param woa_params list of additional parameters for WOA. Default is NULL.
#'        List can include: trace (logical), batch_size (numeric for batch processing)
#' @param lbfgsb_params list of additional parameters for L-BFGS-B. Default is NULL.
#'        List can include: factr (numeric), pgtol (numeric), trace (numeric), REPORT (integer)
#' @importFrom GA ga
#' @importFrom stats na.omit
#' @return list with optimization results
#' @export
minimize_loss <- function(lidar_footprints, input_rast, minimizing_method, target_variable, buf, crs_code, lidar_value, 
                          optimization_method = "ga", lower_bounds = c(-30, -30), upper_bounds = c(30, 30), 
                          pop_size = 50, max_iter = 100, parallel = FALSE, pso_params = NULL, woa_params = NULL,
                          lbfgsb_params = NULL) {

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
            
    } else if (optimization_method == "woa") {
        # Whale Optimization Algorithm (WOA) implementation
        
        # Default WOA parameters
        woa_control <- list(
            verbose = FALSE,          # Don't print progress
            batch_size = pop_size     # Process all in one batch by default
        )
        
        # Update with user-provided parameters if any
        if (!is.null(woa_params)) {
            for (param_name in names(woa_params)) {
                woa_control[[param_name]] <- woa_params[[param_name]]
            }
        }
        
        # Try standard optimization with L-BFGS-B first for small search spaces
        # This is often much faster for well-behaved functions
        if (all(abs(upper_bounds - lower_bounds) <= 60)) {
            # Standard objective function
            woa_objective <- function(x) {
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
            
            # Try standard optimization first
            optim_result_std <- try({
                optim_std <- stats::optim(
                    par = c(0, 0),
                    fn = woa_objective,
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
            # If standard optimization failed, fall back to WOA
        }
        
        # Batch evaluation function
        batch_evaluate <- function(positions) {
            n_positions <- nrow(positions)
            results <- numeric(n_positions)
            
            for (i in 1:n_positions) {
                tryCatch({
                    results[i] <- objective_function(positions[i, 1], positions[i, 2])
                    if(is.na(results[i]) || !is.numeric(results[i])) {
                        results[i] <- 1000  # Penalty for invalid solutions
                    }
                }, error = function(e) {
                    results[i] <- 1000  # Penalty for errors
                })
            }
            
            return(results)
        }
        
        # Dimensions of the problem
        n_dim <- length(lower_bounds)
        
        # Initialize whale positions randomly within bounds
        whales <- matrix(runif(pop_size * n_dim), pop_size, n_dim)
        for (i in 1:n_dim) {
            whales[, i] <- lower_bounds[i] + whales[, i] * (upper_bounds[i] - lower_bounds[i])
        }
        
        # Create batches for evaluation
        batch_size <- min(woa_control$batch_size, pop_size)
        batch_indices <- split(1:pop_size, ceiling(seq_along(1:pop_size) / batch_size))
        
        # Initialize fitness values
        fitness <- numeric(pop_size)
        
        # Evaluate initial positions in batches
        for (b in 1:length(batch_indices)) {
            indices <- batch_indices[[b]]
            batch_positions <- whales[indices, , drop = FALSE]
            fitness[indices] <- batch_evaluate(batch_positions)
        }
        
        # Find the best whale
        best_idx <- which.min(fitness)
        best_pos <- whales[best_idx, ]
        best_fitness <- fitness[best_idx]
        
        # Initialize convergence tracking
        convergence <- numeric(max_iter)
        convergence[1] <- best_fitness
        
        # Main WOA loop
        for (iter in 1:max_iter) {
            # Update parameters
            a <- 2 - iter * (2 / max_iter)  # Linearly decreased from 2 to 0
            
            # Create matrix for new positions
            new_whales <- matrix(0, pop_size, n_dim)
            
            # Update each whale's position
            for (i in 1:pop_size) {
                # Random parameters
                r1 <- runif(1)
                r2 <- runif(1)
                A <- 2 * a * r1 - a
                C <- 2 * r2
                l <- runif(1, -1, 1)
                p <- runif(1)
                
                # Either follow prey or search randomly
                if (p < 0.5) {
                    # Follow prey (exploitation)
                    if (abs(A) < 1) {
                        D <- abs(C * best_pos - whales[i, ])
                        new_whales[i, ] <- best_pos - A * D
                    }
                    # Search randomly (exploration)
                    else {
                        random_idx <- sample(1:pop_size, 1)
                        random_whale <- whales[random_idx, ]
                        D <- abs(C * random_whale - whales[i, ])
                        new_whales[i, ] <- random_whale - A * D
                    }
                }
                # Bubble-net attacking (spiral updating)
                else {
                    D <- abs(best_pos - whales[i, ])
                    new_whales[i, ] <- D * exp(l * 2 * pi) * cos(2 * pi * l) + best_pos
                }
                
                # Enforce bounds
                for (j in 1:n_dim) {
                    new_whales[i, j] <- max(lower_bounds[j], min(upper_bounds[j], new_whales[i, j]))
                }
            }
            
            # Replace old positions with new ones
            whales <- new_whales
            
            # Evaluate new positions in batches
            for (b in 1:length(batch_indices)) {
                indices <- batch_indices[[b]]
                batch_positions <- whales[indices, , drop = FALSE]
                fitness[indices] <- batch_evaluate(batch_positions)
            }
            
            # Update best whale
            curr_best_idx <- which.min(fitness)
            if (fitness[curr_best_idx] < best_fitness) {
                best_idx <- curr_best_idx
                best_pos <- whales[best_idx, ]
                best_fitness <- fitness[best_idx]
            }
            
            # Track convergence
            convergence[iter] <- best_fitness
            
            # Optional progress output
            if (woa_control$verbose && (iter %% 10 == 0 || iter == 1 || iter == max_iter)) {
                cat("WOA Iteration", iter, ": Best fitness =", best_fitness, "\n")
            }
        }
        
        # Return results in the expected format
        optim_result <- data.frame(
            best_x = best_pos[1],
            best_y = best_pos[2],
            best_value = best_fitness)
            
    } else if (optimization_method == "lbfgsb") {
        # L-BFGS-B optimization (direct)
        
        # Standard objective function
        lbfgs_objective <- function(x) {
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
        
        # Default control parameters
        control_params <- list(
            maxit = max_iter,     # Maximum iterations
            factr = 1e7,          # Default convergence factor
            pgtol = 1e-5,         # Gradient projection tolerance
            trace = 0             # No tracing
        )
        
        # Update with user-provided parameters if any
        if (!is.null(lbfgsb_params)) {
            for (param_name in names(lbfgsb_params)) {
                control_params[[param_name]] <- lbfgsb_params[[param_name]]
            }
        }
        
        # Run L-BFGS-B optimization
        optim_result_lbfgs <- stats::optim(
            par = c(0, 0),              # Starting position
            fn = lbfgs_objective,       # Objective function
            method = "L-BFGS-B",        # Method
            lower = lower_bounds,       # Lower bounds
            upper = upper_bounds,       # Upper bounds
            control = control_params    # Control parameters
        )
        
        # Return results in the expected format
        optim_result <- data.frame(
            best_x = optim_result_lbfgs$par[1],
            best_y = optim_result_lbfgs$par[2],
            best_value = optim_result_lbfgs$value)
            
    } else {
        stop("Invalid optimization_method. Use 'ga', 'pso', 'woa', or 'lbfgsb'.")
    }

    #print(paste("DONE... x:", optim_result$best_x, "y:", optim_result$best_y, "value:", optim_result$best_value))

    return(optim_result)
}
