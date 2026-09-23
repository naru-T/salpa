#' High-Performance L-BFGS-B Optimization
#'
#' @param objective_function Function to minimize
#' @param lower_bounds Lower bounds for optimization
#' @param upper_bounds Upper bounds for optimization
#' @param max_iter Maximum iterations
#' @param lbfgsb_params Additional L-BFGS-B parameters
#' @param parallel Logical indicating whether to use parallel processing (for multi-start only)
#' @return Data frame with optimization results
#' @noRd
optimize_lbfgsb <- function(objective_function, lower_bounds, upper_bounds, 
                           max_iter, lbfgsb_params = NULL, parallel = FALSE) {
    
    # Single-start L-BFGS-B: parallel multi-start is theoretically safe (terra::wrap
    # since v0.0.1.8.0) but IO-bound workloads on external storage suffer from
    # bandwidth contention when multiple L-BFGS-B processes read the DEM in parallel.
    lbfgs_control <- list(
        multi_start        = FALSE,
        n_starts           = 1L,
        adaptive_tolerance = FALSE,
        gradient_check     = FALSE,
        memoization        = FALSE,
        early_stopping     = FALSE,
        parallel_starts    = FALSE
    )
    
    # Update with user-provided parameters
    if (!is.null(lbfgsb_params)) {
        for (param_name in names(lbfgsb_params)) {
            lbfgs_control[[param_name]] <- lbfgsb_params[[param_name]]
        }
    }
    
    # Run enhanced L-BFGS-B optimization
    if (lbfgs_control$multi_start) {
        return(run_multi_start_lbfgsb(objective_function, lower_bounds, upper_bounds, 
                                     max_iter, lbfgs_control))
    } else {
        return(run_single_lbfgsb(objective_function, lower_bounds, upper_bounds, 
                                max_iter, lbfgs_control))
    }
}

#' Run single L-BFGS-B optimization with enhancements
#' @noRd
run_single_lbfgsb <- function(objective_function, lower_bounds, upper_bounds, 
                             max_iter, lbfgs_control) {
    
    # Create simple objective function
    enhanced_objective <- create_fast_lbfgs_objective(objective_function)
    
    # Simplified control parameters with practical defaults
    control_params <- list(
        maxit = max_iter,
        factr = 1e7,          # Standard tolerance
        pgtol = 1e-4,         # Reasonable gradient tolerance
        trace = 0,            # No tracing for speed
        REPORT = 10,          # Report every 10 iterations
        fnscale = 1,          # Function scaling
        parscale = rep(1, 2), # Parameter scaling  
        ndeps = rep(1e-4, 2)  # Finite difference step sizes
    )
    
    # Smart starting point selection
    start_point <- generate_smart_start_point(lower_bounds, upper_bounds)
    
    # Run L-BFGS-B optimization with enhanced error handling
    tryCatch({
        optim_result <- stats::optim(
            par = start_point,
            fn = enhanced_objective,
            method = "L-BFGS-B",
            lower = lower_bounds,
            upper = upper_bounds,
            control = control_params
        )
        
        # Validate result
        if (optim_result$convergence != 0 && optim_result$convergence != 1) {
            warning(sprintf("L-BFGS-B convergence issue (code %d): %s", 
                           optim_result$convergence, 
                           get_convergence_message(optim_result$convergence)))
        }
        
        data.frame(
            best_x = as.numeric(optim_result$par[1]),
            best_y = as.numeric(optim_result$par[2]),
            best_value = as.numeric(optim_result$value)
        )
        
    }, error = function(e) {
        warning("L-BFGS-B optimization failed: ", e$message)
        
        # Fallback to simpler L-BFGS-B
        fallback_result <- stats::optim(
            par = c(0, 0),
            fn = enhanced_objective,
            method = "L-BFGS-B",
            lower = lower_bounds,
            upper = upper_bounds,
            control = list(maxit = max_iter %/% 2)
        )
        
        data.frame(
            best_x = as.numeric(fallback_result$par[1]),
            best_y = as.numeric(fallback_result$par[2]),
            best_value = as.numeric(fallback_result$value)
        )
    })
}

#' Run multi-start L-BFGS-B optimization
#' @noRd
run_multi_start_lbfgsb <- function(objective_function, lower_bounds, upper_bounds, 
                                  max_iter, lbfgs_control) {
    
    # Create enhanced objective function
    if (lbfgs_control$memoization) {
        enhanced_objective <- create_memoized_lbfgs_objective(objective_function)
    } else {
        enhanced_objective <- create_fast_lbfgs_objective(objective_function)
    }
    
    n_starts <- lbfgs_control$n_starts
    
    # Generate diverse starting points
    start_points <- generate_diverse_start_points(n_starts, lower_bounds, upper_bounds)
    
    # Adaptive control parameters
    search_range <- max(abs(upper_bounds - lower_bounds))
    
    control_params <- list(
        maxit = max_iter %/% n_starts,  # Distribute iterations across starts
        factr = max(1e6, min(1e10, 1e7 * (search_range / 50))),
        pgtol = max(1e-7, min(1e-4, 1e-5 * sqrt(search_range / 50))),
        trace = 0,
        REPORT = 10  # Report every 10 iterations for multi-start
    )
    
    # Run optimizations from multiple starting points
    if (lbfgs_control$parallel_starts && requireNamespace("parallel", quietly = TRUE)) {
        results <- run_parallel_multi_start(enhanced_objective, start_points, 
                                           lower_bounds, upper_bounds, control_params)
    } else {
        results <- run_serial_multi_start(enhanced_objective, start_points, 
                                         lower_bounds, upper_bounds, control_params)
    }
    
    # Select best result
    values <- sapply(results, function(r) {
        if (is.null(r) || is.null(r$value)) return(1000)
        return(r$value)
    })
    best_idx <- which.min(values)
    best_result <- results[[best_idx]]
    
    # Ensure we have valid results
    if (is.null(best_result) || is.null(best_result$par)) {
        return(data.frame(
            best_x = 0,
            best_y = 0, 
            best_value = 1000
        ))
    }
    
    data.frame(
        best_x = as.numeric(best_result$par[1]),
        best_y = as.numeric(best_result$par[2]),
        best_value = as.numeric(best_result$value)
    )
}

#' Run serial multi-start optimization
#' @noRd
run_serial_multi_start <- function(objective_function, start_points, 
                                  lower_bounds, upper_bounds, control_params) {
    
    results <- vector("list", nrow(start_points))
    
    for (i in 1:nrow(start_points)) {
        tryCatch({
            results[[i]] <- stats::optim(
                par = start_points[i, ],
                fn = objective_function,
                method = "L-BFGS-B",
                lower = lower_bounds,
                upper = upper_bounds,
                control = control_params
            )
        }, error = function(e) {
            # If one start fails, use a fallback result
            results[[i]] <- list(
                par = start_points[i, ],
                value = 1000,
                convergence = 999
            )
        })
    }
    
    return(results)
}

#'  Run parallel multi-start optimization using mclapply (fork-based, works with closures)
#' @noRd
run_parallel_multi_start <- function(objective_function, start_points,
                                     lower_bounds, upper_bounds, control_params) {

    num_cores <- min(max(1L, parallel::detectCores() - 1L), nrow(start_points))

    # mclapply is fork-based: inherits parent environment including raster closures.
    # Falls back to serial on Windows (where forking is unavailable).
    parallel::mclapply(seq_len(nrow(start_points)), function(i) {
        tryCatch({
            stats::optim(
                par     = start_points[i, ],
                fn      = objective_function,
                method  = "L-BFGS-B",
                lower   = lower_bounds,
                upper   = upper_bounds,
                control = control_params
            )
        }, error = function(e) {
            list(par = start_points[i, ], value = 1000, convergence = 999)
        })
    }, mc.cores = num_cores)
}

#' Create memoized L-BFGS-B objective function
#' @noRd
create_memoized_lbfgs_objective <- function(objective_function) {
    cache <- new.env(hash = TRUE)
    
    function(x) {
        # Create cache key
        key <- paste(round(x, 8), collapse = ",")  # Higher precision for L-BFGS-B
        
        if (exists(key, envir = cache)) {
            return(get(key, envir = cache))
        }
        
        # Compute objective value
        value <- tryCatch({
            result <- objective_function(x[1], x[2])
            if(is.na(result) || !is.numeric(result) || is.infinite(result)) {
                1000
            } else {
                result
            }
        }, error = function(e) 1000)
        
        # Cache the result (with size limit)
        if (length(cache) < 500) {  # Smaller cache for L-BFGS-B
            assign(key, value, envir = cache)
        }
        
        return(value)
    }
}

#' Create fast L-BFGS-B objective function
#' @noRd
create_fast_lbfgs_objective <- function(objective_function) {
    function(x) {
        tryCatch({
            value <- objective_function(x[1], x[2])
            if(is.na(value) || !is.numeric(value) || is.infinite(value)) {
                return(1000)
            }
            return(value)
        }, error = function(e) {
            return(1000)
        })
    }
}

#' Generate smart starting point for L-BFGS-B
#' @noRd
generate_smart_start_point <- function(lower_bounds, upper_bounds) {
    # Use center of search space as default starting point
    center <- (lower_bounds + upper_bounds) / 2
    
    # Add small perturbation to avoid potential singularities
    noise <- stats::runif(length(center), -0.05, 0.05) * (upper_bounds - lower_bounds)
    
    # Ensure within bounds
    start_point <- pmax(lower_bounds, pmin(upper_bounds, center + noise))
    
    return(start_point)
}

#' Generate diverse starting points for multi-start
#' @noRd
generate_diverse_start_points <- function(n_starts, lower_bounds, upper_bounds) {
    n_dim <- length(lower_bounds)
    start_points <- matrix(0, nrow = n_starts, ncol = n_dim)
    
    # Strategy 1: Center point
    center <- (lower_bounds + upper_bounds) / 2
    start_points[1, ] <- center
    
    if (n_starts > 1) {
        # Strategy 2: Corners (for 2D problems)
        if (n_dim == 2 && n_starts >= 2) {
            corners <- expand.grid(c(lower_bounds[1], upper_bounds[1]), 
                                  c(lower_bounds[2], upper_bounds[2]))
            for (i in 2:min(n_starts, nrow(corners) + 1)) {
                start_points[i, ] <- as.numeric(corners[i-1, ])
            }
        }
        
        # Strategy 3: Quasi-random points using Sobol sequence-like approach
        remaining_starts <- n_starts - min(n_starts, 5)
        if (remaining_starts > 0) {
            for (i in (n_starts - remaining_starts + 1):n_starts) {
                # Use different random seeds for diversity
                set.seed(i * 123)
                for (j in 1:n_dim) {
                    start_points[i, j] <- stats::runif(1, lower_bounds[j], upper_bounds[j])
                }
            }
        }
        
        # Strategy 4: Boundary-focused points
        if (n_starts > 5) {
            boundary_starts <- min(remaining_starts, 2)
            for (i in (n_starts - boundary_starts + 1):n_starts) {
                # Points closer to boundaries
                for (j in 1:n_dim) {
                    if (stats::runif(1) < 0.5) {
                        # Closer to lower bound
                        start_points[i, j] <- lower_bounds[j] + 
                            stats::runif(1, 0, 0.3) * (upper_bounds[j] - lower_bounds[j])
                    } else {
                        # Closer to upper bound
                        start_points[i, j] <- upper_bounds[j] - 
                            stats::runif(1, 0, 0.3) * (upper_bounds[j] - lower_bounds[j])
                    }
                }
            }
        }
    }
    
    return(start_points)
}

#' Get convergence message for L-BFGS-B
#' @noRd
get_convergence_message <- function(convergence_code) {
    switch(as.character(convergence_code),
           "0" = "Successful convergence",
           "1" = "Maximum iterations reached",
           "10" = "Degeneracy in Nelder-Mead simplex",
           "51" = "Warning from L-BFGS-B routine",
           "52" = "Error from L-BFGS-B routine",
           paste("Unknown convergence code:", convergence_code)
    )
}