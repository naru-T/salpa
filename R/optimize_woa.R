#' High-Performance Whale Optimization Algorithm
#'
#' @param objective_function Function to minimize
#' @param lower_bounds Lower bounds for optimization
#' @param upper_bounds Upper bounds for optimization
#' @param pop_size Population (whale pod) size
#' @param max_iter Maximum iterations
#' @param woa_params Additional WOA parameters
#' @param parallel Logical indicating whether to use parallel fitness evaluation
#' @return Data frame with optimization results
#' @noRd
optimize_woa <- function(objective_function, lower_bounds, upper_bounds,
                         pop_size, max_iter, woa_params = NULL, parallel = FALSE,
                         batch_fn = NULL) {

    # Enhanced WOA parameters with performance optimizations
    woa_control <- list(
        verbose              = FALSE,
        early_stop           = TRUE,
        convergence_threshold = 1e-8,
        stagnation_limit     = max_iter %/% 5,
        memoization          = TRUE
    )

    # Update with user-provided parameters
    if (!is.null(woa_params)) {
        for (param_name in names(woa_params)) {
            woa_control[[param_name]] <- woa_params[[param_name]]
        }
    }

    return(run_serial_woa(objective_function, lower_bounds, upper_bounds,
                          pop_size, max_iter, woa_control, parallel,
                          batch_fn = batch_fn))
}

#' Run WOA with C++ position updates and optional parallel fitness evaluation
#' @noRd
run_serial_woa <- function(objective_function, lower_bounds, upper_bounds,
                           pop_size, max_iter, woa_control, parallel = FALSE,
                           batch_fn = NULL) {

    # Create cached objective function
    if (woa_control$memoization) {
        cached_objective <- create_memoized_objective(objective_function)
    } else {
        cached_objective <- create_fast_woa_objective(objective_function)
    }

    n_dim <- length(lower_bounds)
    b     <- 1.0   # Spiral shape constant (Mirjalili 2016)

    # Initialise whale positions
    whales <- initialize_smart_whale_positions(pop_size, n_dim, lower_bounds, upper_bounds)

    # Evaluate initial fitness (batch path skips per-particle fork overhead)
    fitness <- if (!is.null(batch_fn)) batch_fn(whales) else
                   eval_population_fitness(whales, cached_objective, parallel)

    best_idx     <- which.min(fitness)
    best_pos     <- whales[best_idx, ]
    best_fitness <- fitness[best_idx]

    convergence_history <- numeric(max_iter)
    stagnation_count    <- 0L
    previous_best       <- Inf

    for (iter in seq_len(max_iter)) {

        # Bug 3 fix: a decreases linearly 2->0; a2 decreases -1->-2
        a  <- 2.0 - iter * (2.0 / max_iter)
        a2 <- -1.0 + iter * (-1.0 / max_iter)

        # Pre-generate all random numbers for this iteration
        r1_vec <- stats::runif(pop_size)
        r2_vec <- stats::runif(pop_size)
        # Bug 3 fix: l range is [(a2-1), 1], adapts adaptively
        l_vec    <- (a2 - 1.0) * stats::runif(pop_size) + 1.0
        p_vec    <- stats::runif(pop_size)
        # rand_idx: 0-indexed random whale indices for exploration phase
        rand_idx <- sample.int(pop_size, pop_size, replace = TRUE) - 1L

        # Bug 2 fix: correct spiral formula via C++ (b=1, cos(2*pi*l))
        whales <- woa_update_positions_cpp(
            whales, best_pos, rand_idx,
            r1_vec, r2_vec, l_vec, p_vec,
            a, b, 0.5,
            lower_bounds, upper_bounds
        )

        # Evaluate fitness (batch or parallel scalar)
        fitness <- if (!is.null(batch_fn)) batch_fn(whales) else
                       eval_population_fitness(whales, cached_objective, parallel)

        # Update best
        curr_best_idx <- which.min(fitness)
        if (fitness[curr_best_idx] < best_fitness) {
            best_idx     <- curr_best_idx
            best_pos     <- whales[best_idx, ]
            best_fitness <- fitness[curr_best_idx]
        }

        convergence_history[iter] <- best_fitness

        # Early stopping
        if (woa_control$early_stop) {
            improvement <- abs(previous_best - best_fitness)
            if (improvement < woa_control$convergence_threshold) {
                stagnation_count <- stagnation_count + 1L
            } else {
                stagnation_count <- 0L
            }
            if (stagnation_count >= woa_control$stagnation_limit) {
                if (woa_control$verbose) {
                    message(sprintf("WOA early stopping at iteration %d", iter))
                }
                break
            }
            previous_best <- best_fitness
        }

        if (woa_control$verbose &&
            (iter %% 20 == 0 || iter == 1L || iter == max_iter)) {
            cat(sprintf("WOA Iteration %d: Best fitness = %.8f\n", iter, best_fitness))
        }
    }

    data.frame(
        best_x     = best_pos[1],
        best_y     = best_pos[2],
        best_value = best_fitness
    )
}


#' Create memoized objective function for WOA
#' @noRd
create_memoized_objective <- function(objective_function) {
    cache <- new.env(hash = TRUE)

    function(x, y) {
        key <- paste(round(c(x, y), 6), collapse = ",")

        if (exists(key, envir = cache)) {
            return(get(key, envir = cache))
        }

        value <- tryCatch({
            result <- objective_function(x, y)
            if (is.na(result) || !is.numeric(result) || is.infinite(result)) 1000
            else result
        }, error = function(e) 1000)

        if (length(cache) < 5000) {
            assign(key, value, envir = cache)
        }

        return(value)
    }
}

#' Create fast WOA objective function (no caching)
#' @noRd
create_fast_woa_objective <- function(objective_function) {
    function(x, y) {
        tryCatch({
            value <- objective_function(x, y)
            if (is.na(value) || !is.numeric(value) || is.infinite(value)) 1000
            else value
        }, error = function(e) 1000)
    }
}


#' Initialise whale positions with mixed strategy
#' @noRd
initialize_smart_whale_positions <- function(pop_size, n_dim, lower_bounds, upper_bounds) {
    n_center  <- max(1L, pop_size %/% 4L)
    n_corners <- min(2^n_dim, pop_size %/% 4L)

    whales <- matrix(0, pop_size, n_dim)
    idx    <- 1L

    # Centre-based
    center <- (lower_bounds + upper_bounds) / 2
    for (i in seq_len(n_center)) {
        noise         <- stats::runif(n_dim, -0.2, 0.2) * (upper_bounds - lower_bounds)
        whales[idx, ] <- pmax(lower_bounds, pmin(upper_bounds, center + noise))
        idx           <- idx + 1L
    }

    # Corner-based (2-D only)
    if (n_corners > 0 && n_dim == 2L) {
        corners <- expand.grid(c(lower_bounds[1], upper_bounds[1]),
                               c(lower_bounds[2], upper_bounds[2]))
        for (i in seq_len(min(n_corners, nrow(corners)))) {
            whales[idx, ] <- as.numeric(corners[i, ])
            idx           <- idx + 1L
        }
    }

    # Random fill
    while (idx <= pop_size) {
        for (j in seq_len(n_dim)) {
            whales[idx, j] <- stats::runif(1, lower_bounds[j], upper_bounds[j])
        }
        idx <- idx + 1L
    }

    return(whales)
}
