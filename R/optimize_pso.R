#' High-Performance Particle Swarm Optimization
#'
#' @param objective_function Function to minimize
#' @param lower_bounds Lower bounds for optimization
#' @param upper_bounds Upper bounds for optimization
#' @param pop_size Population (swarm) size
#' @param max_iter Maximum iterations
#' @param pso_params Additional PSO parameters
#' @param parallel Logical indicating whether to use parallel fitness evaluation
#' @return Data frame with optimization results
#' @noRd
optimize_pso <- function(objective_function, lower_bounds, upper_bounds,
                         pop_size, max_iter, pso_params = NULL, parallel = FALSE,
                         batch_fn = NULL) {

    return(run_enhanced_pso(objective_function, lower_bounds, upper_bounds,
                            pop_size, max_iter, pso_params, parallel = parallel,
                            batch_fn = batch_fn))
}

#' Internal PSO loop using C++ velocity/position updates
#' @noRd
run_enhanced_pso <- function(objective_function, lower_bounds, upper_bounds,
                             pop_size, max_iter, pso_params, parallel = FALSE,
                             batch_fn = NULL) {

    # PSO hyperparameters
    w   <- 0.7   # Inertia
    c1  <- 1.5   # Cognitive coefficient
    c2  <- 1.5   # Social coefficient
    search_range <- abs(upper_bounds - lower_bounds)
    v_max <- max(search_range) * 0.2

    # Convergence settings
    abstol           <- 1e-3
    stagnation_limit <- max(25L, max_iter %/% 5L)
    verbose          <- FALSE

    # Override with user-provided parameters
    if (!is.null(pso_params)) {
        if (!is.null(pso_params$w))                w                <- pso_params$w
        if (!is.null(pso_params$c1))               c1               <- pso_params$c1
        if (!is.null(pso_params$c2))               c2               <- pso_params$c2
        if (!is.null(pso_params$v_max))            v_max            <- pso_params$v_max
        if (!is.null(pso_params$abstol))           abstol           <- pso_params$abstol
        if (!is.null(pso_params$stagnation_limit)) stagnation_limit <- pso_params$stagnation_limit
        if (!is.null(pso_params$trace))            verbose          <- isTRUE(pso_params$trace)
    }

    # Memoized objective wrapper
    pso_obj <- create_fast_objective(objective_function)

    n_dim <- length(lower_bounds)

    # Initialise positions and velocities
    positions <- matrix(
        stats::runif(pop_size * n_dim,
                     rep(lower_bounds, each = pop_size),
                     rep(upper_bounds, each = pop_size)),
        pop_size, n_dim)
    velocities <- matrix(0, pop_size, n_dim)

    # Evaluate initial fitness (batch path skips per-particle fork overhead)
    fitness <- if (!is.null(batch_fn)) batch_fn(positions) else
                   eval_population_fitness(positions, pso_obj, parallel)

    pbest_pos <- positions
    pbest_fit <- fitness

    gbest_idx <- which.min(fitness)
    gbest_pos <- positions[gbest_idx, ]
    gbest_fit <- fitness[gbest_idx]

    stagnation_count <- 0L
    previous_best    <- Inf

    for (iter in seq_len(max_iter)) {

        # Generate flat random matrices (pop_size * n_dim each)
        r1_flat <- stats::runif(pop_size * n_dim)
        r2_flat <- stats::runif(pop_size * n_dim)

        # C++ velocity and position update
        upd <- pso_update_positions_cpp(
            positions, velocities, pbest_pos, gbest_pos,
            r1_flat, r2_flat,
            w, c1, c2,
            lower_bounds, upper_bounds, v_max)
        positions  <- upd$positions
        velocities <- upd$velocities

        # Evaluate fitness (batch or parallel scalar)
        fitness <- if (!is.null(batch_fn)) batch_fn(positions) else
                       eval_population_fitness(positions, pso_obj, parallel)

        # Update personal bests
        improved <- fitness < pbest_fit
        if (any(improved)) {
            pbest_pos[improved, ] <- positions[improved, , drop = FALSE]
            pbest_fit[improved]   <- fitness[improved]
        }

        # Update global best
        curr_best_idx <- which.min(pbest_fit)
        if (pbest_fit[curr_best_idx] < gbest_fit) {
            gbest_pos <- pbest_pos[curr_best_idx, ]
            gbest_fit <- pbest_fit[curr_best_idx]
        }

        if (verbose && (iter %% 20 == 0 || iter == 1L || iter == max_iter)) {
            cat(sprintf("PSO Iteration %d: Best fitness = %.8f\n", iter, gbest_fit))
        }

        # Early stopping
        improvement <- abs(previous_best - gbest_fit)
        if (improvement < abstol) {
            stagnation_count <- stagnation_count + 1L
        } else {
            stagnation_count <- 0L
        }
        if (stagnation_count >= stagnation_limit) {
            if (verbose) message(sprintf("PSO early stopping at iteration %d", iter))
            break
        }
        previous_best <- gbest_fit
    }

    data.frame(
        best_x     = gbest_pos[1],
        best_y     = gbest_pos[2],
        best_value = gbest_fit
    )
}


#' Create fast objective function with memoization
#' @noRd
create_fast_objective <- function(objective_function) {
    cache <- new.env(hash = TRUE)

    function(x, y) {
        key <- paste(round(c(x, y), 6), collapse = ",")

        if (exists(key, envir = cache)) {
            return(get(key, envir = cache))
        }

        value <- tryCatch({
            v <- objective_function(x, y)
            if (is.na(v) || !is.numeric(v) || is.infinite(v)) 1000 else v
        }, error = function(e) 1000)

        if (length(cache) < 5000) {
            assign(key, value, envir = cache)
        }

        return(value)
    }
}


#' Generate smart initial position based on search space
#' @noRd
generate_smart_initial_position <- function(lower_bounds, upper_bounds) {
    center <- (lower_bounds + upper_bounds) / 2
    noise  <- stats::runif(length(center), -0.1, 0.1) * (upper_bounds - lower_bounds)
    pmax(lower_bounds, pmin(upper_bounds, center + noise))
}
