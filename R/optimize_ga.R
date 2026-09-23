#' High-Performance Genetic Algorithm Optimization (custom implementation)
#'
#' Uses an explicit GA loop with C++ population evolution (ga_evolve_population_cpp)
#' and fork-safe parallel fitness evaluation via eval_population_fitness().
#' This replaces the GA package dependency and enables parallel=TRUE safely.
#'
#' @param objective_function Function to minimize
#' @param lower_bounds Lower bounds for optimization
#' @param upper_bounds Upper bounds for optimization
#' @param pop_size Population size
#' @param max_iter Maximum iterations
#' @param ga_params Additional GA parameters (pcrossover, pmutation, run, elitism, tournament_k)
#' @param parallel Logical indicating whether to use parallel fitness evaluation
#' @return Data frame with optimization results
#' @noRd
optimize_ga <- function(objective_function, lower_bounds, upper_bounds,
                        pop_size, max_iter, ga_params = NULL, parallel = FALSE,
                        batch_fn = NULL) {

    return(run_custom_ga(objective_function, lower_bounds, upper_bounds,
                         pop_size, max_iter, ga_params, parallel,
                         batch_fn = batch_fn))
}


#' Custom GA loop: C++ population evolution + fork-safe parallel fitness eval
#' @noRd
run_custom_ga <- function(objective_function, lower_bounds, upper_bounds,
                          pop_size, max_iter, ga_params = NULL, parallel = FALSE,
                          batch_fn = NULL) {

    # --- Hyperparameters (overridable via ga_params) ---
    pcrossover   <- 0.8
    pmutation    <- 0.1
    elitism_k    <- max(1L, pop_size %/% 10L)
    tournament_k <- max(2L, pop_size %/% 10L)
    early_stop   <- max(25L, max_iter %/% 8L)

    if (!is.null(ga_params)) {
        if (!is.null(ga_params$pcrossover))   pcrossover   <- ga_params$pcrossover
        if (!is.null(ga_params$pmutation))    pmutation    <- ga_params$pmutation
        if (!is.null(ga_params$run))          early_stop   <- ga_params$run
        if (!is.null(ga_params$elitism))      elitism_k    <- as.integer(ga_params$elitism)
        if (!is.null(ga_params$tournament_k)) tournament_k <- as.integer(ga_params$tournament_k)
    }

    cat("GA bounds:", lower_bounds, "to", upper_bounds, "\n")
    cat("GA pop_size:", pop_size, "max_iter:", max_iter,
        "parallel:", parallel, "\n")

    # --- Memoized fitness wrapper (GA maximises, so negate objective) ---
    # create_memoized_fitness returns function(x_vec) where x_vec = c(x, y).
    # eval_population_fitness calls function(x, y) with two scalars, so we
    # wrap it to match the two-argument interface.
    memoized_fn    <- create_memoized_fitness(objective_function)
    cached_fitness <- function(x, y) memoized_fn(c(x, y))

    # --- Initialise population with smart suggestions + random fill ---
    suggestions <- generate_smart_initial_population(pop_size, lower_bounds, upper_bounds)
    n_sugg <- if (!is.null(suggestions)) nrow(suggestions) else 0L
    n_rand <- pop_size - n_sugg

    rand_pop <- matrix(
        c(stats::runif(n_rand, lower_bounds[1], upper_bounds[1]),
          stats::runif(n_rand, lower_bounds[2], upper_bounds[2])),
        nrow = n_rand, ncol = 2L)

    pop <- if (n_sugg > 0L) rbind(suggestions, rand_pop) else rand_pop

    # --- Initial fitness evaluation ---
    # batch_fn returns loss (minimize); GA maximises, so negate.
    fitness <- if (!is.null(batch_fn)) -batch_fn(pop) else
                   eval_population_fitness(pop, cached_fitness, parallel)

    best_fitness   <- max(fitness)
    no_improve_cnt <- 0L

    # --- Evolution loop ---
    for (iter in seq_len(max_iter)) {

        # C++ one-generation evolution (selection + crossover + mutation + elitism)
        pop <- ga_evolve_population_cpp(
            population   = pop,
            fitness      = fitness,
            lower_bounds = lower_bounds,
            upper_bounds = upper_bounds,
            pcrossover   = pcrossover,
            pmutation    = pmutation,
            elitism_k    = elitism_k,
            tournament_k = tournament_k
        )

        fitness <- if (!is.null(batch_fn)) -batch_fn(pop) else
                       eval_population_fitness(pop, cached_fitness, parallel)

        cur_best <- max(fitness)
        if (cur_best > best_fitness + 1e-6) {
            best_fitness   <- cur_best
            no_improve_cnt <- 0L
        } else {
            no_improve_cnt <- no_improve_cnt + 1L
            if (no_improve_cnt >= early_stop) break
        }
    }

    best_idx <- which.max(fitness)
    data.frame(
        best_x     = pop[best_idx, 1L],
        best_y     = pop[best_idx, 2L],
        best_value = -fitness[best_idx]   # convert back to minimisation value
    )
}


#' Create memoized fitness function for GA (negates objective for maximisation)
#' @noRd
create_memoized_fitness <- function(objective_function) {
    cache      <- new.env(hash = TRUE)
    eval_count <- 0L

    function(x) {
        key <- paste(round(x, 6), collapse = ",")

        if (exists(key, envir = cache)) {
            return(get(key, envir = cache))
        }

        fitness_value <- tryCatch({
            eval_count <<- eval_count + 1L
            if (eval_count %% 100L == 1L) {
                cat("[GA fitness eval #", eval_count, "] x=(",
                    x[1], ",", x[2], ") ", sep = "")
            }
            obj_value <- objective_function(x[1], x[2])
            if (eval_count %% 100L == 1L) {
                cat("-> obj_value =", obj_value, "\n")
            }
            if (is.na(obj_value) || !is.numeric(obj_value) || is.infinite(obj_value)) {
                cat("[GA fitness invalid] x=(", x[1], ",", x[2],
                    ") obj_value =", obj_value, "\n")
                -1000
            } else {
                -obj_value   # GA maximises, so negate
            }
        }, error = function(e) {
            cat("[GA fitness error] x=(", x[1], ",", x[2], ") error:", e$message, "\n")
            -1000
        })

        if (length(cache) < 1500L) {
            assign(key, fitness_value, envir = cache)
        }

        return(fitness_value)
    }
}


#' Generate smart initial population for GA
#' @noRd
generate_smart_initial_population <- function(pop_size, lower_bounds, upper_bounds) {
    n_suggestions <- min(5L, pop_size %/% 4L)
    if (n_suggestions < 1L) return(NULL)

    suggestions <- matrix(0, nrow = n_suggestions, ncol = length(lower_bounds))

    center <- (lower_bounds + upper_bounds) / 2
    suggestions[1, ] <- center

    if (n_suggestions > 1L) {
        if (length(lower_bounds) == 2L && n_suggestions >= 2L) {
            corners <- expand.grid(c(lower_bounds[1], upper_bounds[1]),
                                   c(lower_bounds[2], upper_bounds[2]))
            for (i in 2L:min(n_suggestions, nrow(corners) + 1L)) {
                suggestions[i, ] <- as.numeric(corners[i - 1L, ])
            }
        }

        remaining <- n_suggestions - min(n_suggestions, 5L)
        if (remaining > 0L) {
            for (i in (n_suggestions - remaining + 1L):n_suggestions) {
                for (j in seq_along(lower_bounds)) {
                    suggestions[i, j] <- stats::runif(1, lower_bounds[j], upper_bounds[j])
                }
            }
        }
    }

    return(suggestions)
}
