#' Optimization factory function
#'
#' Factory function to select and run the appropriate optimization method.
#'
#' @param optimization_method Character string for optimization method
#' @param objective_function Function to minimize
#' @param lower_bounds Lower bounds for optimization
#' @param upper_bounds Upper bounds for optimization
#' @param pop_size Population size
#' @param max_iter Maximum iterations
#' @param pso_params Additional PSO parameters
#' @param woa_params Additional WOA parameters
#' @param lbfgsb_params Additional L-BFGS-B parameters
#' @param parallel Logical indicating whether to use parallel processing
#' @return Data frame with optimization results
#' @noRd
optimization_factory <- function(optimization_method, objective_function, lower_bounds,
                                upper_bounds, pop_size, max_iter,
                                pso_params = NULL, woa_params = NULL, lbfgsb_params = NULL,
                                parallel = FALSE, batch_fn = NULL) {

    # Validate optimization method
    validate_optimization_method(optimization_method)

    # Validate optimization parameters
    params <- validate_optimization_params(pop_size, max_iter, lower_bounds, upper_bounds)

    # Dispatch to appropriate optimization function
    switch(optimization_method,
        "ga" = optimize_ga(
            objective_function = objective_function,
            lower_bounds = params$lower_bounds,
            upper_bounds = params$upper_bounds,
            pop_size = params$pop_size,
            max_iter = params$max_iter,
            parallel = parallel,
            batch_fn = batch_fn
        ),
        "pso" = optimize_pso(
            objective_function = objective_function,
            lower_bounds = params$lower_bounds,
            upper_bounds = params$upper_bounds,
            pop_size = params$pop_size,
            max_iter = params$max_iter,
            pso_params = pso_params,
            parallel = parallel,
            batch_fn = batch_fn
        ),
        "woa" = optimize_woa(
            objective_function = objective_function,
            lower_bounds = params$lower_bounds,
            upper_bounds = params$upper_bounds,
            pop_size = params$pop_size,
            max_iter = params$max_iter,
            woa_params = woa_params,
            parallel = parallel,
            batch_fn = batch_fn
        ),
        "lbfgsb" = optimize_lbfgsb(
            objective_function = objective_function,
            lower_bounds = params$lower_bounds,
            upper_bounds = params$upper_bounds,
            max_iter = params$max_iter,
            lbfgsb_params = lbfgsb_params,
            parallel = parallel
            # lbfgsb uses gradient descent: no batch_fn needed
        ),
        stop("Invalid optimization_method. This should not happen as validation should catch this.")
    )
}