#' Utility functions for the salpa package
#'
#' This file contains shared utility functions used across the package.

#' Validate CRS code
#'
#' @param crs_code Numeric EPSG code to validate
#' @return Validated CRS object
#' @noRd
validate_crs_code <- function(crs_code) {
    if (!is.numeric(crs_code)) {
        stop("crs_code must be a numeric EPSG code")
    }
    
    crs_obj_from_code <- try(sf::st_crs(crs_code), silent = TRUE)
    if (inherits(crs_obj_from_code, "try-error")) {
        stop(paste("Invalid crs_code:", crs_code, "- sf::st_crs() failed to interpret it."))
    }
    
    if (is.na(crs_obj_from_code$epsg)) {
        stop(paste("Invalid crs_code:", crs_code, "- Resulting CRS has no valid EPSG."))
    }
    
    return(crs_obj_from_code)
}

#' Validate buffer size
#'
#' @param buf Buffer size to validate
#' @return TRUE if valid, throws error if not
#' @noRd
validate_buffer <- function(buf) {
    if (buf <= 0) {
        stop("Buffer size must be a positive number")
    }
    return(TRUE)
}

#' Validate optimization method
#'
#' @param optimization_method Character string to validate
#' @return TRUE if valid, throws error if not
#' @noRd
validate_optimization_method <- function(optimization_method) {
    valid_methods <- c("ga", "pso", "woa", "lbfgsb")
    if (!(optimization_method %in% valid_methods)) {
        stop(paste("optimization_method must be one of:", paste(valid_methods, collapse = ", ")))
    }
    return(TRUE)
}

#' Convert bounds from meters to degrees for geographic coordinates
#'
#' @param bounds Numeric vector of bounds in meters
#' @return Numeric vector of bounds in degrees
#' @noRd
convert_bounds_to_degrees <- function(bounds) {
    # Approximate conversion: 1 degree ≈ 111,320 meters at equator
    return(bounds / 111320)
}

#' Create high-performance objective function for optimization
#'
#' @param lidar_footprints sf object containing satellite LiDAR footprints
#' @param input_rast Referenced raster object for extraction
#' @param buf Buffer size
#' @param minimizing_method Distance method to use
#' @param target_variable Target variable for extraction
#' @param crs_code CRS code
#' @param lidar_value Name of LiDAR value column
#' @return Function that takes x, y coordinates and returns loss value
#' @noRd
create_objective_function <- function(lidar_footprints, input_rast, buf, minimizing_method,
                                     target_variable, crs_code, lidar_value) {

    # Pre-extract LiDAR values once (these never change during optimization)
    lidar_val <- as.numeric(sf::st_drop_geometry(lidar_footprints[[lidar_value]]))

    # Wrap SpatRaster for fork-safe parallel evaluation (terra::wrap serialises it).
    # Each forked worker will call terra::unwrap() to get its own fresh SpatRaster,
    # avoiding the silent NULL-pointer crash that occurs when a terra object is
    # shared across fork() boundaries (macOS + exactextractr).
    wrapped_rast <- terra::wrap(input_rast)

    # PID-keyed cache so that:
    #   - serial runs:   unwrap once, reuse in every call (zero extra I/O)
    #   - parallel runs: each forked child gets its own fresh SpatRaster on first call
    rast_cache <- new.env(parent = emptyenv())
    rast_cache$rast <- NULL
    rast_cache$pid  <- NA_integer_

    function(x, y) {
        # Invalidate cache when we detect a new process (fork boundary)
        cur_pid <- Sys.getpid()
        if (is.null(rast_cache$rast) || rast_cache$pid != cur_pid) {
            rast_cache$rast <- terra::unwrap(wrapped_rast)
            rast_cache$pid  <- cur_pid
        }
        rast <- rast_cache$rast

        adjusted_geom <- sf::st_geometry(lidar_footprints) + c(x, y)
        sf::st_crs(adjusted_geom) <- sf::st_crs(lidar_footprints)
        temp_footprints <- sf::st_set_geometry(lidar_footprints, adjusted_geom)
        buffered_data <- sf::st_buffer(temp_footprints, buf)

        ref_val <- exactextractr::exact_extract(rast, buffered_data,
                                              target_variable, progress = FALSE)
        ref_val <- as.numeric(ref_val)

        valid_indices <- !is.na(lidar_val) & !is.na(ref_val) &
                         is.finite(lidar_val) & is.finite(ref_val)

        if (sum(valid_indices) < 1) return(1000)

        distance_result <- perform_distance(lidar_val[valid_indices],
                                          ref_val[valid_indices],
                                          method = minimizing_method)

        if (is.na(distance_result$distance) || !distance_result$valid) return(1000)

        return(distance_result$distance)
    }
}

#' Handle optimization failure
#'
#' @param error Error object from failed optimization
#' @return Data frame with default values
#' @noRd
handle_optimization_failure <- function(error) {
    warning("Optimization failed: ", error$message)
    return(data.frame(best_x = 0, best_y = 0, best_value = NA))
}

#' Validate and process optimization parameters
#'
#' @param pop_size Population size
#' @param max_iter Maximum iterations  
#' @param lower_bounds Lower bounds vector
#' @param upper_bounds Upper bounds vector
#' @return List of validated parameters
#' @noRd
validate_optimization_params <- function(pop_size, max_iter, lower_bounds, upper_bounds) {
    if (pop_size <= 0) {
        stop("Population size must be positive")
    }
    if (max_iter <= 0) {
        stop("Maximum iterations must be positive")
    }
    if (length(lower_bounds) != length(upper_bounds)) {
        stop("Lower and upper bounds must have same length")
    }
    if (any(lower_bounds >= upper_bounds)) {
        stop("Lower bounds must be less than upper bounds")
    }
    
    return(list(
        pop_size = as.integer(pop_size),
        max_iter = as.integer(max_iter),
        lower_bounds = as.numeric(lower_bounds),
        upper_bounds = as.numeric(upper_bounds)
    ))
}

#' Create batch objective function for population-based optimizers
#'
#' Returns a function(positions_matrix) that evaluates ALL particles in a single
#' exactextractr call, eliminating per-particle setup overhead.
#' Uses per-position memoization to skip re-evaluation of already-seen offsets.
#'
#' @param lidar_footprints sf object of LiDAR footprints
#' @param input_rast SpatRaster reference DEM
#' @param buf Buffer radius in metres
#' @param minimizing_method Distance method (euclidean, manhattan, etc.)
#' @param target_variable exactextractr summary stat (e.g. "mean")
#' @param crs_code EPSG code (unused directly; CRS taken from lidar_footprints)
#' @param lidar_value Column name of LiDAR measurement
#' @return Function(positions_matrix) → numeric vector of length pop_size
#' @noRd
create_batch_objective_function <- function(lidar_footprints, input_rast, buf,
                                            minimizing_method, target_variable,
                                            crs_code, lidar_value) {

    lidar_val    <- as.numeric(sf::st_drop_geometry(lidar_footprints[[lidar_value]]))
    wrapped_rast <- terra::wrap(input_rast)
    n_shots      <- nrow(lidar_footprints)
    base_geom    <- sf::st_geometry(lidar_footprints)
    fp_crs       <- sf::st_crs(lidar_footprints)

    rast_cache <- new.env(parent = emptyenv())
    rast_cache$rast <- NULL
    rast_cache$pid  <- NA_integer_

    memo_cache <- new.env(hash = TRUE, parent = emptyenv())

    function(positions_matrix) {
        pop_size <- nrow(positions_matrix)

        # Per-position memoization check
        keys   <- apply(positions_matrix, 1, function(r) paste(round(r, 6L), collapse = ","))
        cached <- vapply(keys, exists, logical(1L), envir = memo_cache, USE.NAMES = FALSE)

        loss_vals <- rep(NA_real_, pop_size)
        cached_idx <- which(cached)
        if (length(cached_idx) > 0L) {
            loss_vals[cached_idx] <- vapply(keys[cached_idx], get,
                                           numeric(1L), envir = memo_cache,
                                           USE.NAMES = FALSE)
        }

        new_idx <- which(!cached)
        if (length(new_idx) == 0L) return(loss_vals)

        # PID-keyed raster cache (fork-safe)
        cur_pid <- Sys.getpid()
        if (is.null(rast_cache$rast) || rast_cache$pid != cur_pid) {
            rast_cache$rast <- terra::unwrap(wrapped_rast)
            rast_cache$pid  <- cur_pid
        }
        rast    <- rast_cache$rast
        new_pos <- positions_matrix[new_idx, , drop = FALSE]
        n_new   <- nrow(new_pos)

        # Build ALL shifted+buffered geometries for new positions at once.
        # lapply avoids repeated st_crs assignment inside a loop; CRS set once at end.
        all_geom_list <- lapply(seq_len(n_new), function(ii) {
            sf::st_buffer(base_geom + c(new_pos[ii, 1L], new_pos[ii, 2L]), buf)
        })
        all_geoms <- do.call(c, all_geom_list)
        sf::st_crs(all_geoms) <- fp_crs   # set CRS once for the combined sfc

        # Single exactextractr call for all n_new * n_shots geometries
        all_ref <- as.numeric(exactextractr::exact_extract(
            rast, all_geoms, target_variable, progress = FALSE))

        # Compute loss per new position and populate cache
        for (ii in seq_len(n_new)) {
            i       <- new_idx[ii]
            start_i <- (ii - 1L) * n_shots + 1L
            end_i   <- ii * n_shots
            ref_val <- all_ref[start_i:end_i]

            valid <- !is.na(lidar_val) & !is.na(ref_val) &
                     is.finite(lidar_val) & is.finite(ref_val)

            lv <- if (sum(valid) < 1L) {
                1000
            } else {
                dr <- perform_distance(lidar_val[valid], ref_val[valid],
                                       method = minimizing_method)
                if (is.na(dr$distance) || !dr$valid) 1000 else dr$distance
            }

            loss_vals[i] <- lv
            if (length(memo_cache) < 10000L) assign(keys[i], lv, envir = memo_cache)
        }

        loss_vals
    }
}


#' Evaluate population fitness, optionally in parallel
#'
#' Shared helper used by WOA, PSO, and GA to evaluate all individuals.
#' On macOS/Linux uses \code{parallel::mclapply}; on Windows falls back to
#' serial evaluation (fork-based parallelism is unavailable there).
#'
#' @param positions Numeric matrix, pop_size x n_dim
#' @param objective_fn Function taking (x, y) and returning a numeric scalar
#' @param parallel Logical; if TRUE attempt parallel evaluation
#' @return Numeric vector of length pop_size
#' @noRd
eval_population_fitness <- function(positions, objective_fn, parallel = FALSE) {
    pop_size <- nrow(positions)
    if (parallel && requireNamespace("parallel", quietly = TRUE) &&
        .Platform$OS.type != "windows") {
        n_cores <- if (is.numeric(parallel) && parallel > 1L) {
            as.integer(parallel)
        } else {
            max(1L, parallel::detectCores() - 1L)
        }
        unlist(parallel::mclapply(seq_len(pop_size),
            function(i) objective_fn(positions[i, 1], positions[i, 2]),
            mc.cores = n_cores))
    } else {
        vapply(seq_len(pop_size),
            function(i) objective_fn(positions[i, 1], positions[i, 2]),
            numeric(1))
    }
}