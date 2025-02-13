#' Linear Alignment of Points
#'
#' This function aligns satellite LiDAR footprints to a best fit line using genetic algorithm optimization.
#' The sf object should be transformed to the specified (projected) CRS before using this function.
#' The fitting is done by minimizing the perpendicular distances between the points and the line using a genetic algorithm (GA) optimization.
#'
#' @param sf_object An sf object containing spatial features.
#' @param crs_code An integer representing the coordinate reference system (CRS) code to transform the sf object.
#'
#' @return An sf object with adjusted coordinates that best fit a line.
#'
#' @importFrom sf st_transform st_coordinates st_as_sf st_crs
#' @importFrom GA ga
#' @importFrom dplyr bind_rows
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(GA)
#'
#' # Example usage
#' sf_object <- st_read("path_to_your_shapefile.shp")
#' crs_code <- 3857
#' adjusted_sf <- linear_alignment(sf_object, crs_code)
#' }
#'
#' @export
linear_alignment <- function(sf_object, crs_code) {
    sf_object <- st_transform(sf_object, crs_code)
    if (!"shot_number" %in% colnames(sf_object)) {
        stop("The sf_object does not contain a 'shot_number' column, which is required for the function to work correctly.")
    }
    shot_number <- as.character(sf_object$shot_number)
    shot_number <- substr(shot_number, 1, 12)
    # remove duplicated shot number
    shot_number <- unique(shot_number)
    # Initialize an empty list to store adjusted sf objects
    adjusted_sf_list <- list()
    # Function to calculate total error between line and points
    calculate_error <- function(params, x_points, y_points) {
        slope <- params[1]
        offset <- params[2]

        # Calculate perpendicular distances (total error)
        distances <- abs(-slope * x_points + y_points - offset) /
                     sqrt(slope^2 + 1)

        # Return sum of squared distances
        return(sum(distances^2))
    }

    # Initial line parameters from two points
    get_initial_params <- function(x, y) {
        # Sort points
        sorted_idx <- order(x)
        x_sorted <- x[sorted_idx]
        y_sorted <- y[sorted_idx]

        # Get first and last points
        x1 <- x_sorted[1]
        y1 <- y_sorted[1]
        x2 <- x_sorted[length(x_sorted)]
        y2 <- y_sorted[length(y_sorted)]

        # Calculate initial slope and offset
        slope_0 <- (y2 - y1) / (x2 - x1)
        offset_0 <- y1 - slope_0 * x1

        return(c(slope_0, offset_0))
    }

    # Optimization using GA
    optimize_line <- function(x, y, slope_0, offset_0) {
        ga_result <- ga(
            type = "real-valued",
            fitness = function(params) -calculate_error(params, x, y), # Negative because GA maximizes
            lower = c(slope_0 - 1, offset_0 - 10), # Adjust ranges as needed
            upper = c(slope_0 + 1, offset_0 + 10),
            popSize = 50,
            maxiter = 100,
            monitor = FALSE
        )

        return(ga_result@solution)
    }

    # Loop through each unique shot number
    for (shot in shot_number) {
        # Subset the sf object for the current shot number
        sf_object_subset <- sf_object[substr(sf_object$shot_number, 1, 12) == shot, ]
        # Extract coordinates from the subset sf object
        coords <- st_coordinates(sf_object_subset)
        coords_x <- coords[, 1]
        coords_y <- coords[, 2]

        # Get initial parameters
        initial_params <- get_initial_params(coords_x, coords_y)
        slope_0 <- initial_params[1]
        offset_0 <- initial_params[2]

        # Optimize parameters
        optimized_params <- optimize_line(coords_x, coords_y, slope_0, offset_0)

        slope <- optimized_params[1]
        intercept <- optimized_params[2]

        # Project points onto line using perpendicular projection formula
        x_proj <- (coords_x + slope * coords_y - slope * intercept) / (1 + slope^2)
        y_proj <- intercept + slope * x_proj

        # Create adjusted sf object
        sf_object_subset$orig_x <- coords[,1]
        sf_object_subset$orig_y <- coords[,2]
        sf_object_subset$coords_x <- x_proj
        sf_object_subset$coords_y <- y_proj
        adjusted_sf <- st_as_sf(st_drop_geometry(sf_object_subset), coords = c("coords_x", "coords_y"), crs = crs_code)

        # Append to the list
        adjusted_sf_list[[length(adjusted_sf_list) + 1]] <- adjusted_sf
    }

    # Merge all adjusted sf objects into one
    final_adjusted_sf <- dplyr::bind_rows(adjusted_sf_list)

    return(final_adjusted_sf)
}
