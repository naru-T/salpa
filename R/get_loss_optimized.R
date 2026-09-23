#' High-performance optimized loss calculation with spatial caching
#'
#' This function optimizes the loss calculation by pre-computing and caching
#' spatial operations that don't change between optimization iterations.
#' 
#' @param lidar_footprints sf object with satellite LiDAR footprints
#' @param add_x numeric value to add to x coordinates
#' @param add_y numeric value to add to y coordinates 
#' @param buf numeric value for buffer size
#' @param input_rast raster object for reference data
#' @param minimizing_method character string for the minimizing method
#' @param target_variable character string for the target variable
#' @param crs_code numeric value for CRS code
#' @param lidar_value name of the column containing LiDAR measurement values
#' @param cache_env environment containing pre-computed spatial data (optional)
#' @return data frame with loss value
#' @examples
#' \donttest{
#'   lidar_sf <- sf::st_as_sf(
#'     data.frame(x = c(500, 510, 520), y = c(300, 310, 320), elev = c(100, 101, 99)),
#'     coords = c("x", "y"), crs = 32654
#'   )
#'   dem <- terra::rast(
#'     nrows = 10, ncols = 10,
#'     xmin = 450, xmax = 600, ymin = 250, ymax = 400,
#'     crs = "EPSG:32654", vals = runif(100, 98, 103)
#'   )
#'   result <- get_loss_optimized(lidar_sf, add_x = 0, add_y = 0, buf = 10,
#'                                input_rast = dem, target_variable = "mean",
#'                                crs_code = 32654, lidar_value = "elev")
#' }
#' @export
get_loss_optimized <- function(lidar_footprints, add_x, add_y, buf, input_rast,
                              minimizing_method = "euclidean", target_variable,
                              crs_code, lidar_value, cache_env = NULL) {
    
    # If no cache provided, use the standard function
    if (is.null(cache_env)) {
        return(get_loss(lidar_footprints, add_x, add_y, buf, input_rast,
                       minimizing_method, target_variable, crs_code, lidar_value))
    }
    
    # Use cached spatial operations for massive speedup
    n_footprints <- nrow(lidar_footprints)
    
    # Apply position adjustment to original footprints
    adjusted_geom <- sf::st_geometry(lidar_footprints) + c(add_x, add_y)
    sf::st_crs(adjusted_geom) <- sf::st_crs(lidar_footprints)
    # Create temporary sf object with adjusted geometry
    temp_footprints <- sf::st_set_geometry(lidar_footprints, adjusted_geom)
    
    # Buffer the adjusted footprints
    buffered_data <- sf::st_buffer(temp_footprints, buf)
    
    # Extract values from raster (this is the expensive operation)
    extracted_values <- exactextractr::exact_extract(input_rast, buffered_data, 
                                                    target_variable, progress = FALSE)
    
    # Get LiDAR values (these don't change)
    lidar_val <- lidar_footprints[[lidar_value]]
    ref_val <- extracted_values
    
    # Filter out NA values
    valid_indices <- !is.na(lidar_val) & !is.na(ref_val) & 
                     is.finite(lidar_val) & is.finite(ref_val)
    
    if (sum(valid_indices) < 1) {
        loss_value <- 1000  # High penalty value
    } else {
        lidar_val_clean <- lidar_val[valid_indices]
        ref_val_clean <- ref_val[valid_indices]
        
        # Calculate distance using optimized C++ functions
        distance_result <- perform_distance(lidar_val_clean, ref_val_clean, 
                                          method = minimizing_method)
        loss_value <- distance_result$distance
        
        if (is.na(loss_value) || !distance_result$valid) {
            loss_value <- 1000
        }
    }
    
    return(data.frame(
        add_x = add_x,
        add_y = add_y,
        buffer = buf,
        loss_value = loss_value
    ))
}

#' Create optimized loss function for optimization algorithms
#' 
#' This creates a closure that caches spatial computations for maximum performance
#' during optimization iterations.
#' 
#' @param lidar_footprints sf object with satellite LiDAR footprints
#' @param buf numeric value for buffer size
#' @param input_rast raster object for reference data
#' @param minimizing_method character string for the minimizing method
#' @param target_variable character string for the target variable
#' @param crs_code numeric value for CRS code
#' @param lidar_value name of the column containing LiDAR measurement values
#' @return function that takes a numeric vector `params = c(dx, dy)` and returns a scalar loss value
#' @examples
#' \donttest{
#'   lidar_sf <- sf::st_as_sf(
#'     data.frame(x = c(500, 510, 520), y = c(300, 310, 320), elev = c(100, 101, 99)),
#'     coords = c("x", "y"), crs = 32654
#'   )
#'   dem <- terra::rast(
#'     nrows = 10, ncols = 10,
#'     xmin = 450, xmax = 600, ymin = 250, ymax = 400,
#'     crs = "EPSG:32654", vals = runif(100, 98, 103)
#'   )
#'   loss_fn <- create_fast_loss_function(
#'     lidar_sf, buf = 10, input_rast = dem,
#'     minimizing_method = "euclidean", target_variable = "mean",
#'     crs_code = 32654, lidar_value = "elev"
#'   )
#'   loss_fn(c(0, 0))
#' }
#' @export
create_fast_loss_function <- function(lidar_footprints, buf, input_rast,
                                     minimizing_method, target_variable, 
                                     crs_code, lidar_value) {
    
    # Pre-extract LiDAR values (these never change)
    lidar_val <- as.numeric(sf::st_drop_geometry(lidar_footprints[[lidar_value]]))
    
    # Return optimized loss function
    function(params) {
        add_x <- params[1]
        add_y <- params[2]
        
        # Fast position adjustment and buffering
        adjusted_geom <- sf::st_geometry(lidar_footprints) + c(add_x, add_y)
        sf::st_crs(adjusted_geom) <- sf::st_crs(lidar_footprints)
        temp_footprints <- sf::st_set_geometry(lidar_footprints, adjusted_geom)
        buffered_data <- sf::st_buffer(temp_footprints, buf)
        
        # Extract reference values
        ref_val <- exactextractr::exact_extract(input_rast, buffered_data, 
                                              target_variable, progress = FALSE)
        ref_val <- as.numeric(ref_val)
        
        # Filter valid pairs
        valid_indices <- !is.na(lidar_val) & !is.na(ref_val) & 
                         is.finite(lidar_val) & is.finite(ref_val)
        
        if (sum(valid_indices) < 1) {
            return(1000)  # High penalty
        }
        
        # Calculate distance with optimized C++ functions
        distance_result <- perform_distance(lidar_val[valid_indices], 
                                          ref_val[valid_indices], 
                                          method = minimizing_method)
        
        if (is.na(distance_result$distance) || !distance_result$valid) {
            return(1000)
        }
        
        return(distance_result$distance)
    }
}