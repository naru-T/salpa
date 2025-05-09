#' Linear Alignment of Points
#'
#' This function aligns satellite LiDAR footprints to a best-fit line using an analytical total least squares approach.
#' The method computes the principal component of the coordinates for each shot group to minimize the perpendicular distances
#' between the points and the line, and then projects the points onto this line.
#'
#' @param sf_object An sf object containing spatial features. It must include a column named 'shot_number'.
#' @param crs_code An integer representing the coordinate reference system (CRS) code to which the sf object will be transformed.
#'
#' @return An sf object with adjusted coordinates that best fit a line, including original coordinates stored in 'orig_x' and 'orig_y'.
#'
#' @importFrom sf st_transform st_coordinates st_as_sf
#' @importFrom dplyr bind_rows
#' @importFrom stats cov
#'
#' @examples
#' \dontrun{
#' library(sf)
#' 
#' # Example usage
#' sf_object <- st_read('path_to_your_shapefile.shp')
#' crs_code <- 3857
#' adjusted_sf <- linear_alignment(sf_object, crs_code)
#' }
#'
#' @export
linear_alignment <- function(sf_object, crs_code) {
    # Validate CRS code
    if (!is.numeric(crs_code)) {
        stop("crs_code must be a numeric EPSG code")
    }
    
    # Check if the CRS exists
    tryCatch({
        crs <- sf::st_crs(crs_code)
        if (is.na(crs)) {
            stop("Invalid crs_code: ", crs_code, ". Please provide a valid EPSG code.")
        }
    }, error = function(e) {
        stop("Invalid crs_code: ", crs_code, ". Please provide a valid EPSG code.")
    })
    
    # Check for shot_number column
    if (!"shot_number" %in% colnames(sf_object)) {
        stop("The sf_object does not contain a 'shot_number' column, which is required for the function to work correctly.")
    }
    
    # Transform to specified CRS
    sf_object <- st_transform(sf_object, crs_code)
    
    # Extract unique shot groups (using the first 10 characters of shot_number)
    shot_ids <- unique(substr(as.character(sf_object$shot_number), 1, 10))
    
    adjusted_sf_list <- list()
    
    for (shot in shot_ids) {
        # Subset the sf object for the current shot group
        sf_subset <- sf_object[substr(as.character(sf_object$shot_number), 1, 10) == shot, ]
        
        # Extract coordinates
        coords <- st_coordinates(sf_subset)
        
        # If there's only one point, use the point as is
        if (nrow(coords) == 1) {
            sf_subset$orig_x <- coords[1]
            sf_subset$orig_y <- coords[2]
            sf_subset$coords_x <- coords[1]
            sf_subset$coords_y <- coords[2]
        } else {
            # Compute the centroid
            centroid <- colMeans(coords)
            
            # Center the coordinates
            centered <- sweep(coords, 2, centroid)
            
            # Compute the covariance matrix and perform eigen decomposition
            cov_mat <- cov(centered)
            eig <- eigen(cov_mat)
            
            # Get the principal eigenvector (direction of maximum variance)
            direction <- eig$vectors[,1]
            
            # Project each point onto the line passing through the centroid in the direction of 'direction'
            projections <- t(apply(coords, 1, function(pt) {
                proj <- centroid + sum((pt - centroid) * direction) * direction
                return(proj)
            }))
            
            # Store original and projected coordinates
            sf_subset$orig_x <- coords[,1]
            sf_subset$orig_y <- coords[,2]
            sf_subset$coords_x <- projections[,1]
            sf_subset$coords_y <- projections[,2]
        }
        
        # Create an adjusted sf object with the new projected coordinates
        adjusted_sf <- st_as_sf(st_drop_geometry(sf_subset), coords = c("coords_x", "coords_y"), crs = crs_code)
        
        # Append to the list
        adjusted_sf_list[[length(adjusted_sf_list) + 1]] <- adjusted_sf
    }
    
    # Merge all adjusted sf objects into one
    final_adjusted_sf <- bind_rows(adjusted_sf_list)
    
    return(final_adjusted_sf)
}
