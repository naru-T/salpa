#' Positional adjustment of satellite LiDAR footprints
#'
#' This function adjusts the point locations of satellite LiDAR footprints with specified offsets.
#'
#' @param lidar_footprints Input satellite LiDAR footprints as an sf object
#' @param add_x X-axis offset
#' @param add_y Y-axis offset
#' @param crs_code Coordinate reference system code
#' @importFrom sf  st_transform
#' @importFrom sf  st_coordinates
#' @importFrom sf  st_centroid
#' @importFrom sf  st_as_sf
#' @return Adjusted satellite LiDAR footprints
#' @export
position_adjustment <- function(lidar_footprints, add_x, add_y, crs_code){    
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
    
    # Transform to specified CRS
    lidar_footprints <- st_transform(lidar_footprints, crs_code)
    coords <- st_coordinates(st_centroid(lidar_footprints))

    coords_x <- coords[,1] + add_x
    coords_y <- coords[,2] + add_y

    lidar_footprints$orig_x <- coords[,1]
    lidar_footprints$orig_y <- coords[,2]
    lidar_footprints$coords_x <- coords_x
    lidar_footprints$coords_y <- coords_y

    sf_newdata <- st_as_sf(st_drop_geometry(lidar_footprints), coords = c("coords_x", "coords_y"), crs = crs_code)
    return(sf_newdata)
}