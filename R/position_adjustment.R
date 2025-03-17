#' Positional adjustment of satellite LiDAR footprints
#'
#' This function adjusts the point locations of satellite LiDAR footprints with specified offsets.
#'
#' @param lidar_footprints Input satellite LiDAR footprints as an sf object
#' @param add_x X-axis offset
#' @param add_y Y-axis offset
#' @param crs_code Coordinate reference system code
#' @importFrom sf st_transform st_coordinates st_centroid st_geometry st_crs st_geometry<- st_as_sf
#' @return Adjusted satellite LiDAR footprints as an sf object with updated geometry
#' @export
position_adjustment <- function(lidar_footprints, add_x, add_y, crs_code){
    # Ensure input is an sf object with valid geometry
    if (!inherits(lidar_footprints, "sf")) {
        stop("Input lidar_footprints must be an sf object")
    }
    if (is.null(sf::st_geometry(lidar_footprints))) {
        stop("Input lidar_footprints does not have any geometry")
    }

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

    # Transform to specified CRS if needed
    if (sf::st_crs(lidar_footprints)$epsg != crs_code) {
        lidar_footprints <- st_transform(lidar_footprints, crs_code)
    }

    # Get original coordinates
    coords <- st_coordinates(st_centroid(lidar_footprints))

    # Create new coordinates with offsets
    new_coords <- data.frame(
        coords_x = coords[,1] + add_x,
        coords_y = coords[,2] + add_y
    )

    # Store original coordinates
    lidar_footprints$orig_x <- coords[,1]
    lidar_footprints$orig_y <- coords[,2]

    # Create new sf object with adjusted coordinates
    adjusted_footprints <- st_as_sf(
        cbind(st_drop_geometry(lidar_footprints), new_coords),
        coords = c("coords_x", "coords_y"),
        crs = crs_code
    )

    return(adjusted_footprints)
}
