#' Positional adjustment of satellite LiDAR footprints
#'
#' This function adjusts the point locations of satellite LiDAR footprints with specified offsets.
#'
#' @param lidar_footprints Input satellite LiDAR footprints as an sf object
#' @param add_x X-axis offset
#' @param add_y Y-axis offset
#' @param buf Buffer size
#' @param crs_code Coordinate reference system code
#' @importFrom sf  st_transform
#' @importFrom sf  st_coordinates
#' @importFrom sf  st_centroid
#' @importFrom sf  st_as_sf
#' @importFrom sf  st_buffer
#' @return Adjusted satellite LiDAR footprints
#' @export
positional_correction <- function(lidar_footprints, add_x, add_y, buf, crs_code){
    lidar_footprints <- st_transform(lidar_footprints, crs_code)
    coords <- st_coordinates(st_centroid(lidar_footprints))

    coords_x <- coords[,1] + add_x
    coords_y <- coords[,2] + add_y

    lidar_footprints$orig_x <- coords[,1]
    lidar_footprints$orig_y <- coords[,2]
    lidar_footprints$coords_x <- coords_x
    lidar_footprints$coords_y <- coords_y

    sf_newdata <- st_as_sf(st_drop_geometry(lidar_footprints), coords = c("coords_x", "coords_y"), crs = crs_code)
    sf_newdata <- st_buffer(sf_newdata, buf)
    return(sf_newdata)
}