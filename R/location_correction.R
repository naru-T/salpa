#' Correct the location of the data
#'
#' @param gedi_l2a sf object of gedi l2a
#' @param add_x numeric value to add to x coordinates
#' @param add_y numeric value to add to y coordinates
#' @param buf numeric value for buffer size
#' @param crs_code numeric value for CRS code
#' @importFrom sf  st_transform
#' @importFrom sf  st_coordinates
#' @importFrom sf  st_centroid
#' @importFrom sf  st_as_sf
#' @importFrom sf  st_buffer
#' @return sf object with corrected location
#' @export
location_correction <- function(gedi_l2a, add_x, add_y, buf, crs_code){
    gedi_l2a <- st_transform(gedi_l2a, crs_code)
    coords <- st_coordinates(st_centroid(gedi_l2a))

    coords_x <- coords[,1] + add_x
    coords_y <- coords[,2] + add_y

    gedi_l2a$orig_x <- coords[,1]
    gedi_l2a$orig_y <- coords[,2]
    gedi_l2a$coords_x <- coords_x
    gedi_l2a$coords_y <- coords_y

    sf_newdata <- st_as_sf(gedi_l2a |> st_drop_geometry(), coords = c("coords_x", "coords_y"), crs = crs_code)
    sf_newdata <- st_buffer(sf_newdata, buf)
    return(sf_newdata)
}