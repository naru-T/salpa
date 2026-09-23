#' Extract values of input raster data in the buffered data
#'
#' @param buffered_data sf object with buffered data
#' @param input_rast raster object for extraction
#' @param target_variable character string for the target variable (min, max, mean, median, etc.). See exact_extract
#' @return sf object with an added `ref_val` column containing the extracted raster summary
#' @importFrom exactextractr exact_extract
#' @importFrom dplyr mutate
#' @examples
#' \donttest{
#'   pts <- sf::st_as_sf(
#'     data.frame(x = c(500, 510), y = c(300, 310)),
#'     coords = c("x", "y"), crs = 32654
#'   )
#'   buffered <- sf::st_buffer(pts, dist = 10)
#'   dem <- terra::rast(
#'     nrows = 10, ncols = 10,
#'     xmin = 450, xmax = 600, ymin = 250, ymax = 400,
#'     crs = "EPSG:32654", vals = runif(100, 98, 103)
#'   )
#'   result <- extract_values(buffered, dem, target_variable = "mean")
#' }
#' @export
extract_values <- function(buffered_data, input_rast, target_variable){
    extracted_values <- exact_extract(input_rast, buffered_data, target_variable, progress = FALSE)
    out <- buffered_data |>
            mutate(
                ref_val = extracted_values
            )
    return(out)
}