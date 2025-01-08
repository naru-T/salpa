#' Extract values of input raster data in the buffered data
#'
#' @param buffered_data sf object with buffered data
#' @param input_raster raster object for extraction
#' @param target_variable character string for the target variable (min, max, mean, median, etc.). See exact_extract
#' @return sf object with summary statistics
#' @importFrom exactextractr exact_extract
#' @importFrom dplyr mutate
#' @export
extract_values <- function(buffered_data, input_raster, target_variable){
    extracted_values <- exact_extract(input_raster, buffered_data, target_variable, progress = FALSE)
    out <- buffered_data |>
            mutate(
                ref_val = extracted_values
            )
    return(out)
}