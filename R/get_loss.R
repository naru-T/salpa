#' Get loss value to be minimized
#'
#' @param lidar_footprints sf object with satellite LiDAR footprints
#' @param add_x numeric value to add to x coordinates
#' @param add_y numeric value to add to y coordinates
#' @param buf numeric value for buffer size
#' @param input_rast raster object for reference data
#' @param minimizing_method character string for the minimizing method. Default is "euclidean".
#'        Options are "dtw", "euclidean", "manhattan", "correlation", "area", "hausdorff".
#' @param target_variable character string for the target variable
#' @param crs_code numeric value for CRS code
#' @param lidar_value name of the column containing LiDAR measurement values
#' @importFrom sf st_drop_geometry
#' @importFrom sf st_buffer
#' @return data frame with loss value for the given offset
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
#'   result <- get_loss(lidar_sf, add_x = 0, add_y = 0, buf = 10,
#'                      input_rast = dem, target_variable = "mean",
#'                      crs_code = 32654, lidar_value = "elev")
#' }
#' @export
get_loss <- function(lidar_footprints, add_x, add_y, buf, input_rast,
                    minimizing_method = "euclidean", target_variable,
                    crs_code, lidar_value) {
    posit_corrected_data <- position_adjustment(lidar_footprints, add_x, add_y, crs_code)
    buffered_data <- st_buffer(posit_corrected_data, buf)
    out <- extract_values(buffered_data, input_rast, target_variable)

    # Extract LiDAR and reference values, ensuring they correspond to same footprints
    lidar_val <- lidar_footprints[[lidar_value]] |> st_drop_geometry() |> as.numeric()
    ref_val <- out[["ref_val"]] |> st_drop_geometry() |> as.numeric()
    
    # Filter out NA values to ensure both sequences have same length and valid data
    valid_indices <- !is.na(lidar_val) & !is.na(ref_val) & is.finite(lidar_val) & is.finite(ref_val)
    
    if (sum(valid_indices) < 1) {
        # No valid data pairs available
        loss_value <- 1000  # High penalty value
        message("[get_loss] No valid data pairs: add_x=", add_x, " add_y=", add_y)
    } else {
        lidar_val_clean <- lidar_val[valid_indices]
        ref_val_clean <- ref_val[valid_indices]
        
        # Calculate distance using cleaned, matched sequences
        distance_result <- perform_distance(lidar_val_clean, ref_val_clean, method = minimizing_method)
        loss_value <- distance_result$distance
        
        # Handle case where distance calculation fails
        if (is.na(loss_value) || !distance_result$valid) {
            loss_value <- 1000  # High penalty value
            message("[get_loss] Distance calculation failed: add_x=", add_x, " add_y=", add_y)
        }
    }
    
    loss_ <- data.frame(loss_value = loss_value)

    out <- data.frame(
        add_x = add_x,
        add_y = add_y,
        buffer = buf,
        loss_value = loss_$loss_value
    )

    return(out)
}
