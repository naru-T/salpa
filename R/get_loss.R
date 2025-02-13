#' Get loss value to be minimized
#'
#' @param lidar_footprints sf object with satellite LiDAR footprints
#' @param add_x numeric value to add to x coordinates
#' @param add_y numeric value to add to y coordinates
#' @param buf numeric value for buffer size
#' @param input_raster raster object for extraction
#' @param minimizing_method character string for the minimizing method
#' @param target_variable character string for the target variable
#' @param crs_code numeric value for CRS code
#' @param lidar_value name of the column containing LiDAR measurement values
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr summarise
#' @return data frame with summary statistics
#' @export
get_loss <- function(lidar_footprints, add_x, add_y, buf, input_raster, minimizing_method, target_variable, crs_code, lidar_value){
    buffered_data <- positional_correction(lidar_footprints, add_x, add_y, buf, crs_code)
    out <- extract_values(buffered_data, input_raster, target_variable)
    
    lidar_val <- lidar_footprints[[lidar_value]] |> st_drop_geometry() |> as.numeric()
    ref_val <- out[["ref_val"]] |> st_drop_geometry() |> as.numeric()
    #WIP: option to apply different methods according to minimizing_method
    loss_ <- out |>
                    st_drop_geometry() |>
                    summarise(
                        loss_value = perform_dtw(lidar_val, ref_val)
                    )

    out <- data.frame(
        add_x = add_x,
        add_y = add_y,
        buffer = buf,
        loss_value = loss_$loss_value
    )
    
    return(out)
}