test_that("positional_correction works with valid input", {
  # Load sample data
  gpkg_path <- system.file("extdata", "gedi_l2a_shot_sample.gpkg", package = "salpa")
  lidar_footprints <- sf::st_read(gpkg_path, quiet = TRUE)
  
  # Create a dummy raster for testing
  input_rast <- terra::rast(matrix(1:100, nrow = 10, ncol = 10))
  terra::ext(input_rast) <- terra::ext(lidar_footprints)
  terra::crs(input_rast) <- "EPSG:3857"
  
  # Run correction
  result <- positional_correction(lidar_footprints, input_rast, target_variable = "mean", crs_code = 3857, lidar_value = "rh98")
  
  # Test output structure
  expect_s3_class(result$corrected_footprints, "sf")
  expect_true(sf::st_geometry_type(result$corrected_footprints)[1] == "POINT")
  expect_equal(nrow(result$corrected_footprints), nrow(lidar_footprints))
})

test_that("positional_correction handles edge cases", {
  # Load sample data
  gpkg_path <- system.file("extdata", "gedi_l2a_shot_sample.gpkg", package = "salpa")
  lidar_footprints <- sf::st_read(gpkg_path, quiet = TRUE)
  
  # Create a dummy raster for testing
  input_rast <- terra::rast(matrix(1:100, nrow = 10, ncol = 10))
  terra::ext(input_rast) <- terra::ext(lidar_footprints)
  terra::crs(input_rast) <- "EPSG:3857"
  
  # Test with default parameters
  result_default <- positional_correction(lidar_footprints, input_rast, target_variable = "mean", crs_code = 3857, lidar_value = "rh98")
  expect_s3_class(result_default$corrected_footprints, "sf")
  
  # Test negative buffer
  expect_error(positional_correction(lidar_footprints, input_rast, target_variable = "mean", buf = -1, crs_code = 3857, lidar_value = "rh98"),
               "Buffer size must be a positive number")
  
  # Test zero buffer
  expect_error(positional_correction(lidar_footprints, input_rast, target_variable = "mean", buf = 0, crs_code = 3857, lidar_value = "rh98"),
               "Buffer size must be a positive number")
  
  # Test non-numeric CRS
  expect_error(positional_correction(lidar_footprints, input_rast, target_variable = "mean", crs_code = "invalid", lidar_value = "rh98"),
               "crs_code must be a numeric EPSG code")
               
  # Test invalid numeric CRS
  expect_error(positional_correction(lidar_footprints, input_rast, target_variable = "mean", crs_code = 999999, lidar_value = "rh98"),
               "Invalid crs_code: 999999 - Resulting CRS has no valid EPSG.")
}) 