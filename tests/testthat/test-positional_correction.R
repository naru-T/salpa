test_that("positional_correction works with valid input", {
  # Load sample data
  gpkg_path <- system.file("extdata", "gedi_l2a_shot_sample.gpkg", package = "salpa")
  lidar_footprints <- sf::st_read(gpkg_path, quiet = TRUE)
  
  # Run correction
  result <- positional_correction(lidar_footprints, 10, 10, 12.5, 3857)
  
  # Test output structure
  expect_s3_class(result, "sf")
  expect_true(sf::st_geometry_type(result)[1] == "POLYGON")
  expect_equal(nrow(result), nrow(lidar_footprints))
})

test_that("positional_correction handles edge cases", {
  # Load sample data
  gpkg_path <- system.file("extdata", "gedi_l2a_shot_sample.gpkg", package = "salpa")
  lidar_footprints <- sf::st_read(gpkg_path, quiet = TRUE)
  
  # Test zero offset
  result_zero <- positional_correction(lidar_footprints, 0, 0, 12.5, 3857)
  expect_s3_class(result_zero, "sf")
  
  # Test negative buffer
  expect_error(positional_correction(lidar_footprints, 10, 10, -1, 3857),
               "Buffer size must be a positive number")
  
  # Test zero buffer
  expect_error(positional_correction(lidar_footprints, 10, 10, 0, 3857),
               "Buffer size must be a positive number")
  
  # Test non-numeric CRS
  expect_error(positional_correction(lidar_footprints, 10, 10, 12.5, "invalid"),
               "crs_code must be a numeric EPSG code")
               
  # Test invalid numeric CRS
  expect_error(positional_correction(lidar_footprints, 10, 10, 12.5, 999999),
               "Invalid crs_code: 999999")
}) 