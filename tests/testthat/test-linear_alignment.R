test_that("linear_alignment works with valid input", {
  # Load sample data
  gpkg_path <- system.file("extdata", "gedi_l2a_shot_sample.gpkg", package = "salpa")
  lidar_footprints <- sf::st_read(gpkg_path, quiet = TRUE)
  
  # Run alignment
  result <- linear_alignment(lidar_footprints, 3857)
  
  # Test output structure
  expect_s3_class(result, "sf")
  expect_true("orig_x" %in% names(result))
  expect_true("orig_y" %in% names(result))
  expect_equal(nrow(result), nrow(lidar_footprints))
})

test_that("linear_alignment handles invalid input appropriately", {
  # Test missing shot_number column
  gpkg_path <- system.file("extdata", "gedi_l2a_shot_sample.gpkg", package = "salpa")
  lidar_footprints <- sf::st_read(gpkg_path, quiet = TRUE)
  bad_data <- lidar_footprints[, !names(lidar_footprints) %in% "shot_number"]
  
  expect_error(linear_alignment(bad_data, 3857),
               "The sf_object does not contain a 'shot_number' column")
  
  # Test non-numeric CRS
  expect_error(linear_alignment(lidar_footprints, "invalid"),
               "crs_code must be a numeric EPSG code")
               
  # Test invalid numeric CRS
  expect_error(linear_alignment(lidar_footprints, 999999),
               "Invalid crs_code: 999999")
}) 