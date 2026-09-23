test_that("validate_crs_code works correctly", {
  # Valid CRS codes should not throw errors
  expect_no_error(validate_crs_code(4326))
  expect_no_error(validate_crs_code(3857))
  
  # Invalid CRS codes should throw errors
  expect_error(validate_crs_code("invalid"))
  expect_error(validate_crs_code(999999))
})

test_that("validate_buffer works correctly", {
  # Valid buffer sizes should not throw errors
  expect_true(validate_buffer(1))
  expect_true(validate_buffer(12.5))
  expect_true(validate_buffer(100))
  
  # Invalid buffer sizes should throw errors
  expect_error(validate_buffer(0))
  expect_error(validate_buffer(-1))
})

test_that("validate_optimization_method works correctly", {
  # Valid methods should not throw errors
  expect_true(validate_optimization_method("ga"))
  expect_true(validate_optimization_method("pso"))
  expect_true(validate_optimization_method("woa"))
  expect_true(validate_optimization_method("lbfgsb"))
  
  # Invalid methods should throw errors
  expect_error(validate_optimization_method("invalid"))
  expect_error(validate_optimization_method("genetic"))
})

test_that("convert_bounds_to_degrees works correctly", {
  # Test conversion
  result <- convert_bounds_to_degrees(c(1000, 2000))
  expected <- c(1000, 2000) / 111320
  expect_equal(result, expected)
  
  # Test with negative values
  result_neg <- convert_bounds_to_degrees(c(-1000, -2000))
  expected_neg <- c(-1000, -2000) / 111320
  expect_equal(result_neg, expected_neg)
})

test_that("validate_optimization_params works correctly", {
  # Valid parameters should not throw errors
  result <- validate_optimization_params(50, 100, c(-30, -30), c(30, 30))
  expect_type(result, "list")
  expect_equal(result$pop_size, 50L)
  expect_equal(result$max_iter, 100L)
  
  # Invalid parameters should throw errors
  expect_error(validate_optimization_params(0, 100, c(-30, -30), c(30, 30)))  # pop_size <= 0
  expect_error(validate_optimization_params(50, 0, c(-30, -30), c(30, 30)))   # max_iter <= 0
  expect_error(validate_optimization_params(50, 100, c(-30), c(30, 30)))      # length mismatch
  expect_error(validate_optimization_params(50, 100, c(30, 30), c(-30, -30))) # lower >= upper
})