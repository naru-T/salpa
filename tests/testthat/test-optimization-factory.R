test_that("optimization_factory validates inputs correctly", {
  # Create a simple objective function for testing
  test_objective <- function(x, y) { (x - 1)^2 + (y - 2)^2 }
  
  # Valid inputs should not throw validation errors
  expect_no_error({
    optimization_factory(
      optimization_method = "lbfgsb",
      objective_function = test_objective,
      lower_bounds = c(-10, -10),
      upper_bounds = c(10, 10),
      pop_size = 20,
      max_iter = 10
    )
  })
  
  # Invalid optimization method should throw error
  expect_error(
    optimization_factory(
      optimization_method = "invalid_method",
      objective_function = test_objective,
      lower_bounds = c(-10, -10),
      upper_bounds = c(10, 10),
      pop_size = 20,
      max_iter = 10
    )
  )
  
  # Invalid parameters should throw error
  expect_error(
    optimization_factory(
      optimization_method = "ga",
      objective_function = test_objective,
      lower_bounds = c(-10, -10),
      upper_bounds = c(10, 10),
      pop_size = 0,  # Invalid
      max_iter = 10
    )
  )
})

test_that("optimization_factory returns expected structure", {
  # Create a simple objective function
  test_objective <- function(x, y) { (x - 1)^2 + (y - 2)^2 }
  
  # Test L-BFGS-B (fastest, should converge quickly)
  result <- optimization_factory(
    optimization_method = "lbfgsb",
    objective_function = test_objective,
    lower_bounds = c(-10, -10),
    upper_bounds = c(10, 10),
    pop_size = 20,
    max_iter = 50
  )
  
  # Check result structure
  expect_s3_class(result, "data.frame")
  expect_true("best_x" %in% names(result))
  expect_true("best_y" %in% names(result))
  expect_true("best_value" %in% names(result))
  
  # Check that it found a reasonable solution (should be close to x=1, y=2)
  # More tolerant bounds for different optimization algorithms
  expect_true(abs(result$best_x - 1) < 2.0)
  expect_true(abs(result$best_y - 2) < 2.0)
  expect_true(result$best_value < 10)  # Should be reasonably low
})