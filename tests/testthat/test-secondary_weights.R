# Test secondary_weights() without checking actual values in output
test_that("Data structure of secondary_weights is normal", {

  # Run secondary_weights
  output <- secondary_weights(cropland_nj_2015)

  # Expect the correct column names
  expect_equal(names(output), c("x", "y", "weight"))

  # Expect the correct number of rows
  expect_equal(nrow(output), 120)

  # Expect that area = length * width
  expect_equal(length(unique(output$x)) * length(unique(output$y)),
               120)

  # Expect that weights are between 0 and 1
  expect_equal(sum(output$weight > 1 | output$weight < 0))



  # !!! THIS LAST EXPECTATION FAILS !!! The edges of the raster can sometimes,
  # due to bilinear interpolation, be negative :(
})

