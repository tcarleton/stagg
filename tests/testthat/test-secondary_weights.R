# Test secondary_weights()
test_that("secondary_weights outputs are normal", {

  # Run secondary_weights
  output <- secondary_weights(cropland_nj_2015)

  # Expect the correct column names
  expect_equal(names(output), c("x", "y", "weight"))

  # Expect the correct number of rows
  expect_equal(nrow(output), 182)

  # Expect that area = length * width
  expect_equal(length(unique(output$x)) * length(unique(output$y)),
               182)

  # Expect no NAs with normal output
  expect_true(all(!is.na(output$weight)))

})

test_that("secondary_weights warnings", {

  # Full extent of interest
  nj_extent <- c(-76, -73, 38, 42)

  # Modify the cropland data to test warnings
  ## Make it smaller than extent of interest
  crop_small <- raster::crop(cropland_nj_2015, c(-75, -74, 39, 40))
  ## Add NAs
  crop_na <- raster::reclassify(cropland_nj_2015, cbind(0,0.9,NA))

  # Expect a warning that secondary weights doesn't overlap area of interest
  expect_warning(secondary_weights(secondary_raster = crop_small,
                                   grid = era5_grid,
                                   extent = nj_extent),
                 "Warning: the secondary raster does not fully overlap with the user-specified extent. Resulting data frame will not fully cover user-defined extent.")

  # Expect a warning if secondary weights raster contains NAs
  expect_warning(secondary_weights(secondary_raster = crop_na,
                                   grid = era5_grid,
                                   extent = "full"),
                 "Warning: secondary raster contains NA values. NAs will be returned for weights.")

})

test_that("secondary_weights errors", {

  # Extent type that is not compatible
  nj_counties <- tigris::counties("nj")
  extent_poly <- sf::st_as_sfc(sf::st_bbox(nj_counties))

  # Expect error when extent cannot be used with raster::crop()
  expect_error(secondary_weights(secondary_raster = crop_na,
                                 grid = era5_grid,
                                 extent = extent_poly),
               regexp = "User-defined extent not compatible with raster.")
})
