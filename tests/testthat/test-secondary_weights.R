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



################################################################################
# Tests to keep after refactoring
################################################################################

# Fix issue 50
test_that('secondary_weights works with mismatched projections', {

  # Create test data
  test_secondary_raster <- terra::rast(cropland_nj_2015) |>
    terra::project("ESRI:54008")

  # Run secondary_weights on test data
  test_output <- secondary_weights(test_secondary_raster)

  # Get validation data
  validation_secondary_raster <- terra::rast(cropland_nj_2015)
  validation_output <- secondary_weights(validation_secondary_raster, grid = era5_grid)

  # Join
  joined_output <- dplyr::full_join(
    dplyr::rename(test_output, test_weight = weight),
    dplyr::rename(validation_output, validation_weight = weight),
    by = c('x', 'y')
  ) |>
    dplyr::mutate(diff = test_weight - validation_weight) |>
    dplyr::filter(!is.na(test_weight), !is.na(validation_weight))


  # ggplot2::ggplot(joined_output) +
  #   ggplot2::geom_tile(ggplot2::aes(x,y, fill = test_weight))
  #
  # ggplot2::ggplot(joined_output) +
  #   ggplot2::geom_tile(ggplot2::aes(x,y, fill = validation_weight))
  #
  # ggplot2::ggplot(joined_output) +
  #   ggplot2::geom_histogram(ggplot2::aes(x = diff))
  #
  # ggplot2::ggplot(joined_output) +
  #   ggplot2::geom_point(ggplot2::aes(x = validation_weight, y = test_weight)) +
  #   ggplot2::scale_color_viridis_c()



  # Compute regression and make sure coefficients are within reason
  reg <- lm(test_weight ~ 0 + validation_weight + x + y, joined_output)

  # Location should play little role
  expect_lt(reg$coefficients[['x']], 0.1)
  expect_lt(reg$coefficients[['y']], 0.1)

  # validation_weight should be near 1
  expect_gt(reg$coefficients[['validation_weight']], 0.9)
  expect_lt(reg$coefficients[['validation_weight']], 1.1)

  # Should be fairly strong correlation
  r_squared <- cor(joined_output$validation_weight, joined_output$test_weight)^2
  expect_gt(r_squared, 0.8)
})
