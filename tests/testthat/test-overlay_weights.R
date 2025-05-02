# Test overlay weights
test_that("overlay_weights outputs are normal", {

  nj_counties <- tigris::counties("nj")

  # Run secondary_weights
  crop_weights <- secondary_weights(cropland_nj_2015)

  # Run overlay_weights normal
  normal_output <- overlay_weights(polygons = nj_counties,
                                   polygon_id_col = "COUNTYFP",
                                   grid = era5_grid,
                                   secondary_weights = crop_weights)
  # And sum by polygon to check outputs
  sum_normal <- normal_output %>%
    dplyr::group_by(poly_id) %>%
    dplyr::summarize(w_area = sum(w_area),
              weight = sum(weight)) %>%
    dplyr::ungroup()

  # Expect the correct column names
  expect_equal(names(normal_output), c("x", "y", "poly_id", "w_area", "weight"))

  # Expect the correct number of rows
  expect_equal(nrow(normal_output), 141)

  # Expect the correct number of unique counties
  expect_equal(length(unique(normal_output$poly_id)), 21)

  # Expect that area weights are between 0 and 1
  ## Smallest value in the area weights column > 0
  expect_gt(min(normal_output$w_area), 0)
  ## Largest value in the area weights column <= 1
  expect_lte(max(normal_output$w_area), 1)

  # Expect that area + secondary weights are between 0 and 1
  ## Smallest value in the weights column > 0
  expect_gte(min(normal_output$weight), 0)
  ## Largest value in the weights column <= 1
  expect_lte(max(normal_output$weight), 1)

  # Expect that area weights sum to 1 for each county
  expect_true(all(round(sum_normal$w_area) == 1))

  # Expect that area + secondary weights sum to 1 by county
  expect_true(all(round(sum_normal$weight) == 1))

})

test_that("overlay_weights works when some secondary weights are NA", {

  nj_counties <- tigris::counties("nj")

  # Run secondary_weights
  crop_weights_na <- secondary_weights(cropland_nj_2015)  %>%
    # Set all odd rows to NA
    dplyr::mutate(weight = ifelse(x==-74.5, NA, weight))

  # Run overlay_weights with NAs
  options(warn=-1)  # This won't actually run within the testing code because of the warning so we suppress it
  na_output <- overlay_weights(polygons = nj_counties,
                               polygon_id_col = "COUNTYFP",
                               grid = era5_grid,
                               secondary_weights = crop_weights_na)

  sum_na <- na_output %>%
    dplyr::group_by(poly_id) %>%
    dplyr::summarize(w_area = sum(w_area),
              weight = sum(weight)) %>%
    dplyr::ungroup()

  # Expect that area weights sum to 1 for each county (rounded because expect true requires identical matching)
  expect_true(all(round(sum_na$w_area) == 1))

  # Expect that area + secondary weights are either 1 or NA
  expect_true(all(round(sum_na$weight) == 1 | is.na(sum_na$weight)))

})

test_that("overlay_weights warnings", {

  nj_counties <- tigris::counties("nj")

  # Run secondary_weights
  crop_weights_na <- secondary_weights(cropland_nj_2015) %>%
    # Set all odd rows to NA
    dplyr::mutate(weight = ifelse(x==-74.5, NA, weight))

  # Make the extent of secondary weights smaller than polygons
  crop_weights_small <- secondary_weights(cropland_nj_2015) %>%
    dplyr::filter(y >= 39 & y< 41)

  # Expect a warning that secondary weights contains one more NAs
  expect_warning(overlay_weights(polygons = nj_counties,
                                 polygon_id_col = "COUNTYFP",
                                 grid = era5_grid,
                                 secondary_weights = crop_weights_na),
                 "Warning: secondary weight values contain one or more NAs. The resulting weights for x,y coordinates with NA secondary weight values will be NAs.")

  # Expect a warning that secondary weights doesn't overlap fully with polygons
  expect_warning(overlay_weights(polygons = nj_counties,
                                 polygon_id_col = "COUNTYFP",
                                 grid = era5_grid,
                                 secondary_weights = crop_weights_small),
                 "Warning: secondary weights do not fully overlap with the administrative regions. Resulting weights will contain NAs.")

})

test_that("overlay_weights errors", {

  nj_counties <- tigris::counties("nj")

  # Shift the nj_countiespolygons so they are 0-360
  nj_shift <- nj_counties %>% sf::st_shift_longitude()

  # Errors: 3 stops in the overlay_weight function
  ## Error if polygons are in 0-360
  expect_error(overlay_weights(polygons = nj_shift,
                               polygon_id_col = "COUNTYFP",
                               grid = era5_grid))

  ## Error if area weights or weights don't sum to 1: can't break this since it
  ## requires changing the function itself not the input data
  ## But checking outputs sum to 1 is included above


})

