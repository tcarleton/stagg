# Test overlay weights
test_that("Data structure of overlay_weights is normal", {

  load("data/overlay_weights_nj.rda")

  # Expect the correct column names
  expect_equal(names(overlay_weights_nj), c("x", "y", "poly_id", "w_area", "weight"))

  # Expect the correct number of rows
  expect_equal(nrow(overlay_weights_nj), 141)

  # Expect the correct number of unique counties
  expect_equal(length(unique(overlay_weights_nj$poly_id)), 21)

})

test_that("overlay_weights outputs are normal", {

  # Load data to run in overlay_weights
  load("data/era5_grid.rda")
  load("data/cropland_nj_2015.rda")
  nj <- tigris::counties(state="New Jersey")

  # Run secondary_weights
  crop_weights <- secondary_weights(cropland_nj_2015)
  crop_weights_na <- crop_weights |>
    # Set all odd rows to NA
    dplyr::mutate(weight = ifelse(x==-74.5, NA, weight))

  # Run overlay_weights normal
  normal_output <- overlay_weights(polygons = nj,
                                   polygon_id_col = "COUNTYFP",
                                   grid = era5_grid,
                                   secondary_weights = crop_weights)

  sum_normal <- normal_output |>
    group_by(poly_id) |>
    summarize(w_area = sum(w_area),
              weight = sum(weight)) |>
    ungroup()

  # Expect that area weights are between 0 and 1
  ## Smallest value in the area weights column > 0
  expect_gt(min(normal_output$w_area), 0)
  ## Largest value in the area weights column <= 1
  expect_lte(max(normal_output$w_area), 1)

  # Expect that area + secondary weights are between 0 and 1
  ## Smallest value in the weights column > 0
  expect_gt(min(normal_output$weight), 0)
  ## Largest value in the weights column <= 1
  expect_lte(max(normal_output$weight), 1)

  # Normal output
  ## Expect that area weights sum to 1 for each county
  expect_true(all(round(sum_normal$w_area) == 1))

  ## Expect that area + secondary weights sum to 1 by county
  expect_true(all(round(sum_normal$weight) == 1))

})

test_that("overlay_weights works when some secondary weights are NA", {

  # Load data to run in overlay_weights
  load("data/era5_grid.rda")
  load("data/cropland_nj_2015.rda")
  nj <- tigris::counties(state="New Jersey")

  # Run secondary_weights
  crop_weights <- secondary_weights(cropland_nj_2015)
  crop_weights_na <- crop_weights |>
    # Set all odd rows to NA
    dplyr::mutate(weight = ifelse(x==-74.5, NA, weight))

  # Expect a warning that secondary weights contains one more NAs
  expect_warning(overlay_weights(polygons = nj,
                                 polygon_id_col = "COUNTYFP",
                                 grid = era5_grid,
                                 secondary_weights = crop_weights_na))

  # Run overlay_weights with NAs
  options(warn=-1)  # This won't actually run within the testing code because of the warning so we supress it
  na_output <- overlay_weights(polygons = nj,
                               polygon_id_col = "COUNTYFP",
                               grid = era5_grid,
                               secondary_weights = crop_weights_na)

  sum_na <- na_output |>
    group_by(poly_id) |>
    summarize(w_area = sum(w_area),
              weight = sum(weight)) |>
    ungroup()

  # Expect that area weights sum to 1 for each county (rounded because expect true requires identical matching)
  expect_true(all(round(sum_na$w_area) == 1))

  # Expect that area + secondary weights are either 1 or NA
  expect_true(all(round(sum_na$weight) == 1 | is.na(sum_na$weight)))

})

test_that("overlay_weights errors", {

  # Load data to run in overlay_weights
  load("data/era5_grid.rda")
  load("data/cropland_nj_2015.rda")
  nj <- tigris::counties(state="New Jersey")

  # Shift the nj polygons so they are 0-360
  nj_shift <- nj |> sf::st_shift_longitude()

  # Errors: 3 stops in the overlay_weight function
  ## Error if polygons are in 0-360
  expect_error(overlay_weights(polygons = nj_shift,
                               polygon_id_col = "COUNTYFP",
                               grid = era5_grid))

  ## Error if area weights or weights don't sum to 1: can't break this since it
  ## requires changing the function itself not the input data
  ## But checking outputs sum to 1 is included above


})

