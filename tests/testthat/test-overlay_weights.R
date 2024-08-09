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

test_that("Data values of overlay_weights are normal", {

  load("data/overlay_weights_nj.rda")
  sum_by_poly <- overlay_weights_nj |>
    group_by(poly_id) |>
    summarize(w_area = sum(w_area),
              weight = sum(weight)) |>
    ungroup()

  # Expect that area weights are between 0 and 1
  ## Smallest value in the area weights column > 0
  expect_gt(min(overlay_weights_nj$w_area), 0)
  ## Largest value in the area weights column <= 1
  expect_lte(max(overlay_weights_nj$w_area), 1)

  # Expect that area + secondary weights are between 0 and 1
  ## Smallest value in the weights column > 0
  expect_gt(min(overlay_weights_nj$weight), 0)
  ## Largest value in the weights column <= 1
  expect_lte(max(overlay_weights_nj$weight), 1)

  # Expect that area weights sum to 1 for each county
  expect_equal(sum_by_poly$w_area, rep(1, nrow(sum_by_poly)))

  # Expect that area + secondary weights sum to 1 by county
  expect_equal(sum_by_poly$weight, rep(1, nrow(sum_by_poly)))

})

test_that("overlay_weights errors", {

  # Load data to run in overlay_weights
  load("data/era5_grid.rda")
  load("data/cropland_nj_2015.rda")
  nj <- tigris::counties(state="New Jersey")
  # Run secondary_weights
  crop_weights <- secondary_weights(cropland_nj_2015) |>
    # Set all odd rows to NA
    dplyr::mutate(weight = ifelse(x==-74.5, NA, weight))

  # Shift the nj polygons so they are 0-360
  nj_shift <- nj |> sf::st_shift_longitude()

  # Errors: 4 stops in the overlay_weight function
  ## Error if polygons are in 0-360
  expect_error(overlay_weights(polygons = nj_shift,
                               polygon_id_col = "COUNTYFP",
                               grid = era5_grid))

  ## Error if weights don't sum to 1 or all cells are NA
  expect_error(overlay_weights(polygons = nj,
                               polygon_id_col = "COUNTYFP",
                               grid = era5_grid,
                               secondary_weights = crop_weights)) # Still doesn't error?
  # Always returns NA for the entire polygon
  # Investing in overlay weights it looks like we change every cell to NA within a certain polygon id
  # if one cell is NA in the weights.

  ## Error if area weights or weights don't sum to 1: can't break this since it
  ## requires changing the function itself not the input data
  ## But checking outputs sum to 1 is included above


})
