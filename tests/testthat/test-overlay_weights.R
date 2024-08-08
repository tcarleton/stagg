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


