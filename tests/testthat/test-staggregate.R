# Test that as.data.table.raster conversion is the same
testthat::test_that("as.data.table.raster works", {


  as_data_table_raster_output <- as.data.table.raster(temp_nj_jun_2024_era5[[1]],
                                                      xy = TRUE) %>%
    dplyr::arrange(x, y)

  as_data_table_raster_key <- readRDS(
    testthat::test_path("fixtures/as_data_table_raster_key.rds")
  ) %>%
    dplyr::arrange(x,y)

  expect_true(all.equal(as_data_table_raster_output, as_data_table_raster_key))


})


# Test staggregate_polynomial()
testthat::test_that("staggregate_polynomial works", {

  staggregate_polynomial_output <- staggregate_polynomial(
    data = temp_nj_jun_2024_era5 - 273.15, # Climate data to transform and aggregate
    overlay_weights = overlay_weights_nj, # Output from overlay_weights()
    daily_agg = "average", # Average hourly values to produce daily values
    # before transformation
    time_agg = "month", # Sum the transformed daily values across months
    start_date = "2020-01-01 00:00:00", # The start date of the supplied data, only required if the layer name format is not compatible with stagg
    time_interval = "1 hour", # The temporal interval of the supplied data, required if daily_agg is not "none" or if the start_date argument is not NA
    degree = 4 # Highest order
  )

  staggregate_polynomial_key <- readRDS(
    testthat::test_path("fixtures/staggregate_polynomial_key.rds")
  )

  expect_equal(staggregate_polynomial_output, staggregate_polynomial_key)
})
