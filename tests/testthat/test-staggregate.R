# Test basic elements of output for helper functions

testthat::test_that(
  "as.data.table.terra has correct output type, class, and column names",
  {
    output <- as_data_table_terra(temp_nj_jun_2024_era5[[1]], xy = TRUE)

    expect_type(output, "list")
    expect_is(output, "data.table")
    expect_setequal(names(output),
                    c("x", "y", names(temp_nj_jun_2024_era5)[1]))
  }
)

testthat::test_that(
  "daily_aggregation has correct output type and class",
  {
    output <- daily_aggregation(temp_nj_jun_2024_era5,
                                overlay_weights_nj,
                                'average',
                                '1 hour')

    # Raster brick element of the list
    expect_type(output[[1]], "S4")
    expect_is(output[[1]], "SpatRaster")

    # Names vector element of the list
    expect_type(output[[2]], "character")
  }
)

testthat::test_that(
  "infer_layer_datetimes has correct output type and class",
  {
    output <- infer_layer_datetimes(temp_nj_jun_2024_era5,
                                    "2024-06-01 00:00:00",
                                    "1 hour")

    expect_type(output, "S4")
    expect_is(output, "SpatRaster")


  }
)


testthat::test_that(
  "polygon_aggregation has correct output type, class, and column names",
  {
    # Make clim_dt input for polygon_aggregation() (the numbers are nonsensical)
    clim_dt_input <- as_data_table_terra(temp_nj_jun_2024_era5, xy = TRUE) %>%
      dplyr::select(c(1,2, 4:10)) %>%
      tidyr::pivot_longer(cols = 3:9, names_to = "date", values_to = "order_1") %>%
      dplyr::mutate(x = x + 360) %>%
      data.table::as.data.table()

    output <- polygon_aggregation(clim_dt_input,
                                  overlay_weights_nj,
                                  c("order_1"),
                                  "hour",
                                  weights_join_tolerance_x = 0,
                                  weights_join_tolerance_y = 0,
                                  na_rm = F)

    expect_type(output, "list")
    expect_is(output, "data.frame")
    expect_setequal(names(output),
                    c("year", "month", "day", "hour", "poly_id", "order_1"))

  }
)


# Test basic elements of user-level functions
testthat::test_that(
  "staggregate_polynomial has correct output type, class, and column names",
  {
    polynomial_output <- staggregate_polynomial(
      data = temp_nj_jun_2024_era5 - 273.15, # Climate data to transform and
      # aggregate
      overlay_weights = overlay_weights_nj, # Output from overlay_weights()
      daily_agg = "average", # Average hourly values to produce daily values
      # before transformation
      time_agg = "hour", # Sum the transformed daily values across months
      start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
      # only required if the layer name
      # format is not compatible with stagg
      time_interval = "1 hour", # The temporal interval of the supplied data,
      # required if daily_agg is not "none" or if the
      # start_date argument is not NA
      degree = 4 # Highest order
    )

    expect_type(polynomial_output, "list")
    expect_is(polynomial_output, "data.frame")
    expect_setequal(names(polynomial_output),
                    c("year", "month", "day", "hour", "poly_id", "order_1",
                      "order_2", "order_3", "order_4"))
  }
)


testthat::test_that(
  "staggregate_spline has correct output type, class, and column names",
  {
    spline_output <- staggregate_spline(
      data = temp_nj_jun_2024_era5 - 273.15, # Climate data to transform and
      # aggregate
      overlay_weights = overlay_weights_nj, # Output from overlay_weights()
      daily_agg = "average", # Average hourly values to produce daily values
      # before transformation
      time_agg = "hour", # Sum the transformed daily values across months
      start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
      # only required if the layer name
      # format is not compatible with stagg
      time_interval = "1 hour", # The temporal interval of the supplied data,
      # required if daily_agg is not "none" or if the
      # start_date argument is not NA
      knot_locs = c(0, 7.5, 12.5, 20) # Highest order
    )

    expect_type(spline_output, "list")
    expect_is(spline_output, "data.frame")
    expect_setequal(names(spline_output),
                    c("year", "month", "day", "hour", "poly_id", "value",
                      "term_1", "term_2"))
  }
)



testthat::test_that(
  "staggregate_bin has correct output type, class, and column names",
  {
    bin_output <- staggregate_bin(
      data = temp_nj_jun_2024_era5 - 273.15, # Climate data to transform and
      # aggregate
      overlay_weights = overlay_weights_nj, # Output from overlay_weights()
      daily_agg = "average", # Average hourly values to produce daily values
      # before transformation
      time_agg = "hour", # Sum the transformed daily values across months
      start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
      # only required if the layer name
      # format is not compatible with stagg
      time_interval = "1 hour", # The temporal interval of the supplied data,
      # required if daily_agg is not "none" or if the
      # start_date argument is not NA
      bin_breaks = c(0, 2.5, 5, 7.5, 10) # Highest order
    )

    expect_type(bin_output, "list")
    expect_is(bin_output, "data.frame")
    expect_setequal(names(bin_output),
                    c('year', 'month', 'day', 'hour', 'poly_id',
                      'bin_ninf_to_0', 'bin_0_to_2.5', 'bin_2.5_to_5',
                      'bin_5_to_7.5', 'bin_7.5_to_10', 'bin_10_to_inf'))
  }
)


testthat::test_that(
  "staggregate_degree_days has correct output type, class, and column names",
  {
    degree_days_output <- staggregate_degree_days(
      data = temp_nj_jun_2024_era5 - 273.15, # Climate data to transform and
      # aggregate
      overlay_weights = overlay_weights_nj,
      time_agg = "day", # Sum the transformed daily values across months
      start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
      # only required if the layer name
      # format is not compatible with stagg
      time_interval = "1 hour", # The temporal interval of the supplied data,
      # required if daily_agg is not "none" or if the
      # start_date argument is not NA
      thresholds = c(0, 10, 20) # Highest order
    )

    expect_type(degree_days_output, "list")
    expect_is(degree_days_output, "data.frame")
    expect_setequal(names(degree_days_output),
                    c('year', 'month', 'day', 'poly_id', 'threshold_ninf_to_0', 'threshold_0_to_10', 'threshold_10_to_20', 'threshold_20_to_inf'))
  }
)
