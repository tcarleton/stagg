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





################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


#                           CODE WORK AHEAD
# !!!       Unverified Tests for the New Code Written Below                  !!!


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


# ==============================================================================
# Internal Helper Functions
# ==============================================================================

# 1. Input Validation and Error Catching
# ______________________________________________________________________________

#   a) validate_data
#   -----------------------------------

# Test that the data types we expect are correctly handled by validate_data()
test_that('validate_data() correctly handles SpatRaster stacks', {

  # Generate testing SpatRaster stack
  test_rast <- terra::rast(vals = 1, ncol = 10, nrow = 10, nlyr = 5)

  # Run validate_data()
  validated_rast <- validate_data(test_rast)

  # Make sure object is still a SpatRaster
  expect_true(inherits(validated_rast, 'SpatRaster'))

  # Make sure layer names preserved
  expect_equal(names(test_rast), names(validated_rast))

  # Make sure values are preserved
  expect_equal(terra::values(test_rast), terra::values(validated_rast))
})


test_that('validate_data() correctly handles raster bricks', {

  # Generate testing raster brick and raster stack
  test_brick <- raster::brick(ncol = 10, nrow = 10, nl = 5)
  raster::values(test_brick) <- 1

  # Run validate_data()
  validated_rast <- validate_data(test_brick)

  # Make sure object was coerced to a SpatRaster
  expect_true(inherits(validated_rast, 'SpatRaster'))

  # Make sure layer names preserved
  expect_equal(names(test_brick), names(validated_rast))

  # Make sure values are preserved
  expect_equal(raster::values(test_brick), terra::values(validated_rast))

})

test_that('validate_data() correctly handles raster stacks', {

  # Generate testing raster brick and raster stack
  test_stack <- raster::brick(ncol = 10, nrow = 10, nl = 5) |>
    raster::stack()
  raster::values(test_stack) <- 1

  # Run validate_data()
  validated_rast <- validate_data(test_stack)

  # Make sure object was coerced to a SpatRaster
  expect_true(inherits(validated_rast, 'SpatRaster'))

  # Make sure layer names preserved
  expect_equal(names(test_stack), names(validated_rast))

  # Make sure values are preserved
  expect_equal(raster::values(test_stack), terra::values(validated_rast))
})


test_that('validate_data() errors when expected', {

  # Vector
  expect_error(
    validate_data(c(1, 2, 3))
  )

  # Missing
  expect_error(
    validate_data()
  )

  # Data Frame
  df <- data.frame(x = c(1,2), y = c(1, 2), value = c(3,4))

  expect_error(
    validate_data(df)
  )

})

#   b) validate_overlay_weights
#   -----------------------------------

#   c) validate_daily_agg
#   -----------------------------------

#   d) validate_time_agg
#   -----------------------------------

#   e) validate_transformations
#   -----------------------------------

#   f) validate_result_cols
#   -----------------------------------

#   g) validate_start_date
#   -----------------------------------

#   h) validate_time_interval
#   -----------------------------------

#   i) validate_weights_join_tolerance
#   -----------------------------------

#   j) validate_na_rm
#   -----------------------------------

# 2. Crop Data to Weights Extent
# ______________________________________________________________________________

#   a) look_for_poly_split
#   -----------------------------------
test_that("look_for_poly_split correctly identifies polygons split by prime merridean", {

  # Create polygons spanning prime meridian
  two_triangles <- gen_two_triangles(edge_case = 'span_pm')

  # Create a grid in climate coordinates with 1x1 resolution
  climate_grid <- terra::rast(
    xmin = 0,
    xmax = 360,
    ncol = 360,
    ymin = -90,
    ymax = 90,
    nrow = 180
  )

  # # Visualize post alignment overlay
  # climate_grid_rotated <- terra::rotate(climate_grid) |>
  #   terra::crop(c(-1, 1, 0, 2))
  # terra::values(climate_grid_rotated) <- c(1,2,3,4)
  # plot_overlay(climate_grid_rotated, two_triangles)

  # Run overlay_weights to generate test input to look_for_poly_split
  overlay_weights_output <- overlay_weights(
    polygons = two_triangles,
    polygon_id_col = 'poly',
    grid = climate_grid
  )

  # Run look_for_poly_split
  polygons_split <- look_for_poly_split(
    data = climate_grid,
    overlay_weights = overlay_weights_output,
    coord_alignment = 'climate'
  )

  expect_true(polygons_split)

})

test_that("look_for_poly_split correctly identifies polygons split by date line", {

  # Create polygons spanning international date line
  two_triangles <- gen_two_triangles(edge_case = 'span_dl')

  # Create a grid in standard coordinates with 1x1 resolution
  standard_grid <- terra::rast(
    xmin = -180,
    xmax = 180,
    ncol = 360,
    ymin = -90,
    ymax = 90,
    nrow = 180
  )

  # # Visualize post alignment overlay
  # standard_grid_rotated <- terra::rotate(standard_grid, left = F) |>
  #   terra::crop(c(179, 181, 0, 2))
  # terra::values(standard_grid_rotated) <- c(1,2,3,4)
  # plot_overlay(standard_grid_rotated, sf::st_shift_longitude(two_triangles))

  # Run overlay_weights to generate test input to look_for_poly_split
  overlay_weights_output <- overlay_weights(
    polygons = two_triangles,
    polygon_id_col = 'poly',
    grid = standard_grid
  )

  # Run look_for_poly_split
  polygons_split <- look_for_poly_split(
    data = standard_grid,
    overlay_weights = overlay_weights_output,
    coord_alignment = 'standard'
  )

  expect_true(polygons_split)

})


test_that("look_for_poly_split correctly returns false when it should", {

  # Create polygons in standard coordinates
  two_triangles <- gen_two_triangles(edge_case = 'none')

  # Create a grid in standard coordinates with 1x1 resolution
  standard_grid <- terra::rast(
    xmin = -180,
    xmax = 180,
    ncol = 360,
    ymin = -90,
    ymax = 90,
    nrow = 180
  )

  # Visualize Overlay
  standard_grid_cropped <- terra::crop(standard_grid, c(10, 12, 0, 2))
  terra::values(standard_grid_cropped) <- c(1,2,3,4)
  plot_overlay(standard_grid_cropped, two_triangles)

  # Run overlay weights to generate test input to look for poly split
  # NOTE: eventually, may want to replace this call to overlay_weights(), in
  # this test and the two above, with the calculated, expected data. That way,
  # if overlay_weights() breaks, it doesn't trigger these tests too, making it
  # easier to pinpoint where the failure occurred. But until overlay_weights()
  # is refactored and unit tested, let's leave as is.
  overlay_weights_output <- overlay_weights(
    polygons = two_triangles,
    polygon_id_col = 'poly',
    grid = standard_grid
  )

  polygons_split <- look_for_poly_split(
    data = standard_grid,
    overlay_weights = overlay_weights_output,
    coord_alignment = 'standard'
  )

  expect_false(polygons_split)

})




