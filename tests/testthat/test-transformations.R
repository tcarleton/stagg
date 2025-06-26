# ==============================================================================
# Comparing Old and New Staggregate Outputs (will delete after verifying)
# ==============================================================================

# Staggregate Polynomial
test_that('old staggregate_polynomial output is the same as new', {

  # Run old example
  polynomial_output <- staggregate_polynomial(
    data = terra::rast(temp_nj_jun_2024_era5) - 273.15, # Climate data to transform and
    # aggregate
    overlay_weights = overlay_weights_nj, # Output from overlay_weights()
    daily_agg = "average", # Average hourly values to produce daily values
    # before transformation
    time_agg = "month", # Sum the transformed daily values across months
    start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
    # only required if the layer name
    # format is not compatible with stagg
    time_interval = "1 hour", # The temporal interval of the supplied data,
    # required if daily_agg is not "none" or if the
    # start_date argument is not NA
    degree = 4 # Highest order
  )

  # Run new example
  polynomial_output_new <- staggregate_polynomial_new(
    data = terra::rast(temp_nj_jun_2024_era5) - 273.15, # Climate data to transform and
    # aggregate
    overlay_weights = overlay_weights_nj, # Output from overlay_weights()
    daily_agg = "average", # Average hourly values to produce daily values
    # before transformation
    time_agg = "month", # Sum the transformed daily values across months
    start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
    # only required if the layer name
    # format is not compatible with stagg
    time_interval = "1 hour", # The temporal interval of the supplied data,
    # required if daily_agg is not "none" or if the
    # start_date argument is not NA
    degree = 4 # Highest order
  )

  expect_true(
    all.equal(polynomial_output, polynomial_output_new)
  )
})

# Staggregate Spline
test_that('old staggregate_spline output is the same as new', {

  # Run old example
  spline_output <- staggregate_spline(
    data = terra::rast(temp_nj_jun_2024_era5) - 273.15, # Climate data to transform and
    # aggregate
    overlay_weights = overlay_weights_nj, # Output from overlay_weights()
    daily_agg = "average", # Average hourly values to produce daily values
    # before transformation
    time_agg = "month", # Sum the transformed daily values across months
    start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
    # only required if the layer name format
    # is not compatible with stagg
    time_interval = "1 hour", # The temporal interval of the supplied data,
    # required if daily_agg is not "none" or if the
    # start_date argument is not NA
    knot_locs = c(0, 7.5, 12.5, 20) # Where to place knots
  ) |>
    # Modify some of the names
    dplyr::rename(spline_term_1 = term_1, spline_term_2 = term_2, orig_value = value)

  # Run new example
  spline_output_new <- staggregate_spline_new(
    data = terra::rast(temp_nj_jun_2024_era5) - 273.15, # Climate data to transform and
    # aggregate
    overlay_weights = overlay_weights_nj, # Output from overlay_weights()
    daily_agg = "average", # Average hourly values to produce daily values
    # before transformation
    time_agg = "month", # Sum the transformed daily values across months
    start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
    # only required if the layer name format
    # is not compatible with stagg
    time_interval = "1 hour", # The temporal interval of the supplied data,
    # required if daily_agg is not "none" or if the
    # start_date argument is not NA
    knot_locs = c(0, 7.5, 12.5, 20) # Where to place knots
  )

  expect_true(
    all.equal(spline_output, spline_output_new)
  )
})

# Staggregate bin
test_that('old staggregate_bin output is the same as new',{

  # Run old example
  bin_output <- staggregate_bin(
    data = terra::rast(temp_nj_jun_2024_era5) - 273.15, # Climate data to transform and
    # aggregate
    overlay_weights = overlay_weights_nj, # Output from overlay_weights()
    daily_agg = "average", # Average hourly values to produce daily values
    # before transformation
    time_agg = "month", # Sum the transformed daily values across months
    start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
    # only required if the layer name
    # format is not compatible with stagg
    time_interval = "1 hour", # The temporal interval of the supplied data,
    # required if daily_agg is not "none" or if the
    # start_date argument is not NA
    bin_breaks = c(0, 2.5, 5, 7.5, 10) # Draw 6 bins from ninf to 0, 0 to 2.5,
    # 2.5 to 5, 5 to 7.5, 7.5 to 10, 10 to
    # inf
  )

  # Run new example
  bin_output_new <- staggregate_bin_new(
    data = terra::rast(temp_nj_jun_2024_era5) - 273.15, # Climate data to transform and
    # aggregate
    overlay_weights = overlay_weights_nj, # Output from overlay_weights()
    daily_agg = "average", # Average hourly values to produce daily values
    # before transformation
    time_agg = "month", # Sum the transformed daily values across months
    start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
    # only required if the layer name
    # format is not compatible with stagg
    time_interval = "1 hour", # The temporal interval of the supplied data,
    # required if daily_agg is not "none" or if the
    # start_date argument is not NA
    bin_breaks = c(0, 2.5, 5, 7.5, 10) # Draw 6 bins from ninf to 0, 0 to 2.5,
    # 2.5 to 5, 5 to 7.5, 7.5 to 10, 10 to
    # inf
  )

  expect_true(
    all.equal(bin_output, bin_output_new)
  )
})

# Staggregate degree days
test_that('old staggregate_degree_days output is the same as new', {

  # Run old example
  degree_days_output <- staggregate_degree_days(
    data = terra::rast(temp_nj_jun_2024_era5) - 273.15, # Climate data to transform and
    # aggregate
    overlay_weights = overlay_weights_nj, # Output from overlay_weights()
    time_agg = "month", # Sum the transformed daily values across months
    start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
    # only required if the layer name
    # format is not compatible with stagg
    time_interval = "1 hour", # The temporal interval of the supplied data,
    # only required if the start_date is not NA
    thresholds = c(0, 10, 20) # Calculate degree days above 0, 10, and 20
    # degrees Celsius
  )

  # Run new example
  degree_days_output_new <- staggregate_degree_days_new(
    data = terra::rast(temp_nj_jun_2024_era5) - 273.15, # Climate data to transform and
    # aggregate
    overlay_weights = overlay_weights_nj, # Output from overlay_weights()
    time_agg = "month", # Sum the transformed daily values across months
    start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
    # only required if the layer name
    # format is not compatible with stagg
    time_interval = "1 hour", # The temporal interval of the supplied data,
    # only required if the start_date is not NA
    thresholds = c(0, 10, 20) # Calculate degree days above 0, 10, and 20
    # degrees Celsius
  )

  expect_true(
    all.equal(degree_days_output, degree_days_output_new)
  )

})






