# Generate the expected output and store as a static object to test against
# Last run: August 15, 2024

# Expected secondary_weights output
secondary_weights_full_extent_key <- secondary_weights(
  secondary_raster = cropland_nj_2015,
  grid = era5_grid,
  extent = "full"
)

saveRDS(object = secondary_weights_full_extent_key,
        file = testthat::test_path("fixtures/secondary_weights_full_extent_key.rds"))

secondary_weights_cropped_key <- secondary_weights(
  secondary_raster = cropland_nj_2015,
  grid = era5_grid,
  extent = c(-75.75, -73.5, 38.75, 41.25)
)

saveRDS(object = secondary_weights_cropped_key,
        file = testthat::test_path("fixtures/secondary_weights_cropped_key.rds"))

# Expected overlay_weights output
overlay_weights_cropland_key <- overlay_weights(
  polygons = nj_counties,
  polygon_id_col = "COUNTYFP",
  grid = era5_grid,
  secondary_weights = cropland_world_2015_era5
)

saveRDS(object = overlay_weights_cropland_key,
        file = testthat::test_path("fixtures/overlay_weights_cropland_key.rds"))

overlay_weights_area_key <- overlay_weights(
  polygons = nj_counties,
  polygon_id_col = "COUNTYFP",
  grid = era5_grid
)

saveRDS(object = overlay_weights_area_key,
        file = testthat::test_path("fixtures/overlay_weights_area_key.rds"))



# Generate staggregate_* keys
staggregate_polynomial_key <- staggregate_polynomial(
  data = temp_nj_jun_2024_era5 - 273.15, # Climate data to transform and aggregate
  overlay_weights = overlay_weights_nj, # Output from overlay_weights()
  daily_agg = "average", # Average hourly values to produce daily values
  # before transformation
  time_agg = "month", # Sum the transformed daily values across months
  start_date = "2020-01-01 00:00:00", # The start date of the supplied data, only required if the layer name format is not compatible with stagg
  time_interval = "1 hour", # The temporal interval of the supplied data, required if daily_agg is not "none" or if the start_date argument is not NA
  degree = 4 # Highest order
)

saveRDS(object = staggregate_polynomial_key,
        file = testthat::test_path("fixtures/staggregate_polynomial_key.rds"))



staggregate_spline_key <- staggregate_spline(
  data = temp_nj_jun_2024_era5 - 273.15, # Climate data to transform and aggregate
  overlay_weights = overlay_weights_nj, # Output from overlay_weights()
  daily_agg = "average", # Average hourly values to produce daily values before
  # transformation
  time_agg = "month", # Sum the transformed daily values across months
  start_date = "2020-01-01 00:00:00", # The start date of the supplied data, only required if the layer name format is not compatible with stagg
  time_interval = "1 hour", # The temporal interval of the supplied data, required if daily_agg is not "none" or if the start_date argument is not NA
  knot_locs = c(0, 7.5, 12.5, 20) # Where to place knots
)



saveRDS(object = staggregate_spline_key,
        file = testthat::test_path("fixtures/staggregate_spline_key.rds"))


staggregate_bin_key <- staggregate_bin(
  data = temp_nj_jun_2024_era5 - 273.15, # Climate data to transform and aggregate
  overlay_weights = overlay_weights_nj, # Output from overlay_weights()
  daily_agg = "average", # Average hourly values to produce daily values
  # before transformation
  time_agg = "month", # Sum the transformed daily values across months
  start_date = "2020-01-01 00:00:00", # The start date of the supplied data, only required if the layer name format is not compatible with stagg
  time_interval = "1 hour", # The temporal interval of the supplied data, required if daily_agg is not "none" or if the start_date argument is not NA
  bin_breaks = c(0, 2.5, 5, 7.5, 10) # Draw 6 bins from ninf to 0, 0 to 2.5,
  # 2.5 to 5, 5 to 7.5, 7.5 to 10, 10 to inf
)


saveRDS(object = staggregate_bin_key,
        file = testthat::test_path("fixtures/staggregate_bin_key.rds"))



staggregate_degree_days_key <- staggregate_degree_days(
  data = temp_nj_jun_2024_era5 - 273.15, # Climate data to transform and aggregate
  overlay_weights = overlay_weights_nj, # Output from overlay_weights()
  time_agg = "month", # Sum the transformed daily values across months
  start_date = "2020-01-01 00:00:00", # The start date of the supplied data, only required if the layer name format is not compatible with stagg
  time_interval = "1 hour", # The temporal interval of the supplied data, only required if the start_date is not NA
  thresholds = c(0, 10, 20) # Calculate degree days above 0, 10, and 20
  # degrees Celsius
)

saveRDS(object = staggregate_degree_days_key,
        file = testthat::test_path("fixtures/staggregate_degree_days_key.rds"))



as_data_table_raster_key <- as.data.table.raster(temp_nj_jun_2024_era5[[1]],
                                                 xy = TRUE)

saveRDS(object = as_data_table_raster_key,
        file = testthat::test_path("fixtures/as_data_table_raster_key.rds"))

