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
