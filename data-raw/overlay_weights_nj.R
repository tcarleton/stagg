# code to prepare `overlay_weights_nj` dataset

# last run: August 14, 2024

# load county polygons
nj_counties <- tigris::counties("nj")

# run overlay_weights on nj data
overlay_weights_nj <- overlay_weights(
  polygons = nj_counties,
  polygon_id_col = "COUNTYFP",
  grid = era5_grid,
  secondary_weights = cropland_world_2015_era5)

# Use in package
usethis::use_data(overlay_weights_nj, overwrite = TRUE)
