## tracey mangin
## october 2, 2023
## inputs for testing

## secondary weights testsing inputs
small_extent <- c(-95.75, -95.25, 37.25, 37.75)

secondary_raster <- cropland_kansas_2011
grid <- era5_grid
extent <- small_extent

## staggregate polynomial
data <- temp_kansas_jan_2020_era5
overlay_weights <- overlay_weights_kansas
daily_agg <- "average"
time_agg <-  "month"
degree <- 3

## overlay weights

cropland_weights <- dplyr::filter(cropland_world_2015_era5,
                                  x >= -103, x <= -94, y >= 37, y <= 41)

polygons <- tigris::counties("Kansas")
polygon_id_col <- "COUNTYFP"
grid <- era5_grid
secondary_weights <- cropland_weights


test_out <- w_norm %>%
  dplyr::rename(new_w_area = w_area,
         new_weight = weight) %>%
  dplyr::left_join(overlay_weights_kansas) %>%
  dplyr::mutate(diff_w_area = new_w_area - w_area,
         diff_weight = new_weight - weight)
