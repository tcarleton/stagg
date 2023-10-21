## tracey mangin
## october 2, 2023
## inputs for testing

## secondary weights testsing inputs
small_extent <- c(-95.75, -95.25, 37.25, 37.75)

secondary_raster <- cropland_kansas_2011
grid <- era5_grid
extent <- small_extent

## overlay weights
cropland_weights <- dplyr::filter(cropland_world_2015_era5,
                                  x >= -103, x <= -94, y >= 37, y <= 41)

polygons <- tigris::counties("Kansas")
polygon_id_col <- "COUNTYFP"
grid <- era5_grid
secondary_weights <- cropland_weights

## staggregate polynomial
data <- temp_kansas_jan_2020_era5
overlay_weights <- overlay_weights_kansas
daily_agg <- "average"
time_agg <-  "month"
degree <- 3




## jrising data
rr.tas <- raster::stack("R/extra-for-testing-delete/tas_day_BEST_historical_station_19800101-19891231.nc")
rr.tas <- terra::rast("R/extra-for-testing-delete/tas_day_BEST_historical_station_19800101-19891231.nc")

rr.pop <- raster::raster("R/extra-for-testing-delete/usap90ag.nc")

counties <- tigris::counties()
counties$FIPS <- paste0(counties$STATEFP, counties$COUNTYFP)

## secondary weights
secondary_raster <- rr.pop
grid <- rr.tas

## overlay weights
polygons <- counties
polygon_id_col <- "FIPS"
grid <- rr.tas
secondary_weights <- weight_table

## staggregate
data = rr.tas
overlay_weights = w_norm
daily_agg = "average"
time_agg = "month"
degree = 1


# grid.weights <- secondary_weights(secondary_raster=rr.pop, grid=rr.tas)
# county.weights <- overlay_weights(polygons=counties, polygon_id_col="FIPS", grid=rr.tas, secondary_weights=grid.weights)
#
# county.tas <- staggregate_polynomial(data=rr.tas, daily_agg="none", time_agg='year', overlay_weights=county.weights, degree=2)
## Error in .local(x, y, ...) : extents do not overlap
