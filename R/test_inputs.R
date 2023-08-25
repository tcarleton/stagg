## Tracey Mangin
## August 18, 2023
## inputs for testing package updates

##
library(raster)
library(sp)
library(rgdal)
library(tiff)

## path for crop rasters
tracey_laptop = "/Users/traceymangin/Library/CloudStorage/GoogleDrive-tmangin@ucsb.edu/Shared\ drives/emlab/projects/current-projects/climate-data-pipeline/demo-data/"
tracey_desktop = "/Users/tracey/Library/CloudStorage/GoogleDrive-tmangin@ucsb.edu/Shared\ drives/emlab/projects/current-projects/climate-data-pipeline/demo-data/"

## not working
# ## read in cropland demo data
# crop_ne <- raster(paste0(tracey_laptop, "era5_cropland_NE_2019_full.tif"))
# crop_ne <- readTIFF(paste0(tracey_laptop, "era5_cropland_NE_2019_full.tif"))
# crop_nw <- raster(paste0(tracey_laptop, "era5_cropland_NW_2019_full.tif"))

## inputs
kansas_counties <- tigris::counties("Kansas") ## polygons

polygons <- kansas_counties
polygon_id_col <- "COUNTYFP"
grid <- era5_grid
secondary_raster_ks <- cropland_kansas_2011
secondary_weights_world <- cropland_world_2015_era5
secondary_raster_world <- rasterFromXYZ(secondary_weights_world)


small_extent <- c(-95.75, -95.25, 37.25, 37.75)
small_extent_360 <- c(264.25, 264.75, 37.25, 37.75)
extent <- small_extent
pos_grid_extent <- c(-0.125, 179.875, -90.125, 90.125)
neg_grid_extent <- c(180.125, 359.875, -90.125, 90.125)

## kansas era 5
grid_rot <- rotate(era5_grid)

kansas_grid <- raster::crop(grid_rot, st_bbox(kansas_counties))
kansas_grid_pos <- raster::crop(grid, small_extent_360)
pos_grid <- raster::crop(grid, pos_grid_extent)
neg_grid <- raster::crop(grid, neg_grid_extent)

## test 1 secondary weights function:
## 0-360 climate data provided, world cropland, full extent
secondary_raster <- cropland_kansas_2011
grid <- neg_grid
extent <- "full"
# secondary_weights <-



