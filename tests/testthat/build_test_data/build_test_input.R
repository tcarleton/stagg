# Build Test / Example Data Inputs

# Last Run: August 13, 2024

# This is the code that is used / process that is followed to construct the data
# that appears in the tests and examples. It generates the following package
# data:
  # cropland_nj_2015.rda
  # overlay_weights_nj.rda
  # temp_nj_jun_2024_era5.rda

# Code run:
build_test_input(
  "~/Desktop/Global_cropland_NW_2015.tif",
  "~/Desktop/adaptor.mars.internal-1723526351.1015518-19422-3-13be0fd0-ae10-448c-a278-1c8a8cae14cc.nc"
  )



build_test_input <- function(

  # String specifying local location of Global_cropland_2015.tif (outside of
  # stagg directory)
  full_cropland_raster_path,

  # String specifying local location of ERA5 temperature stack (outside of
  # stagg directory)
  nj_era5_temperature_stack_path

  ){

  # ======================================
  # Step 1: Generate cropland_nj_2015.rda
  # ======================================

  # (A) [execute manually] - download cropland data and pass file path to
  # `full_cropland_raster_path` argument

  # (B) - read in and crop
  cropland_nj_2015 <- raster::crop(
    raster::raster(full_cropland_raster_path),
    c(-76.5, -73, 38, 42)
  )

  # (C) - create a blank grid with the proper resolution
  blank_grid <- raster::aggregate(
    raster::setValues(cropland_nj_2015, values = 0),
    100
  )

  # (D) - resample the cropland data to the resolution of the blank grid
  cropland_nj_2015 <- raster::resample(
    cropland_nj_2015,
    blank_grid
  )


  # (E) - save with usethis::use_data() and document
  usethis::use_data(cropland_nj_2015, overwrite = TRUE)


  # ========================================
  # Step 2: Generate overlay_weights_nj.rda
  # ========================================


  # ===========================================
  # Step 3: Generate temp_nj_jun_2024_era5.rda
  # ===========================================

  # (A) [execute manually] - download ERA5 raster brick
  # Visit
  # https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form
  # and download reanalysis, 2m temperature, for all days and hours in June 2024
  # for a subregion from 38 S to 42.25 N and -76.5 W to -73 E, as a NetCDF

  # (B) - save with usethis::use_data
}




