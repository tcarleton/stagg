# code and instructions to prepare `era5_grid` dataset

# Last download and run: October 4, 2024

# Read in data from
# https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download

# Options selected: reanalysis, 2m temperature, for Jan 1, 2024 for whole world
# as a NetCDF, renamed " ... .nc" to
# "ERA5_one_hour.nc"

# Dataset reference: Hersbach, H., Bell, B., Berrisford, P., Biavati, G.,
# Horányi, A., Muñoz Sabater, J., Nicolas, J., Peubey, C., Radu, R., Rozum, I.,
# Schepers, D., Simmons, A., Soci, C., Dee, D., Thépaut, J-N. (2023): ERA5
# hourly data on single levels from 1940 to present. Copernicus Climate Change
# Service (C3S) Climate Data Store (CDS), DOI: 10.24381/cds.adbb2d47

# API Query (Not used in actual download, but could be):
# import cdsapi
#
# dataset = "reanalysis-era5-single-levels"
# request = {
#   "product_type": ["reanalysis"],
#   "variable": ["2m_temperature"],
#   "year": ["2024"],
#   "month": ["01"],
#   "day": ["01"],
#   "time": ["00:00"],
#   "data_format": "netcdf",
#   "download_format": "unarchived"
# }
#
# client = cdsapi.Client()
# client.retrieve(dataset, request).download()


library(magrittr)

era5_grid <- raster::raster("data-raw/ERA5_one_hour.nc") %>%
  raster::readAll() %>%
  raster::setValues(values = 0) %>%
  terra::rast() %>%
  terra::rotate() %>%
  raster::raster()

names(era5_grid) <- "X"



# Use in package
usethis::use_data(era5_grid, overwrite = TRUE)


