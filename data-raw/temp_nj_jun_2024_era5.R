# code and instructions to prepare `temp_nj_jun_2024_era5` dataset

# Last download: August 14, 2024
# Last code edit and run: October 4, 2024

# Read in data from
# https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download

# Options selected: reanalysis, 2m temperature, for all days and hours in June
# 2024 for a subregion from 38 S to 42.25 N and -76.5 W to -73 E, as a NetCDF,
# renamed "adaptor.mars.internal- ... .nc" to "NJ_temperature_jun_2024.nc"

# Dataset reference: Hersbach, H., Bell, B., Berrisford, P., Biavati, G.,
# Horányi, A., Muñoz Sabater, J., Nicolas, J., Peubey, C., Radu, R., Rozum, I.,
# Schepers, D., Simmons, A., Soci, C., Dee, D., Thépaut, J-N. (2023): ERA5
# hourly data on single levels from 1940 to present. Copernicus Climate Change
# Service (C3S) Climate Data Store (CDS), DOI: 10.24381/cds.adbb2d47

# API Query (Not used in actual download, but could be):
# import cdsapi
#
# c = cdsapi.Client()
#
# c.retrieve(
#   'reanalysis-era5-single-levels',
#   {
#     'product_type': 'reanalysis',
#     'variable': '2m_temperature',
#     'year': '2024',
#     'month': '06',
#     'day': [
#       '01', '02', '03',
#       '04', '05', '06',
#       '07', '08', '09',
#       '10', '11', '12',
#       '13', '14', '15',
#       '16', '17', '18',
#       '19', '20', '21',
#       '22', '23', '24',
#       '25', '26', '27',
#       '28', '29', '30',
#     ],
#     'time': [
#       '00:00', '01:00', '02:00',
#       '03:00', '04:00', '05:00',
#       '06:00', '07:00', '08:00',
#       '09:00', '10:00', '11:00',
#       '12:00', '13:00', '14:00',
#       '15:00', '16:00', '17:00',
#       '18:00', '19:00', '20:00',
#       '21:00', '22:00', '23:00',
#     ],
#     'area': [
#       42.25, -76.5, 38,
#       -73,
#     ],
#     'format': 'netcdf',
#   },
#   'download.nc')
library(magrittr)

temp_nj_jun_2024_era5 <- raster::brick("data-raw/NJ_temperature_jun_2024.nc") %>%
  raster::readAll()

# Use in package
usethis::use_data(temp_nj_jun_2024_era5, overwrite = TRUE)
