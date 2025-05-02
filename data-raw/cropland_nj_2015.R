# code to prepare `cropland_nj_2015` dataset

# last run: October 4, 2024

# Read in raw data, downloaded from https://glad.umd.edu/dataset/croplands and
# crop

# Dataset Reference: P. Potapov, S. Turubanova, M.C. Hansen, A. Tyukavina, V.
# Zalles, A. Khan, X.-P. Song, A. Pickens, Q. Shen, J. Cortez. (2021) Global
# maps of cropland extent and change show accelerated cropland expansion in the
# twenty-first century. Nature Food. https://doi.org/10.1038/s43016-021-00429-z

library(magrittr)

# Crop
cropland_nj_2015 <- raster::crop(
  raster::raster("data-raw/Global_cropland_3km_2015.tif"),
  c(-76.5, -73, 38.25, 42)) %>%
  raster::readAll()

# Turn percentage to decimal
cropland_nj_2015 <- cropland_nj_2015 / 100


# Use in package
usethis::use_data(cropland_nj_2015, overwrite = TRUE)
