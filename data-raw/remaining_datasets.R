# This file documents how we obtained the user level datasets which currently do
# not have their raw forms in data-raw, nor the code used to generate them. This
# was written October 4, 2024

# ==============================================================================
# cropland_world_2015_era5
# ==============================================================================

# Time of creation:
# ------------------
# 2022 (I believe this was done with the early pipeline code, and that my code
# below are to give the idea of what was done but may not be exact. That's also
# the case for the population data, but I think that one was done a little
# later. See
# github.com/tcarleton/climate-aggregate/code/0.5_resample_rasters.R).



# Reason not in data-raw:
# ------------------------
# The raw data file is too large to be stored on git and very difficult to run
# locally. We should write a script to download directly if possible. Even then,
# such a script needs a lot of memory and would likely require remote computing.

# Steps to create:
# -----------------

# 1: Downloaded high resolution (not 3km) global cropland data from
#    https://glad.umd.edu/dataset/croplands in the four quadrants the site
#    provides

# 2: For each quadrant, ran contemporary version of stagg::secondary_weights()
#    with following code where 'downloaded_file' represents one of the four
#.   data quadrants downloaded in step 1:

     # quadrant_XX <- secondary_weights(
     #   secondary_raster = raster::raster('downloaded_file'),
     #   grid = era5_grid,
     #   extent = full
     # )

# 3: Row bound the results from step 2 and averaged any overlap by running

     # cropland_world_2015_era5 <- rbind(
     #   quadrant_NW,
     #   quadrant_NE,
     #   quadrant_SW,
     #  quadrant_SE
     # )
     # keys <- colnames(cropland_world_2015_era5)[!grepl('weight',colnames(full_table))]
     # full_table <- full_table[,list(weight= mean(weight)),keys]

# 4: Saved output from step 2 as package data by running

     # usethis::use_data(cropland_world_2015_era5, overwrite = TRUE)

# ==============================================================================
# pop_world_2015_era5
# ==============================================================================

# Time of creation:
# ------------------
# 2022 (see comment in cropland time section above)

# Reason not in data-raw:
# ------------------------
# The raw data file is too large to be stored on git. We should write a script
# to query the API if possible.

# Steps to create:
# -----------------

# 1:  Downloaded high resolution population data from
#     https://landscan.ornl.gov/

# 2:  Ran contemporary version of stagg::secondary_weights() [pre terra
#     conversion] with following code where 'downloaded_file' represents data
#     downloaded in step 1:

# pop_world_2015_era5 <- secondary_weights(
#   secondary_raster = raster::raster('downloaded_file'),
#   grid = era5_grid,
#   extent = full
# )

# 3: Saved output from step 2 as package data by running
# usethis::use_data(pop_world_2015_era5, overwrite = TRUE)

