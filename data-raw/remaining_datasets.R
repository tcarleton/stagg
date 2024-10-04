# This file documents how we obtained the user level datasets which currently do
# not have their raw forms in data-raw, nor the code used to generate them. This
# was written October 4, 2024

# ==============================================================================
# cropland_world_2015_era5
# ==============================================================================

# Time of creation:
# ------------------
# 2022 (Could have been updated since then, I do not know)

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

# 3: Row bound the results from step 2 by running

     # cropland_world_205_era5 <- rbind(
     #   quadrant_NW,
     #   quadrant_NE,
     #   quadrant_SW,
     #  quadrant_SE
     # )

# 4: Saved output from step 2 as package data by running

     # usethis::use_data(cropland_world_2015_era5, overwrite = TRUE)

# ==============================================================================
# pop_world_2015_era5
# ==============================================================================

# Time of creation:
# ------------------
# August 12, 2024

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

