# code to prepare `nj_polygon` dataset

# last run: June 15, 2025

# Pull tigris supplied state boundaries and simplify for size constraints

# ______________________________________________________________________________

# Load libraries
library(sf)
library(tigris)
library(dplyr)
library(ggplot2)

# Query census state boundaries from tigris
nj_full <- tigris::states() |>

  # Only save NJ
  filter(STUSPS == 'NJ') |>

  # Remove unecessary columns
  select(GEOID, geometry)


nj_simplified <- nj_full |>
  st_simplify(dTolerance = 5000)

# Plot to see difference
ggplot() +
  geom_sf(data = nj_simplified, fill = NA, color = "red") +
  geom_sf(data =nj_full, fill = NA, color = "blue")

# Change name to nj_polygon and save
nj_polygon <- nj_simplified

usethis::use_data(nj_polygon, overwrite = TRUE)
