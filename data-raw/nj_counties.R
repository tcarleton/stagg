# code to prepare `nj_counties` dataset
nj_counties <- tigris::counties("nj")

# Use in package
usethis::use_data(nj_counties, overwrite = TRUE)
