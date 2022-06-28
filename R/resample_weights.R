## Anna Boser
## March 27, 2022
## Make a general function that resamples a raster of interest to the ERA grid
## Processing step for potapov crop weights to turn into a single global map by year

#' Function that resamples a raster of interest to the ERA grid
#'
#' @param weights_raster a raster of a continuous variable, for example cropland coverage or population
#' @param data_source the source of climate data (default is era5)
#' @param extent an optional extent to crop the weights_raster to for faster processing
#'
#' @return a data.table of geoweights (area weighted raster/polygon overlap)

# Package will not write to any data folder, and instead put objects into environment

#' @export
calc_raster_weights <- function(weights_raster, data_source, extent = "full"){

  ## Setup
  ## -----------------------------------------------
  # Removed load packages lines which will be installed by the imports section in description and called using the pkg::fun() syntax. See bottom of user_run_example for full list of packages to be included in imports.

  ## If an extent was included, crop it to the extent to save ram
  ## -----------------------------------------------
  if (!is.character(extent)){
    weights_raster <- raster::crop(weights_raster, extent)
  }


  # Create ERA raster from input raster
  clim_raster <- raster::raster(data_source) # only reads the first band

  ## Raster alignment: make sure clim_raster is -180 to 180 longitude
  ## -----------------------------------------------

  message(crayon::yellow('Checking for raster alignment'))

  poly_xmin <- -180
  poly_xmax <- 180
  rast_xmin <- raster::extent(clim_raster)@xmin
  rast_xmax <- raster::extent(clim_raster)@xmax

  # Rotate raster if initial longitudes don't align
  if(!dplyr::near(poly_xmax, rast_xmax, tol=1.01)) {

    message(crayon::yellow('Adjusting raster longitude from',
                           round(rast_xmin,0), '-', round(rast_xmax,0),
                           'to', round(poly_xmin,0), '-', round(poly_xmax,0)))

    clim_raster <- raster::rotate(clim_raster)

  }

  # Check longitude ranges match (with a tolerance of 1 in case lon +- 179 vs. +-180)
  poly_range <- c(-180, 180)
  rast_range <- c(raster::extent(clim_raster)@xmin, raster::extent(clim_raster)@xmax)

  if(dplyr::near(poly_range[1], rast_range[1], tol=1.01) & dplyr::near(poly_range[2], rast_range[2], tol=1.01)){

    message(crayon::green('Longitude ranges match'))

  } else {

    stop(crayon::red('Raster longitude must be -180 to 180'))

  }

  ## crop the ERA raster to the polygon or at least the raster extent
  ## -----------------------------------------------
  if (!is.character(extent)){
    clim_raster <- raster::crop(clim_raster, extent)
  } else {
    clim_raster <- raster::crop(clim_raster, raster::extent(weights_raster))
  }

  ## Match raster crs
  ## -----------------------------------------------
  # weights_raster <- weights_raster %>%
  #   st_transform(crs = st_crs(clim_raster))
  # this doesn't work and as far as I can tell this isn't super necessary

  ## Make the values of the clim_raster resampled weights
  ## -----------------------------------------------
  resampled_raster = raster::resample(weights_raster, clim_raster, method="bilinear")



  ## Make a data.table of the values of the resampled raster with lat/lon
  ## -----------------------------------------------
  weight_table <- raster::as.data.frame(resampled_raster, xy=TRUE)
  colnames(weight_table) <- c("x", "y", "weight")
  weight_table <- data.table::as.data.table(weight_table)


  ## Return the weights table
  ## -----------------------------------------------
  return(weight_table)

}
