#' Resample a raster of secondary weights
#'
#' The `secondary_weights()` function resamples a raster of interest to the
#' desired resolution and outputs a table of weights for use in the function
#' `overlay_weights()`.
#'
#' @param secondary_raster a raster of a secondary variable, for example
#'   cropland coverage or population
#' @param grid a raster layer with the same spatial resolution as the climate
#'   data
#' @param extent an optional extent to crop the secondary_raster to for faster
#'   processing. Takes a vector in the form `c(min_lon, max_lon, min_lat,
#'   max_lat)`
#'
#' @return a data.table of secondary weights
#'
#' @examples
#' secondary_weights_output <- secondary_weights(
#'   secondary_raster = cropland_kansas_2011, # A raster of cropland to resample
#'                                            # and convert to a data.table
#'   grid = era5_grid, # The grid to resample the secondary_raster to (`era5_grid`
#'                     # is the default)
#'   extent = "full" # The default which resamples the whole secondary raster
#'                   # without cropping (`'full'` is the default)
#'   )
#'
#' head(secondary_weights_output)
#'
#' @export
secondary_weights <- function(secondary_raster, grid = era5_grid, extent = "full"){

  ## If an extent was included, crop it to the extent to save ram
  ## -----------------------------------------------
  if (!is.character(extent)){
    secondary_raster <- raster::crop(secondary_raster, extent)
  }


  # Create ERA raster from input raster
  clim_raster <- raster::raster(grid) # only reads the first band

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
    clim_raster <- raster::crop(clim_raster, raster::extent(secondary_raster))
  }


  ## Make the values of the clim_raster resampled weights
  ## -----------------------------------------------
  message(crayon::green("Resampling secondary_raster"))

  resampled_raster = raster::resample(secondary_raster, clim_raster, method="bilinear")



  ## Make a data.table of the values of the resampled raster with lat/lon
  ## -----------------------------------------------
  message(crayon::green("Creating a table of weights"))

  weight_table <- raster::as.data.frame(resampled_raster, xy=TRUE)
  colnames(weight_table) <- c("x", "y", "weight")
  weight_table <- data.table::as.data.table(weight_table)


  ## Return the weights table
  ## -----------------------------------------------
  return(weight_table)

}
