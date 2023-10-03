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

  ## make sure climate raster is in 0-360
  ## climate raster extent
  c_rast_xmin <- raster::extent(clim_raster)@xmin
  c_rast_res <- raster::xres(clim_raster)

  ## if climate raster is not in 0-360, shift
  clim_raster <- if(c_rast_xmin >= 0 - c_rast_res / 2) {
    clim_raster
  } else {
    raster::shift(clim_raster, dx = 360)}

  ## Raster alignment: make sure secondary raster is in the same coordinates system as climate (0-360)
  ## -----------------------------------------------

  message(crayon::yellow('Checking for raster alignment'))

  ## secondary raster
  s_rast_xmin <- raster::extent(secondary_raster)@xmin
  s_rast_res <- raster::xres(secondary_raster)

  ## check if secondary coordinate system is in 0-360, if no shift raster in 0-360 format
  if(s_rast_xmin < 0 - s_rast_res / 2) {

    message(crayon::yellow('Aligning raster longitudes to 0-360 coordinates.'))

    secondary_raster <- raster::shift(secondary_raster, dx = 360)

  }

  ## crop the ERA/climate raster to the appropriate extent
  ## secondary raster was previously cropped according to user input, so use that
  ## -----------------------------------------------

  clim_raster <- raster::crop(clim_raster, raster::extent(secondary_raster))

  ## set crs of secondary raster to match climate data
  ## -----------------------------------------------
  raster::crs(secondary_raster) <- raster::crs(clim_raster)


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
