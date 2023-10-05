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


  ## Raster alignment: make sure clim_raster is in same coordinate system as secondary
  ## can be in either x coord system
  ## -----------------------------------------------

  message(crayon::yellow('Checking for raster alignment'))

  ## secondary raster
  s_rast_xmax <- raster::extent(secondary_raster)@xmax
  # s_rast_xmin <- raster::extent(secondary_raster)@xmin
  s_rast_res <- raster::xres(secondary_raster)

  ## climate raster
  c_rast_xmax <- raster::extent(clim_raster)@xmax
  # c_rast_xmin <- raster::extent(clim_raster)@xmin
  c_rast_res <- raster::xres(clim_raster)

  ## if secondary raster in -180 to 180 and clim raster 0-360, rotate clim raster
  if(s_rast_xmax <= (180 + s_rast_res / 2) & c_rast_xmax >= (180 + c_rast_res / 2)) {

    message(crayon::yellow('Longitude coordinates do not match. Aligning longitudes to -180 to 180.'))

    ## check if raster needs to be padded, extend if needed
    c_rast_xmin <- raster::extent(clim_raster)@xmin

    if(!dplyr::near(c_rast_xmin, 0, tol = c_rast_res) | !dplyr::near(c_rast_xmax, 360, tol = c_rast_res)) {

      ## create global extent for padding so rotate function can be used
      global_extent <- c(0, 360, -90, 90)

      ## pad
      clim_raster <- raster::extend(clim_raster, global_extent)

    }

    ## rotate
    clim_raster <- raster::rotate(clim_raster)

  }

  ## if secondary raster in 0-360 and clim raster -180 to 180, rotate secondary raster
  if(s_rast_xmax >= (180 + s_rast_res / 2) & c_rast_xmax <= (180 + c_rast_res / 2)) {

    message(crayon::yellow('Longitude coordinates do not match. Aligning longitudes to -180 to 180.'))

    ## check if raster needs to be padded, extend if needed
    s_rast_xmin <- raster::extent(secondary_raster)@xmin

    if(!dplyr::near(s_rast_xmin, 0, tol = s_rast_res) | !dplyr::near(s_rast_xmax, 360, tol = s_rast_res)) {

      ## create global extent for padding so rotate function can be used
      global_extent <- c(0, 360, -90, 90)

      ## pad
      secondary_raster <- raster::extend(secondary_raster, global_extent)

    }

    ## rotate
    secondary_raster <- raster::rotate(secondary_raster)

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
