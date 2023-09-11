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
  ## -----------------------------------------------

  message(crayon::yellow('Checking for raster alignment'))

  ## secondary raster
  s_rast_xmax <- raster::extent(secondary_raster)@xmax
  s_rast_res <- raster::xres(secondary_raster)

  s_rast_coord <- if(s_rast_xmax > 180 + s_rast_res / 2) {
    "climate (0-360)"
  } else {
    "standard (-180 to 180)"}

  message(crayon::yellow('Secondary raster coordinate system:', s_rast_coord))

  ## climate raster
  c_rast_xmax <- raster::extent(clim_raster)@xmax
  c_rast_res <- raster::xres(clim_raster)

  c_rast_coord <- if(c_rast_xmax > 180 + c_rast_res / 2) {
    "climate (0-360)"
  } else {
    "standard (-180 to 180)"}

  message(crayon::yellow('Grid coordinate system:', c_rast_coord))

  ## check if coordinate systems match, if no shift raster in 0-360 format
  if(s_rast_coord != c_rast_coord) {

    message(crayon::yellow('Coordinate systems do not match. Adjusting raster longitude
                           to standard coordiantes between -180 and 180 degrees.'))

    ## create global extent for padding so rotate function can be used
    global_extent <- c(0, 360, -90, 90)

    ## shift the raster that has x coordinates in 0-360 format
    if(s_rast_coord == "standard (-180 to 180)") {

      ## check if raster needs to be padded, extend if needed
      c_rast_xmin <- raster::extent(clim_raster)@xmin

      if(!dplyr::near(c_rast_xmin, 0, tol = c_rast_res) | !dplyr::near(c_rast_xmax, 360, tol = c_rast_res)) {

        clim_raster <- raster::extend(clim_raster, global_extent)

      }

      ## shift raster
      message(crayon::yellow('Adjusting grid longitude from 0 to 360 to
                             standard coordinates between -180 and 180 degrees.'))

      clim_raster <- raster::rotate(clim_raster)

      } else {

        ## check if raster needs to be padded, extend if needed
        s_rast_xmin <- raster::extent(secondary_raster)@xmin

        if(!dplyr::near(s_rast_xmin, 0, tol = s_rast_res) | !dplyr::near(s_rast_xmax, 360, tol = s_rast_res)) {

          secondary_raster <- raster::extend(secondary_raster, global_extent)

        }

        ## shift raster
        message(crayon::yellow('Adjusting secondary raster longitude from 0 to 360 to
                               standard coordinates between -180 and 180 degrees.'))

        secondary_raster <- raster::rotate(secondary_raster)

      }


  } else {

    message(crayon::green('Coordinate systems match.'))

    }


  ## crop the ERA/climate raster to the polygon or at least the raster extent
  ## -----------------------------------------------
  if (!is.character(extent)){
    clim_raster <- raster::crop(clim_raster, extent)
  } else {
    clim_raster <- raster::crop(clim_raster, raster::extent(secondary_raster))
  }

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
