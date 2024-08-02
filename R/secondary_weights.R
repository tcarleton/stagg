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
#'   processing. Format must be compatible with raster::crop()
#'
#' @return a data.table of secondary weights
#'
#' @examples
#' secondary_weights_output <- secondary_weights(
#'   secondary_raster = cropland_nj_2015, # A raster of cropland to resample
#'                                        # and convert to a data.table
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

  ## check if secondary raster fully overlaps with user-defined extent
  if(!is.character(extent)) {

    extent_rect <- raster::extent(extent)

    extent_secondary_r <- raster::extent(secondary_raster)

    # check if the raster extent covers the rectangle extent
    covers <- extent_secondary_r@xmin <= extent_rect@xmin &&
      extent_secondary_r@xmax >= extent_rect@xmax &&
      extent_secondary_r@ymin <= extent_rect@ymin &&
      extent_secondary_r@ymax >= extent_rect@ymax

    if (covers) {
      message(crayon::green('The secondary raster fully overlaps the user-specified extent.'))
    } else {
      message(crayon::yellow('Warning: the secondary raster does not fully overlap with the user-specified extent. Resulting data frame will not fully cover user-defined extent.'))
    }

  }

  # Create ERA raster from input raster
  clim_raster <- raster::raster(grid) # only reads the first band

  ## climate raster information for creating buffer and doing checks/rotations
  c_rast_xmax <- raster::extent(clim_raster)@xmax
  # c_rast_xmin <- raster::extent(clim_raster)@xmin
  c_rast_res <- raster::xres(clim_raster)

  ## add buffer to extent
  buffer_size <- c_rast_res

  ## add buffer to the extent for bbox
  if(is.vector(extent)) {
    extent <- c(
      extent[1] - buffer_size,  # xmin (negative direction)
      extent[2] + buffer_size,  # xmax (positive direction)
      extent[3] - buffer_size,  # ymin (negative direction)
      extent[4] + buffer_size   # ymax (positive direction)
    )
  } else if (inherits(extent, "bbox")) {
    extent <- c(
      extent$xmin - buffer_size,  # xmin (negative direction)
      extent$xmax + buffer_size,  # xmax (positive direction)
      extent$ymin - buffer_size,  # ymin (negative direction)
      extent$ymax + buffer_size   # ymax (positive direction)
    )
  } else if(inherits(extent, "Extent")) {
    extent <- c(
      extent@xmin - buffer_size,  # xmin (negative direction)
      extent@xmax + buffer_size,  # xmax (positive direction)
      extent@ymin - buffer_size,  # ymin (negative direction)
      extent@ymax + buffer_size   # ymax (positive direction)
    )
  } else if(is.character(extent)) {

      extent <- extent

    } else {

      stop("User-defined extent not compatible with raster.")

    }

  ## If an extent was included, crop it to the extent to save ram
  ## -----------------------------------------------
  if (!is.character(extent)){
    secondary_raster <- raster::crop(secondary_raster, extent)
  }

  ## Raster alignment: make sure clim_raster is in same coordinate system as secondary
  ## can be in either x coord system
  ## -----------------------------------------------

  message(crayon::yellow('Checking for raster alignment'))

  ## secondary raster
  s_rast_xmax <- raster::extent(secondary_raster)@xmax
  # s_rast_xmin <- raster::extent(secondary_raster)@xmin
  s_rast_res <- raster::xres(secondary_raster)


  ## if secondary raster in -180 to 180 and clim raster 0-360, rotate clim raster
  if(s_rast_xmax <= (180 + s_rast_res / 2) & c_rast_xmax >= (180 + c_rast_res / 2)) {

    message(crayon::yellow('Longitude coordinates do not match. Aligning longitudes to standard coordinates.'))

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

    message(crayon::yellow('Longitude coordinates do not match. Aligning longitudes to standard coordinates.'))

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
  ## use the extent of the previously user-cropped secondary raster
  ## -----------------------------------------------

  clim_raster <- raster::crop(clim_raster, raster::extent(secondary_raster))

  ## set crs of secondary raster to match climate data
  ## -----------------------------------------------
  raster::crs(secondary_raster) <- raster::crs(clim_raster)


  ## check if the cropped secondary raster contains NA values
  if(isTRUE(any(is.na(values(secondary_raster))))) {

    message(crayon::red("Warning: secondary raster contains NA values. NAs will be returned for weights."))

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
