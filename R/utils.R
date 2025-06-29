#' Check whether data is in climate or standard coordinates
#'
#' A helper function to evaluate whether the passed object, either a spatRast or
#' data.table, is in climate (0 to 360) or standard (-180 to 180) coordinate
#' alignment. This gives us a consistent definition throughout the package.
#'
#' @param data either a spatRast or a data.table with a column 'x' denoting
#' the longitudes of cell centroids
#' @param x_res the cell width
#'
#' @returns either 'climate' or 'standard'
#'
#' @noRd
check_alignment <- function(data, x_res = NA){

  # Get max longitudinal extent
  if(inherits(data, "SpatRaster")){
    max_x <- terra::ext(data)$xmax
    x_res <- terra::xres(data)
  }else if(inherits(data, 'data.table')){
    # Go out half a cell width since this is a table of centroids
    max_x <- max(data[,x]) + (x_res / 2)
  }else{
    stop(crayon::red("Acceptable data type not found for function check_alignment()"))
  }

  # Definition of climate coords: data has cell which is entirely "right" of 180
  if(max_x >= (180 + x_res)){
    alignment <- 'climate'
  } else{
    alignment <- 'standard'
  }

  return(alignment)

}




#' Calculate time steps per day
#'
#' Take a given time interval and calculate the number of intervals within 1 day
#'
#' @param time_interval a string compatible with lubridate::duration()
#'
#' @returns A numeric timesteps per day
#'
#' @noRd
calc_intervals_in_day <- function(time_interval){
  as.numeric(
    (lubridate::duration("1 day") / lubridate::duration(time_interval))
  )
}



#' Layer names to dates
#'
#' Wrapper function around lubridate::as_datetime() to remove the leading
#' non-digit characters
#'
#' @param layer_names vector of strings to coerce to datetime objects
#'
#' @returns A vector of datetimes. Error if not compatible with
#' lubridate::as_datetime()
#'
#' @noRd
layer_names_to_dates <- function(layer_names){


  layer_names |>

    # Remove any non-digit characters from the start of the string
    stringr::str_replace('[^0-9]+', '') |>

    # Coerce to datetimes
    lubridate::as_datetime()
}



