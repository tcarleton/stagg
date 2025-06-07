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
check_alignment <- function(data, x_res){

  # Get max longitudinal extent
  if(inherits(data, "SpatRaster")){
    max_x <- terra::ext(data)$xmax
  }else if(inherits(data, 'data.table')){
    # Go out half a cell width since this is a table of centroids
    max_x <- max(data[,x]) + (x_res / 2)
  }else{
    stop(crayon::red("Acceptable data type not found for function check_alignment()"))
  }

  # Definition of climate coords: data has cell which is entirely "right" of 180
  if(max_x >= 180 + x_res){
    alignment <- 'climate'
  } else{
    alignment <- 'standard'
  }

  return(alignment)

}
