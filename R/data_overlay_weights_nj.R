#' Area and cropland weights in ERA5 grid for New Jersey counties
#'
#' A table of weights returned by running overlay_weights() with nj_counties and
#' secondary cropland weights.
#'
#' @format a data.table with 141 observations of 5 variables
#'   \describe{
#'     \item{x}{longitude, ranging from 0 to 360}
#'     \item{y}{latitude, ranging from -90 to 90}
#'     \item{polygon_id}{unique identifiers for counties in New Jersey}
#'     \item{w_area}{the proportion of the polygon that falls within a cell}
#'     \item{weight}{value to weight by in aggregation, determined by
#'                   multiplying the w_area value of a cell by its secondary
#'                   weight}
#'   }
#'
"overlay_weights_nj"
