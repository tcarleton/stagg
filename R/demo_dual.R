#' A Weight Table Calculated Using Secondary Weights
#'
#' A data table containing results from calc_geoweights() when a secondary
#' weight table is included.
#'
#' @format A data table containing 42 observations of 5 variables:
#' \describe{
#'   \item{x}{longitude of cell based on 0 - 360 crs}
#'   \item{y}{latitude of cell}
#'   \item{poly_id}{unique identifyer for each county}
#'   \item{w_area}{proportion of the county falling within the cell}
#'   \item{weight}{secondary weights normalized by area}
#' }
"demo_dual_weights"
