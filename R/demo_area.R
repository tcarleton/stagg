#' A Weight Table Calculated Using Only SF Object
#'
#' A data table containing results from calc_geoweights() when a secondary
#' weight table is not included.
#'
#' @format A data table containing 42 observations of 5 variables:
#' \describe{
#'   \item{x}{longitude of cell based on 0 - 360 crs}
#'   \item{y}{latitude of cell}
#'   \item{poly_id}{unique identifyer for each county}
#'   \item{w_area}{proportion of the county falling within the cell}
#'   \item{w_sum}{sum of all area weights with that poly_id}
#' }
"demo_area_weights"
