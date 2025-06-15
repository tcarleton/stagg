#' Approximate New Jersey state boundaries
#'
#' A simple features object derived from running tigris::state(), limiting to
#' New Jersey, and running sf::st_simplify to reduce file size. This dataset is
#' what gets substituted for tigris::counties('NJ') is unavailable due to
#' internet.
#'
#' @format a data.table with 141 observations of 5 variables
#'   \describe{
#'     \item{GEOID}{Unique identifier of state of New Jersey}
#'     \item{geometry}{set of simplified geometrey points}
#'   }
#'
"nj_polygon"
