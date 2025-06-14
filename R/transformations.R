# 1. Polynomial Transformation Version of Staggregate
# ______________________________________________________________________________

#' Polynomial transformation and aggregation of climate data
#'
#' The function `staggregate_polynomial()` aggregates climate data to the daily
#' level, raises these daily values to the 1 through nth power, and aggregates
#' the transformed values to the polygon level and desired temporal scale.
#'
#' @param data The spatRasters with the data to be transformed and aggregated
#' @param overlay_weights A table of weights which can be generated using the
#'   function `overlay_weights()`
#' @param daily_agg How to aggregate hourly values to daily values prior to
#'   transformation. Options are `'sum'`, `'average'`, or `'none'` (`'none'`
#'   will transform values without first aggregating to the daily level)
#' @param time_agg the temporal scale to aggregate data to. Options are
#'   `'hour'`, `'day'`, `'month'`, or `'year'` (`'hour'` cannot be selected
#'   unless `daily_agg = 'none'`)
#' @param start_date the date (and time, if applicable) of the first layer in
#'  the stack. To be input in a format compatible with
#'  lubridate::as_datetime(), e.g. `"1991-10-29"` or `"1991-10-29 00:00:00"`.
#'  The default is `NA` since the spatRasters usually already contain temporal
#'  information in the layer names and they do not need to be manually supplied.
#' @param time_interval the time interval between layers in the spatRaster to be
#'  aggregated. To be input in a format compatible with seq(), e.g.
#'  `'1 day'` or `'3 months'`. The default is `'1 hour'` and this argument is
#'  required if daily_agg is not `'none'` or if the `start_date` argument is not
#'  `NA`.
#' @param weights_join_tolerance the tolerance to use when joining
#' overlay_weights with the climate data by the x and y columns. This is useful
#' when the height/width of your data cells expressed in degrees is a very long
#' decimal. The default, `0`, performs a keyed equi-join. Anything other than 0
#' performs a nonequi-join wherein latitudes/longitudes within the specified
#' tolerance (inclusive) are considered a match. Passing a single number results
#' in the tolerance being the same for x and y, but you can also pass a vector
#' of two numbers to have the first specify the x tolerance and second specify
#' the y tolerance.
#' @param na_rm whether to remove NAs in the climate data and adjust the
#' overlay weights so that the weights of non-na cells sum to 1 across each
#' polygon. The default is `FALSE`.
#' @param degree the highest exponent to raise the data to
#'
#' @examples
#' polynomial_output <- staggregate_polynomial_new(
#'   data = terra::rast(temp_nj_jun_2024_era5) - 273.15, # Climate data to transform and
#'                                          # aggregate
#'   overlay_weights = overlay_weights_nj, # Output from overlay_weights()
#'   daily_agg = "average", # Average hourly values to produce daily values
#'                          # before transformation
#'   time_agg = "month", # Sum the transformed daily values across months
#'   start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
#'                                       # only required if the layer name
#'                                       # format is not compatible with stagg
#'   time_interval = "1 hour", # The temporal interval of the supplied data,
#'                             # required if daily_agg is not "none" or if the
#'                             # start_date argument is not NA
#'   degree = 3 # Highest order
#'   )
#'
#' head(polynomial_output)
#'
#' @export
staggregate_polynomial_new <- function(
    data,
    overlay_weights,
    daily_agg = 'none',
    time_agg = 'month',
    degree,
    start_date = NA,
    time_interval = '1 hour',
    weights_join_tolerance = 0,
    na_rm = FALSE
  ){

  # Create transformations list of functions and result_cols list of column name
  # strings to be passed into staggregate_custom()
  transformations <- lapply(
    1:degree,
    function(exponent){force(exponent); \(x) x^exponent}
  )

  result_cols <- paste0('order_', 1:degree)


  # Run staggregate_custom()
  data <- staggregate_custom(
    data = data,
    overlay_weights = overlay_weights,
    daily_agg = daily_agg,
    time_agg = time_agg,
    start_date = start_date,
    time_interval = time_interval,
    weights_join_tolerance = weights_join_tolerance,
    na_rm = na_rm,
    transformations = transformations,
    result_cols = result_cols
  )

  return(data)

}
