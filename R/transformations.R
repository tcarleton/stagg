# Staggregate Transformation Presets
#
# This file contains several different wrapper functions for
# staggregate_custom(). They all implement a different transformation based on
# one or more unique arguments (staggregate_polynomial() takes degree, for
# instance) and then pass that transformation, and all of the arguments supplied
# by the user, to staggregate_custom() and output its result.

# For this reason, running
# staggregate_custom(..., transformations = c(function(x){x}, function(x){x^2}))
# is equivalent to running staggregate_polynomial(..., degree = 2) as long as
# the ... arguments are identical. All staggregate_polynomial is doing is taking
# the information supplied to degree = 2 and interpreting this as a list of two
# transformations, x and x^2. It then passes these to staggregate_custom() along
# with all the ... arguments the user supplied to it, and runs
# staggregate_custom "under the hood".

# ==============================================================================
# Polynomial Transformation Version of Staggregate
# ==============================================================================

# Internal Helper Functions
# ______________________________________________________________________________

#   a) validate_degree
#   -----------------------------------
#' Make Sure Degree is a natural number
#'
#' @param degree the user supplied degree
#'
#' @returns None. Error if degree is not natural number
#'
#' @noRd
validate_degree <- function(degree){
  throw <- any(
    !is.numeric(degree),
    !length(degree) == 1,
    !degree %% 1 == 0,
    !degree > 0
  )

  if(throw){
    stop(crayon::red('degree must be a natural number'))
  }

  return(degree)
}

# Exported Staggregate Function
# ______________________________________________________________________________

#' Polynomial transformation and aggregation of climate data
#'
#' The function `staggregate_polynomial()` aggregates climate data to the daily
#' level, raises these daily values to the 1 through nth power, and aggregates
#' the transformed values to the polygon level and desired temporal scale.
#'
#' @inheritParams staggregate_custom
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
    start_date = NA,
    time_interval = '1 hour',
    weights_join_tolerance = 0,
    na_rm = FALSE,
    degree
  ){

  degree <- validate_degree(degree)

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


# ==============================================================================
# Restricted Cubic Spline Transformation Version of Staggregate
# ==============================================================================

# Internal Helper Functions
# ______________________________________________________________________________

#   a) validate_knot_locs
#   -----------------------------------
#' Make sure knot_locs are in order
#'
#' @param knot_locs The user supplied knot_locs
#'
#' @returns the numeric vector of knot_locs in order from least to greatest
#'
#' @noRd
validate_knot_locs <- function(knot_locs){
  if(!is.numeric(knot_locs) | length(knot_locs < 3)){
    stop(crayon::red('knot_locs must be a numeric vector with length of at least 3'))
  }

  knot_locs <- sort(knot_locs)

  return(knot_locs)
}




#   b) u_pos
#   -----------------------------------
#' Terra compatible eval of  turn to zero if negative (helper for get_spline_eq)
#'
#' @param u value to evaluate and return 0 if negative
#'
#' @returns 0 if u <= 0, u otherwise
#'
#' @noRd
u_pos <- function(u){terra::ifel(u > 0, u, 0)}




#   c) get_spline_eq
#   -----------------------------------
#' Generate restricted cubic spline equation for element in knot_loc vector
#'
#' Function to create 1 f(x)=x function and k - 2 restricted cubic spline
#' functions for k = length(knot_locs) based on the equation on page 3 of
#' https://support.sas.com/resources/papers/proceedings16/5621-2016.pdf.
#'
#' @param t the vector of knot locations, should be passed user supplied
#'   `knot_locs`
#' @param i the index of t (knot locations) to pull. 0 <= i < length(t) - 2
#' @param x the value to run through each function. This should receive the
#'   actual climate data and be the only non-fixed variable in the
#'   transformations list
#'
#' @returns A spline function based on the values of i
#'
#' @noRd
get_spline_eq <- function(t, i, x){

  # Calc k
  k <- length(t)

  # Have the function return x if i == 0, otherwise compute spline equations
  if(i == 0){
    output <- x
  } else{

    # x minus t of i cubed if positive
    part_1 <- u_pos((x - t[i])^3)

    # x minus t of k-1 ...
    part_2 <- u_pos((x - t[k-1])^3) * (t[k] - t[i]) / (t[k] - t[k-1])

    # x minus t of k ...
    part_3 <- u_pos((x - t[k])^3) * (t[k-1] - t[i]) / (t[k] - t[k-1])

    output <- part_1 - part_2 + part_3
  }

  return(output)
}


# Exported Main Function
# ______________________________________________________________________________

#' Restricted cubic spline transformation and aggregation of climate data
#'
#' The function `staggregate_spline()` aggregates climate data to the daily
#' level, performs a restricted cubic spline transformation on these daily
#' values, and aggregates the transformed values to the polygon level and
#' desired temporal scale.
#'
#' @inheritParams staggregate_custom
#'
#' @param knot_locs where to place the knots
#'
#' @examples
#' spline_output <- staggregate_spline(
#' data = terra::rast(temp_nj_jun_2024_era5) - 273.15, # Climate data to transform and
#'                                        # aggregate
#' overlay_weights = overlay_weights_nj, # Output from overlay_weights()
#' daily_agg = "average", # Average hourly values to produce daily values
#'                        # before transformation
#' time_agg = "month", # Sum the transformed daily values across months
#' start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
#'                                     # only required if the layer name format
#'                                     # is not compatible with stagg
#' time_interval = "1 hour", # The temporal interval of the supplied data,
#'                           # required if daily_agg is not "none" or if the
#'                           # start_date argument is not NA
#' knot_locs = c(0, 7.5, 12.5, 20) # Where to place knots
#' )
#'
#' head(spline_output)
#'
#' @export
staggregate_spline_new <- function(
    data,
    overlay_weights,
    daily_agg,
    time_agg = "month",
    start_date = NA,
    time_interval = "1 hour",
    weights_join_tolerance = 0,
    na_rm = FALSE,
    knot_locs
  ){

  # Create 1 function that's just f(x) = x, and then length(knot_locs)-2 spline functions
  transformations <- lapply(
    0:(length(knot_locs) - 2),
    function(knot_locs_index){
      force(knot_locs)
      force(knot_locs_index)
      \(x) get_spline_eq(
        t = knot_locs,
        i = knot_locs_index,
        x = x
      )
    }
  )

  # Create spline column output names
  result_cols <- paste0('spline_term_', 1:(length(knot_locs) - 2))

  # Append original value name to start of result_cols
  result_cols <- c('orig_value', result_cols)

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

# ==============================================================================
# Bin Transformation Version of Staggregate
# ==============================================================================

# Internal Helper Functions
# ______________________________________________________________________________

#   a) validate_bin_breaks
#   -----------------------------------
#' Make sure bin_breaks are in order
#'
#' @param bin_breaks The user supplied bin_breaks
#'
#' @returns the numeric vector of bin_breaks in order from least to greatest
#'
#' @noRd
validate_bin_breaks <- function(bin_breaks){
  if(!is.numeric(bin_breaks)){
    stop(crayon::red('bin_breaks must be a numeric vector'))
  }

  bin_breaks <- sort(bin_breaks)

  return(bin_breaks)
}

#   b) get_bin_funs
#   -----------------------------------
#' Create a list of binning functions from `bin_breaks`
#'
#' @param bin_breaks the ascending numeric vector of bin_breaks
#' @param i the index 0 to length(bin_breaks) to iterate over (corresponds to
#'   bin number)
#' @param x the value to run through each function. This should receive the
#'   actual climate data and be the only non-fixed variable in the
#'   transformations list
#'
#' @returns a list of functions to pass to `transformations`
#'
#' @noRd
get_bin_funs <- function(bin_breaks, i, x){

  if(i == 0){

    # Minimum bin
    output <- terra::ifel(x < min(bin_breaks), 1, 0)

  } else if(i == length(bin_breaks)){

    # Maximum bin
    output <- terra::ifel(x >= max(bin_breaks), 1, 0)

  } else{

    # All other bins
    output <- terra::ifel((bin_breaks[i] <= x) & (x < bin_breaks[i+1]), 1, 0)
  }

  return(output)

}


#   c) get_bin_names
#   -----------------------------------
#' Assign bin column names
#'
#' @param bin_breaks the ascending numeric vector of bin_breaks
#'
#' @returns the bin names to be supplied to result_cols
#'
#' @noRd
get_bin_names <- function(bin_breaks){

  # Assign bin column names
  for(i in 0:length(bin_breaks)){
    if(i == 0){
      result_cols <- paste0('bin_ninf_to_', min(bin_breaks))
    } else if(i == length(bin_breaks)){
      result_cols <- c(
        result_cols,
        paste0('bin_', max(bin_breaks), '_to_inf')
      )
    } else{
      result_cols <- c(
        result_cols,
        paste0('bin_', bin_breaks[i], '_to_', bin_breaks[i+1])
      )
    }
  }

  result_cols <- sub('-', 'n', result_cols)

  return(result_cols)
}

# Exported Main Function
# ______________________________________________________________________________
#
#' Bin transformation and aggregation of climate data
#'
#' The function `staggregate_bin()` aggregates climate data to the daily level,
#' splits these daily values into bins, and aggregates the transformed
#' values to the polygon level and desired temporal scale.
#'
#' @inheritParams staggregate_custom
#'
#' @param bin_breaks A vector of bin boundaries to split the data by
#'
#' @examples
#' bin_output <- staggregate_bin(
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
#'   bin_breaks = c(0, 2.5, 5, 7.5, 10) # Draw 6 bins from ninf to 0, 0 to 2.5,
#'                                      # 2.5 to 5, 5 to 7.5, 7.5 to 10, 10 to
#'                                      # inf
#'   )
#'
#' head(bin_output)
#'
#' @export
staggregate_bin_new <- function(
    data,
    overlay_weights,
    daily_agg,
    time_agg = "month",
    start_date = NA,
    time_interval = '1 hour',
    weights_join_tolerance = 0,
    na_rm = FALSE,
    bin_breaks
){

  # Make sure bin_breaks are ordered vector
  bin_breaks <- validate_bin_breaks(bin_breaks)

  # Create list of functions corresponding to each bin
  transformations <- lapply(
    0:length(bin_breaks),
    function(bin_index){
      force(bin_breaks)
      force(bin_index)
      \(x) get_bin_funs(
        bin_breaks = bin_breaks,
        i = bin_index,
        x = x
      )
    }
  )

  # Assign bin names
  result_cols <- get_bin_names(bin_breaks)

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




# ==============================================================================
# Degree Days Transformation Version of Staggregate
# ==============================================================================

# Internal Helper Functions
# ______________________________________________________________________________

#   a) validate_thresholds
#   -----------------------------------
#' Make sure thresholds are in order
#'
#' @param thresholds The user supplied thresholds
#'
#' @returns the numeric vector of thresholds in order from least to greatest
#'
#' @noRd
validate_thresholds <- function(thresholds){
  if(!is.numeric(thresholds)){
    stop(crayon::red('thresholds must be a numeric vector'))
  }

  thresholds <- sort(thresholds)

  return(thresholds)
}




#   b) get_deg_days_fun
#   -----------------------------------
#' Create a list degree days functions from `thresholds`
#'
#' @param thresholds the ascending numeric vector of thresholds
#' @param i the index 0 to length(thresholds) to iterate over (corresponds to
#'   threshold number)
#' @param x the value to run through each function. This should receive the
#'   actual climate data and be the only non-fixed variable in the
#'   transformations list
#'
#' @returns a list of functions to pass to `transformations`
#'
#' @noRd
get_deg_day_funs <- function(thresholds, i, x){

  if(i == 0){

    # Lowest threshold, threshold - x if x < threshold, 0 otherwise,
    output <- terra::ifel(
      x < min(thresholds),
      min(thresholds) - x,
      0
    )

  } else if(i == length(bin_breaks)){

    # Highest threshold, x - threshold if x > threshold, 0 otherwise
    output <- terra::ifel(
      x > max(thresholds),
      x - max(thresholds),
      0
    )

  } else{
    # Can't use case_when here because it doesn't work on terra objects.
    # Potential suggestion for terra team

    # All other thresholds:
    # 0 if x < threshold,
    # next_threshhold - threshold if x > next_threshold
    # x - threshold otherwise
    output <- terra::ifel(

      # x below range
      x < thresholds[i],
      0,

      # x above range
      terra::ifel(
        x > thresholds[i + 1],
        thresholds[i + 1] - thresholds[i],

        # x in range
        x - thresholds[i]
      )
    )
  }

  return(output)

}


#   c) get_deg_day_names
#   -----------------------------------
#' Assign bin column names
#'
#' @param thresholds the ascending numeric vector of bin_breaks
#'
#' @returns the bin names to be supplied to result_cols
#'
#' @noRd
get_deg_day_names <- function(thresholds){

  # Assign bin column names
  for(i in 0:length(thresholds)){
    if(i == 0){
      result_cols <- paste0('threshold_ninf_to_', min(thresholds))
    } else if(i == length(bin_breaks)){
      result_cols <- c(
        result_cols,
        paste0('threshold_', max(thresholds), '_to_inf')
      )
    } else{
      result_cols <- c(
        result_cols,
        paste0('threshold_', bin_breaks[i], '_to_', bin_breaks[i+1])
      )
    }
  }

  result_cols <- sub('-', 'n', result_cols)

  return(result_cols)
}

# Exported Main Function
# ______________________________________________________________________________
#
#' Degree day transformation and aggregation of climate data
#'
#' The function `staggregate_degree_days()` aggregates climate data to the daily
#' level, performs a degree days transformation on these daily values, and
#' aggregates the transformed values to the polygon level and desired temporal
#' scale
#'
#' @inheritParams staggregate_custom
#'
#' @param thresholds A vector of temperature thresholds critical to a crop
#'
#' @examples
#' degree_days_output <- staggregate_degree_days(
#'   data = terra::rast(temp_nj_jun_2024_era5) - 273.15, # Climate data to transform and
#'                                          # aggregate
#'   overlay_weights = overlay_weights_nj, # Output from overlay_weights()
#'   time_agg = "month", # Sum the transformed daily values across months
#'   start_date = "2024-06-01 00:00:00", # The start date of the supplied data,
#'                                       # only required if the layer name
#'                                       # format is not compatible with stagg
#'   time_interval = "1 hour", # The temporal interval of the supplied data,
#'                             # only required if the start_date is not NA
#'   thresholds = c(0, 10, 20) # Calculate degree days above 0, 10, and 20
#'                             # degrees Celsius
#'   )
#'
#' head(degree_days_output)
#'
#' @export
staggregate_degree_days_new <- function(
    data,
    overlay_weights,
    time_agg = "month",
    start_date = NA,
    time_interval = '1 hour',
    weights_join_tolerance = 0,
    na_rm = FALSE,
    thresholds){


  # Make sure bin_breaks are ordered vector
  thresholds <- validate_bin_breaks(thresholds)

  # Create list of functions corresponding to each bin
  transformations <- lapply(
    0:length(thresholds),
    function(threshold_index){
      force(thresholds)
      force(thresholds_index)
      \(x) get_deg_day_funs(
        thresholds = thresholds,
        i = threshold_index,
        x = x
      )
    }
  )

  # Assign bin names
  result_cols <- get_deg_day_names(thresholds)

  # Automatically supply none to daily_agg
  daily_agg <- 'none'

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


