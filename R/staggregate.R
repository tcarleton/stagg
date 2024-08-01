# Function to convert raster to data.table from https://gist.github.com/etiennebr/9515738
as.data.table.raster <- function(x, row.names = NULL, optional = FALSE, xy=FALSE, inmem = raster::canProcessInMemory(x, 2), ...) {
  if(inmem) {
    v <- data.table::as.data.table(raster::as.data.frame(x, row.names=row.names, optional=optional, xy=xy, ...))
  } else {
    tr <- raster::blockSize(x, n=2)
    l <- lapply(1:tr$n, function(i)
      data.table::as.data.table(raster::as.data.frame(raster::getValues(x,
                                                                        row=tr$row[i],
                                                                        nrows=tr$nrows[i]),
                                                      row.names=row.names, optional=optional, xy=xy, ...)))
    v <- data.table::rbindlist(l)
  }
  coln <- names(x)
  if(xy) coln <- c("x", "y", coln)
  data.table::setnames(v, coln)
  v
}

# Function to convert raster to data.table and aggregate to daily values before transformation
daily_aggregation <- function(data, overlay_weights, daily_agg, time_interval='1 hour'){

  if(!daily_agg %in% c('average', 'sum', 'none')){

    stop(crayon::red("daily_agg must be 'average', 'sum', or 'none'"))
  }

  # Data.table of weights
  weights_dt <- data.table::as.data.table(overlay_weights)

  # Check if raster overlay_weights span prime meridian
  is_pm <- FALSE
  for(i in length(weights_dt$x)){
    if(weights_dt[i,x] < .5 | weights_dt[i,x] > 359.5){
      is_pm = TRUE # Check each x value
      break # If near 0 value found exit loop
    }
  }

  ## read in climate data
  clim_raster <- raster::stack(data)

  ## shift into 0 to 360 if not already in that format
  if(raster::extent(clim_raster)@xmin < 0 - raster::xres(clim_raster) / 2) {

    clim_raster_xmin <- raster::extent(clim_raster)@xmin - raster::xres(clim_raster) / 2
    clim_raster_xmax <- raster::extent(clim_raster)@xmax + raster::xres(clim_raster) / 2
    clim_raster_ymin <- raster::extent(clim_raster)@ymin - raster::yres(clim_raster) / 2
    clim_raster_ymax <- raster::extent(clim_raster)@ymax + raster::yres(clim_raster) / 2

    clim_raster1 <- raster::crop(clim_raster, c(clim_raster_xmin, min(clim_raster_xmax, 0), clim_raster_ymin, clim_raster_ymax))

    clim_raster1 <- raster::shift(clim_raster1, dx = 360)

    if(raster::extent(clim_raster)@xmax > 0) {

      clim_raster2 <- raster::crop(clim_raster, c(0, clim_raster_xmax, clim_raster_ymin, clim_raster_ymax))

      clim_raster <- raster::merge(clim_raster1, clim_raster2)

    } else {

      clim_raster <- clim_raster1
    }

  }

  # If raster overlay_weights does not span prime meridian, crop as usual
  if(is_pm == FALSE){
    # Extent of area weights with 2 cell buffer to make sure all cells are included
    min_x <- min(weights_dt$x) - 2*raster::xres(clim_raster)
    max_x <- max(weights_dt$x) + 2*raster::xres(clim_raster)
    min_y <- min(weights_dt$y) - 2*raster::yres(clim_raster)
    max_y <- max(weights_dt$y) + 2*raster::yres(clim_raster)

    weights_ext <- raster::extent(min_x, max_x, min_y, max_y)

    clim_raster <- raster::crop(clim_raster, weights_ext)

    all_layers <- names(clim_raster)

  }

  else { # If yes, make 2 crops and stitch together
    min_x_left <- min(weights_dt$x[weights_dt$x >= 180]) - 2*raster::xres(data)
    max_x_left <- 360

    min_x_right <- -0.1 # This is necessary because cropping is non-inclusive (min_x_right = 0 excludes the cell at 0)
    max_x_right <- max(weights_dt$x[weights_dt$x < 180]) + 2*raster::xres(data)

    min_y <- min(weights_dt$y) - 2*raster::yres(data)
    max_y <- max(weights_dt$y) + 2*raster::yres(data)

    weights_ext_left <- raster::extent(min_x_left, max_x_left, min_y, max_y)
    weights_ext_right <- raster::extent(min_x_right, max_x_right, min_y, max_y)

    clim_raster_left <- raster::crop(clim_raster, weights_ext_left)
    clim_raster_right <- raster::crop(clim_raster, weights_ext_right)


    clim_raster <- raster::stack(raster::merge(clim_raster_left, clim_raster_right)) # Merge creates raster bricks without proper layer names

    # Get layer names (dates) from clim_raster_left
    all_layers <- names(clim_raster_left)
  }


  ## Load climate data
  ## -----------------------------------------------

  # Pass all layers through if not aggregating to daily level
  if(daily_agg == "none"){
    message(crayon::yellow("Skipping pre-transformation aggregation to daily level"))
    all_names <- names(clim_raster)
    clim_hourly <- clim_raster

    return(list(clim_hourly, all_names))
  }

  # Turn the "time_interval" argument into a number of timesteps per day.
  # Throw errors if:
  # (1) the time interval is 1 day or longer or
  # (2) the number of timesteps in a day is not a whole number

  interval_duration <- duration(time_interval)
  day_duration <- duration("1 day")

  if(interval_duration >= day_duration) {
    stop(crayon::red("The time interval must be less than 1 day in order to perform a daily aggregation. Please set `daily_agg` to `none` to avoid attempting daily aggregation."))
  }

  timesteps_per_day <- as.numeric(day_duration / interval_duration)

  # Check if the number of timesteps in a day is a whole number
  if(timesteps_per_day != as.integer(timesteps_per_day)) {
    stop(crayon::red("The number of timesteps in a day is not a whole number. Please change the `time_interval` argument to reflect a dataset wherein all layers correspond to a specific day in order to perform the daily aggregation."))
  }

  # Check that you have a dataset with a number of layers that is divisible by the number of timesteps in a day
  if(!(raster::nlayers(clim_raster)%%timesteps_per_day == 0)){
    stop(crayon::red(sprintf("The data does not contain a number of layers that is a multiple of %d (the number of timesteps in a day calculated using the `time_interval` argument, currently set to %s). Please use complete data with all timesteps available for each day.", timesteps_per_day, time_interval)))
  }

  layer_names <- all_layers[seq(1, length(all_layers), timesteps_per_day)] # Keep one layer name per day

  ## Aggregate to grid-day level
  ## -----------------------------------------------

  ## Average
  if(daily_agg == 'average'){

    message(crayon::green(sprintf("Averaging over %d layers per day to get daily values", timesteps_per_day)))

    # Average over each set of layers representing one day
    indices<-rep(1:(raster::nlayers(clim_raster)/timesteps_per_day),each=timesteps_per_day)
    clim_daily <- raster::stackApply(clim_raster, indices = indices, fun=mean)
  }

  ## Sum
  if(daily_agg == 'sum'){

    message(crayon::green(sprintf("Summing over %d layers per day to get daily values", timesteps_per_day)))

    # Sum over each set of layers representing one day
    indices<-rep(1:(raster::nlayers(clim_raster)/timesteps_per_day),each=timesteps_per_day)
    clim_daily <- raster::stackApply(clim_raster, indices = indices, fun=sum)
  }



  # Return a list containing, in order, daily aggregated climate data as a raster brick, and the layer_names created.
  return(list(clim_daily, layer_names))
}

# Function to infer date-times for raster layers based on a time interval
infer_layer_datetimes <- function(raster_stack, start_date, time_interval) {
  # Turn the data into a raster stack in case it's a SpatRaster
  raster_stack <- raster::stack(data)

  # Number of layers in the raster stack
  num_layers <- raster::nlayers(raster_stack)

  # Convert start date to POSIXct
  start_date <- lubridate::as_datetime(start_date)

  # Generate the sequence of date-times for each layer
  layer_dates <- seq(start_date, by = time_interval, length.out = num_layers)

  # Make sure the full date shows up in the string every time
  formatted_dates <- format(layer_dates, "X%Y.%m.%d.%H.%M.%S")

  # Assign the inferred date-times to the raster layers
  names(raster_stack) <- as.character(formatted_dates)

  return(raster_stack)
}

# Function to merge with geoweights and aggregate by polygon
polygon_aggregation <- function(clim_dt, weights_dt, list_names, time_agg){

  ## Merge weights with climate raster
  ## -----------------------------------------------

  # Set key column in the climate data.table
  keycols = c("x", "y")
  data.table::setkeyv(clim_dt, keycols)


  # Convert layer names to dates
  clim_dt[, date := stringr::str_replace(date, "^[^0-9]+", "")] # Remove any non-digit characters from the start of the string
  clim_dt[, date := lubridate::as_datetime(date)]

  # Keyed merge on the x/y column
  merged_dt <- clim_dt[weights_dt, allow.cartesian = TRUE] # cols: x, y, date, value cols 1:k, poly_id, w_area, weight (if weights = T)

  ## Multiply weights x climate value (all 1:k values); aggregate by month and polygon
  ## -----------------------------------------------

  # Multiply by secondary weights if weights = TRUE (already normalized by polygon area)
  # Otherwise multiply by just area weights
  if("weight" %in% names(merged_dt)){
    merged_dt[, (list_names) := lapply(list_names, function(x) {get(x) * weight})]
  } else {
    merged_dt[, (list_names) := lapply(list_names, function(x) {get(x) * w_area})]
  }

  # Separate year, month, day, and time columns
  merged_dt[, ':=' (year = lubridate::year(date),
                    month = lubridate::month(date),
                    day = lubridate::day(date),
                    hour = lubridate::hour(date))] # NOTE: if the timestep is smaller than hourly, the minutes and seconds will not be recorded. However, this does not affect the output since the lowest aggregation level possible is 'hour'.

  # Temporal Aggregation
  if(time_agg == "year"){
    message(crayon::green("Aggregating by polygon and year"))

    sum_by_poly <- merged_dt[,  lapply(.SD, sum), by = .(poly_id, year),
                             .SDcols = list_names]

    ## Order columns
    data.table::setcolorder(sum_by_poly, neworder = c('year', 'poly_id', list_names))
  }

  if(time_agg == "month"){
    message(crayon::green("Aggregating by polygon and month"))

    sum_by_poly <- merged_dt[,  lapply(.SD, sum), by = .(poly_id, year, month),
                             .SDcols = list_names]

    ## Order columns
    data.table::setcolorder(sum_by_poly, neworder = c('year', 'month', 'poly_id', list_names))
  }

  if(time_agg == "day"){
    message(crayon::green("Aggregating by polygon"))

    sum_by_poly <- merged_dt[,  lapply(.SD, sum), by = .(poly_id, year, month, day),
                             .SDcols = list_names]

    ## Order columns
    data.table::setcolorder(sum_by_poly, neworder = c('year', 'month', 'day', 'poly_id', list_names))
  }




  if(time_agg == "hour"){
    message(crayon::green("Aggregating by polygon"))
    sum_by_poly <- merged_dt[,  lapply(.SD, sum), by = .(poly_id, year, month, day, hour),
                             .SDcols = list_names]

    ## Order columns
    data.table::setcolorder(sum_by_poly, neworder = c('year', 'month', 'day', 'hour', 'poly_id', list_names))
  }


  ## Return the sums by polygon
  return(sum_by_poly)



}


#=====================================================================================================================================================

#' Polynomial transformation and aggregation of climate data
#'
#' The function `staggregate_polynomial()` aggregates climate data to the daily
#' level, raises these daily values to the 1 through nth power, and aggregates
#' the transformed values to the polygon level and desired temporal scale.
#'
#' @param data The raster brick with the data to be transformed and aggregated
#' @param overlay_weights A table of weights which can be generated using the
#'   function `overlay_weights()`
#' @param daily_agg How to aggregate hourly values to daily values prior to
#'   transformation. Options are `'sum'`, `'average'`, or `'none'` (`'none'`
#'   will transform values without first aggregating to the daily level)
#' @param time_agg the temporal scale to aggregate data to. Options are
#'   `'hour'`, `'day'`, `'month'`, or `'year'` (`'hour'` cannot be selected
#'   unless `daily_agg = 'none'`)
#' @param start_date the date (and time, if applicable) of the first layer in
#'  the raster. To be input in a format compatible with
#'  lubridate::as_datetime(), e.g. `"1991-10-29"` or `"1991-10-29 00:00:00"`.
#'  The default is `NA` since the rasters usually already contain temporal
#'  information in the layer names and they do not need to be manually supplied.
#' @param time_interval the time interval between layers in the raster to be
#'  aggregated. To be input in a format compatible with seq(), e.g.
#'  `'1 day'` or `'3 months'`. The default is `'1 hour'` and this argument is
#'  required if daily_agg is not `'none'` or if the layer name format is not
#'  compatible with stagg.
#' @param degree the highest exponent to raise the data to
#'
#' @examples
#' polynomial_output <- staggregate_polynomial(
#'   data = temp_kansas_jan_2020_era5, # Climate data to transform and aggregate
#'   overlay_weights = overlay_weights_kansas, # Output from overlay_weights()
#'   daily_agg = "average", # Average hourly values to produce daily values
#'                          # before transformation
#'   time_agg = "month", # Sum the transformed daily values across months
#'   start_date = "2020-01-01 00:00:00", # The start date of the supplied data, only required if the layer name format is not compatible with stagg
#'   time_interval = "1 hour", # The temporal interval of the supplied data, required if daily_agg is not "none" or if the layer name format is not compatible with stagg
#'   degree = 4 # Highest order
#'   )
#'
#' head(polynomial_output)
#'
#' @export
staggregate_polynomial <- function(data, overlay_weights, daily_agg, time_agg = "month", start_date = NA, time_interval = "1 hour", degree){

  # If the start date is supplied, overwrite the raster's layer names to reflect the specified temporal metadata
  if(!is.na(start_date)){
    message(crayon::green(sprintf("Rewriting the data's temporal metadata (layer names) to reflect a dataset starting on the supplied start date and with a temporal interval of %s", time_interval)))
    data <- infer_layer_datetimes(data, start_date, time_interval)
  }

  # Change daily_agg to "none" if time_agg is "hour"
  if(time_agg == "hour" & daily_agg != "none"){
    message(crayon::yellow("Hourly output requested. Automatically setting daily_agg to \'none\'"))
    daily_agg = "none"
  }

  # Aggregate climate data to daily values
  setup_list <- daily_aggregation(data, overlay_weights, daily_agg, time_interval)

  clim_daily <- setup_list[[1]] # Pulls the daily aggregated raster brick
  layer_names <- setup_list[[2]] # Pulls the saved layer names

  # Polynomial transformation
  poly_orders <- seq(1:degree) # Compute values from 1 to degree
  list_length <- length(poly_orders) # How many lists are in the final object
  list_names <- sapply(1:list_length, FUN=function(x){paste("order", poly_orders[x], sep="_")})

  message(crayon::green("Executing polynomial transformation"))

  # For each daily layer, raise the value to degree, degree-1, degree-2 etc. until 1
  r <- lapply(poly_orders, FUN=function(x){clim_daily ^ x})


  ## Function: Set names of data.table by month, change from wide to long format, rename based on polynomial orders
  create_dt <- function(x){

    # Should output raster cells x/y with 365 days as column names
    dt <- as.data.table.raster(r[[x]], xy=TRUE)

    # Set column names with months
    new_names <- c('x', 'y', layer_names)
    data.table::setnames(dt, new_names)

    # Change from wide to long format
    dt = data.table::melt(dt, id.vars = c("x", "y"))

    # Update variable names
    var_names <- c('date', list_names[x])
    data.table::setnames(dt, old=c('variable', 'value'), new=var_names)
  }

  # Make each raster layer a data.table
  list_dt <- lapply(1:list_length, create_dt)

  # Merge all data.tables together if there are multiple
  clim_dt <- list_dt[[1]]
  if(list_length > 1 ){
    for(i in 2:list_length){
      dt_m <- list_dt[[i]]
      clim_dt <- merge(clim_dt, dt_m, by=c('x', 'y', 'date'))
    }
  }

  # Aggregate by polygon
  sum_by_poly <- polygon_aggregation(clim_dt, overlay_weights, list_names, time_agg)

  return(sum_by_poly)

}









#==================================================================================================================================================

#' Restricted cubic spline transformation and aggregation of climate data
#'
#' The function `staggregate_spline()` aggregates climate data to the daily
#' level, performs a restricted cubic spline transformation on these daily
#' values, and aggregates the transformed values to the polygon level and
#' desired temporal scale.
#'
#' @param data The raster brick with the data to be transformed and aggregated
#' @param overlay_weights A table of weights which can be generated using the
#'   function `overlay_weights()`
#' @param daily_agg How to aggregate hourly values to daily values prior to
#'   transformation. Options are `'sum'`, `'average'`, or `'none'` (`'none'`
#'   will transform values without first aggregating to the daily level)
#' @param time_agg the temporal scale to aggregate data to. Options are
#'   `'hour'`, `'day'`, `'month'`, or `'year'` (`'hour'` cannot be selected
#'   unless `daily_agg = 'none'`)
#' @param start_date the date (and time, if applicable) of the first layer in
#'  the raster. To be input in a format compatible with
#'  lubridate::as_datetime(), e.g. `"1991-10-29"` or `"1991-10-29 00:00:00"`.
#'  The default is `NA` since the rasters usually already contain temporal
#'  information in the layer names and they do not need to be manually supplied.
#' @param time_interval the time interval between layers in the raster to be
#'  aggregated. To be input in a format compatible with seq(), e.g.
#'  `'1 day'` or `'3 months'`. The default is `'1 hour'` and this argument is
#'  required if daily_agg is not `'none'` or if the layer name format is not
#'  compatible with stagg.
#' @param knot_locs where to place the knots
#'
#' @examples
#' spline_output <- staggregate_spline(
#' data = temp_kansas_jan_2020_era5, # Climate data to transform and aggregate
#' overlay_weights = overlay_weights_kansas, # Output from overlay_weights()
#' daily_agg = "average", # Average hourly values to produce daily values before
#'                        # transformation
#' time_agg = "month", # Sum the transformed daily values across months
#' start_date = "2020-01-01 00:00:00", # The start date of the supplied data, only required if the layer name format is not compatible with stagg
#' time_interval = "1 hour", # The temporal interval of the supplied data, required if daily_agg is not "none" or if the layer name format is not compatible with stagg
#' knot_locs = c(0, 7.5, 12.5, 20) # Where to place knots
#' )
#'
#' head(spline_output)
#'
#' @export
staggregate_spline <- function(data, overlay_weights, daily_agg, time_agg = "month", start_date = NA, time_interval = "1 hour", knot_locs){

  # If the start date is supplied, overwrite the raster's layer names to reflect the specified temporal metadata
  if(!is.na(start_date)){
    message(crayon::green(sprintf("Rewriting the data's temporal metadata (layer names) to reflect a dataset starting on the supplied start date and with a temporal interval of %s", time_interval)))
    data <- infer_layer_datetimes(data, start_date, time_interval)
  }

  # Change daily_agg to "none" if time_agg is "hour"
  if(time_agg == "hour" & daily_agg != "none"){
    message(crayon::yellow("Hourly output requested. Automatically setting daily_agg to \'none\'"))
    daily_agg = "none"
    }

  # Aggregated climate data to daily values
  setup_list <- daily_aggregation(data, overlay_weights, daily_agg, time_interval)

  clim_daily <- setup_list[[1]] # Pulls the daily aggregated raster brick
  layer_names <- setup_list[[2]] # Pulls the saved layer names


  # Spline transformation
  knot_locs <- sort(knot_locs)
  num_knots <- length(knot_locs)
  list_length <- num_knots - 2
  list_names <- sapply(0:list_length, FUN=function(x){if(x == 0){"value"}else{paste("term", x, sep="_")}})


  # Define restricted cubic spline function
  get_spline <- function(x){

    # Make first raster returned just the climate variable to preserve it's column in the resulting data.table
    if(x == 0){
      return(clim_daily)
    }
    # Add in spline terms, all of which are 0 if negative
    else{
      clim_daily_table <- raster::values(clim_daily)

       part1 <- ifelse((clim_daily_table - knot_locs[x]) > 0,
                       (clim_daily_table - knot_locs[x])^3, 0)

       part2 <- (ifelse((clim_daily_table - knot_locs[num_knots - 1]) > 0,
                 (clim_daily_table - knot_locs[num_knots - 1])^3 *
                   ((knot_locs[num_knots] - knot_locs[x]) / (knot_locs[num_knots] - knot_locs[num_knots - 1])), 0))


       part3 <- (ifelse((clim_daily_table - knot_locs[num_knots]) > 0,
                 (clim_daily_table - knot_locs[num_knots])^3 *
                   ((knot_locs[num_knots - 1] - knot_locs[x]) / (knot_locs[num_knots] - knot_locs[num_knots - 1])), 0))

      clim_daily_table <- part1 - part2 + part3

      clim_daily_new <- clim_daily
      raster::values(clim_daily_new) <- clim_daily_table

      return(clim_daily_new)



    }
  }


  message(crayon::green("Executing spline transformation"))

  # For each layer, create new spline variables
  r <- lapply(0:list_length, get_spline)


  ## Function: Set names of data.table by month, change from wide to long format, rename based on polynomial orders
  create_dt <- function(x){

    # Should output raster cells x/y with 365 days as column names
    dt <- as.data.table.raster(r[[x]], xy=TRUE)

    # Set column names with months
    new_names <- c('x', 'y', layer_names)
    data.table::setnames(dt, new_names)

    # Change from wide to long format
    dt = data.table::melt(dt, id.vars = c("x", "y"))

    # Update variable names
    var_names <- c('date', list_names[x])
    data.table::setnames(dt, old=c('variable', 'value'), new=var_names)
  }

  # Make each raster layer a data.table
  list_dt <- lapply(1:(list_length + 1), create_dt)

  # Merge all data.tables together
  clim_dt <- list_dt[[1]]
  for(i in 2:(list_length + 1)){
    dt_m <- list_dt[[i]]
    clim_dt <- merge(clim_dt, dt_m, by=c('x', 'y', 'date'))
  }


  # Aggregate by polygon
  sum_by_poly <- polygon_aggregation(clim_dt, overlay_weights, list_names, time_agg)

  return(sum_by_poly)
}

#==================================================================================================================================================

#' Bin transformation and aggregation of climate data
#'
#' The function `staggregate_bin()` aggregates climate data to the daily level,
#' splits these daily values into bins, and aggregates the transformed
#' values to the polygon level and desired temporal scale.
#'
#' @param data The raster brick with the data to be transformed and aggregated
#' @param overlay_weights A table of weights which can be generated using the
#'   function `overlay_weights()`
#' @param daily_agg How to aggregate hourly values to daily values prior to
#'   transformation. Options are `'sum'`, `'average'`, or `'none'` (`'none'`
#'   will transform values without first aggregating to the daily level)
#' @param time_agg the temporal scale to aggregate data to. Options are
#'   `'hour'`, `'day'`, `'month'`, or `'year'` (`'hour'` cannot be selected
#'   unless `daily_agg = 'none'`)
#' @param start_date the date (and time, if applicable) of the first layer in
#'  the raster. To be input in a format compatible with
#'  lubridate::as_datetime(), e.g. `"1991-10-29"` or `"1991-10-29 00:00:00"`.
#'  The default is `NA` since the rasters usually already contain temporal
#'  information in the layer names and they do not need to be manually supplied.
#' @param time_interval the time interval between layers in the raster to be
#'  aggregated. To be input in a format compatible with seq(), e.g.
#'  `'1 day'` or `'3 months'`. The default is `'1 hour'` and this argument is
#'  required if daily_agg is not `'none'` or if the layer name format is not
#'  compatible with stagg.
#' @param bin_breaks A vector of bin boundaries to split the data by
#'
#' @examples
#' bin_output <- staggregate_bin(
#'   data = temp_kansas_jan_2020_era5, # Climate data to transform and aggregate
#'   overlay_weights = overlay_weights_kansas, # Output from overlay_weights()
#'   daily_agg = "average", # Average hourly values to produce daily values
#'                          # before transformation
#'   time_agg = "month", # Sum the transformed daily values across months
#'   start_date = "2020-01-01 00:00:00", # The start date of the supplied data, only required if the layer name format is not compatible with stagg
#'   time_interval = "1 hour", # The temporal interval of the supplied data, required if daily_agg is not "none" or if the layer name format is not compatible with stagg
#'   bin_breaks = c(0, 2.5, 5, 7.5, 10) # Draw 6 bins from ninf to 0, 0 to 2.5,
#'                                      # 2.5 to 5, 5 to 7.5, 7.5 to 10, 10 to inf
#'   )
#'
#' head(bin_output)
#'
#' @export
staggregate_bin <- function(data, overlay_weights, daily_agg, time_agg = "month", start_date = NA, time_interval = '1 hour', bin_breaks){

  # If the start date is supplied, overwrite the raster's layer names to reflect the specified temporal metadata
  if(!is.na(start_date)){
    message(crayon::green(sprintf("Rewriting the data's temporal metadata (layer names) to reflect a dataset starting on the supplied start date and with a temporal interval of %s", time_interval)))
    data <- infer_layer_datetimes(data, start_date, time_interval)
  }

  # Change daily_agg to "none" if time_agg is "hour"
  if(time_agg == "hour" & daily_agg != "none"){
    message(crayon::yellow("Hourly output requested. Automatically setting daily_agg to \'none\'"))
    daily_agg = "none"
  }

  # Aggregate climate data to daily values
  setup_list <- daily_aggregation(data, overlay_weights, daily_agg, time_interval)
  clim_daily <- setup_list[[1]] # Pulls the daily aggregated raster brick
  layer_names <- setup_list[[2]] # Pulls the saved layer names


  clim_daily_table <- raster::values(clim_daily)

  bin_breaks <- sort(bin_breaks)


  # Create names for new columns
  list_names <- sapply(0:(length(bin_breaks)), FUN=function(x){
    if(x == 0){
      paste("bin", "ninf", "to", sub("-", "n", min(bin_breaks)), sep = "_")
    }
    else if(x == length(bin_breaks)){
      paste("bin", sub("-", "n", max(bin_breaks)), "to", "inf", sep = "_")
    }
    else{
      paste("bin", sub("-", "n", bin_breaks[x]), "to", sub("-", "n", bin_breaks[x+1]), sep = "_")
    }
  })


  # Function check_bins to determine which bins data points fall into
  check_bins <- function(x){
    clim_daily_table <- raster::values(clim_daily)

    if(x == 0){
      clim_daily_table <- ifelse(min(bin_breaks) > clim_daily_table, 1, 0)
    }
    else if(x == length(bin_breaks)){
      clim_daily_table <- ifelse(max(bin_breaks) < clim_daily_table, 1, 0)
    }
    else{
      clim_daily_table <- ifelse(bin_breaks[x] <= clim_daily_table &
                                   bin_breaks[x + 1] > clim_daily_table, 1, 0)
    }


    clim_daily_new <- clim_daily
    raster::values(clim_daily_new) <- clim_daily_table

    return(clim_daily_new)
  }

  message(crayon::green("Executing binning transformation"))

  # For each bin, create new brick of binary values, including edge bins which go from -inf to min, max to inf
  r <- lapply(0:(length(bin_breaks)), FUN = check_bins)


  create_dt <- function(x){

    # Should output raster cells x/y with 365 days as column names
    dt <- as.data.table.raster(r[[x]], xy=TRUE)

    # Set column names with months
    new_names <- c('x', 'y', layer_names)
    data.table::setnames(dt, new_names)

    # Change from wide to long format
    dt = data.table::melt(dt, id.vars = c("x", "y"))

    # Update variable names
    var_names <- c('date', list_names[x])
    data.table::setnames(dt, old=c('variable', 'value'), new=var_names)
  }

  # Make each raster layer a data.table
  list_dt <- lapply(1:(length(bin_breaks) + 1), create_dt)

  # Merge all data.tables together
  clim_dt <- list_dt[[1]]
  for(i in 2:(length(bin_breaks) + 1)){
    dt_m <- list_dt[[i]]
    clim_dt <- merge(clim_dt, dt_m, by=c('x', 'y', 'date'))
  }



  # Aggregate by polygon
  sum_by_poly <- polygon_aggregation(clim_dt, overlay_weights, list_names, time_agg)

  return(sum_by_poly)
}




#================================================================================================================================================

#' Degree days transformation and aggregation of climate data
#'
#' The function `staggregate_degree_days()` aggregates climate data to the daily
#' level, performs a degree days transformation on these daily values, and
#' aggregates the transformed values to the polygon level and desired temporal
#' scale
#'
#' @param data The raster brick with the data to be transformed and aggregated
#' @param overlay_weights A table of weights which can be generated using the
#'   function `overlay_weights()`
#' @param time_agg the temporal scale to aggregate data to. Options are `'day'`,
#'   `'month'`, or `'year'`
#' @param start_date the date (and time, if applicable) of the first layer in
#'  the raster. To be input in a format compatible with
#'  lubridate::as_datetime(), e.g. `"1991-10-29"` or `"1991-10-29 00:00:00"`.
#'  The default is `NA` since the rasters usually already contain temporal
#'  information in the layer names and they do not need to be manually supplied.
#' @param time_interval the time interval between layers in the raster to be
#'  aggregated. To be input in a format compatible with seq(), e.g.
#'  `'1 day'` or `'3 months'`. The default is `'1 hour'` and this argument is
#'  required if the layer name format is not compatible with stagg.
#' @param thresholds A vector of temperature thresholds critical to a crop
#'
#' @examples
#' degree_days_output <- staggregate_degree_days(
#'   data = temp_kansas_jan_2020_era5, # Climate data to transform and aggregate
#'   overlay_weights = overlay_weights_kansas, # Output from overlay_weights()
#'   time_agg = "month", # Sum the transformed daily values across months
#'   start_date = "2020-01-01 00:00:00", # The start date of the supplied data, only required if the layer name format is not compatible with stagg
#'   time_interval = "1 hour", # The temporal interval of the supplied data, only required if the layer name format is not compatible with stagg
#'   thresholds = c(0, 10, 20) # Calculate degree days above 0, 10, and 20
#'                             # degrees Celsius
#'   )
#'
#' head(degree_days_output)
#'
#' @export
staggregate_degree_days <- function(data, overlay_weights, time_agg = "month", start_date = NA, time_interval = '1 hour', thresholds){

  # If the start date is supplied, overwrite the raster's layer names to reflect the specified temporal metadata
  if(!is.na(start_date)){
    message(crayon::green(sprintf("Rewriting the data's temporal metadata (layer names) to reflect a dataset starting on the supplied start date and with a temporal interval of %s", time_interval)))
    data <- infer_layer_datetimes(data, start_date, time_interval)
  }

  # Run climate data through daily_aggregation)() (not actually aggregating to daily values)
  setup_list <- daily_aggregation(data, overlay_weights, daily_agg = "none")
  clim_rast <- setup_list[[1]] # Pulls the raster brick
  layer_names <- setup_list[[2]] # Pulls the saved layer names


  thresholds <- sort(thresholds)

  # Create names for new columns
  list_names <- sapply(0:(length(thresholds)), FUN=function(x){
    if(x == 0){
      paste("threshold", "ninf", "to", sub("-", "n", min(thresholds)), sep = "_")
    }
    else if(x == length(thresholds)){
      paste("threshold", sub("-", "n", max(thresholds)), "to", "inf", sep = "_")
    }
    else{
      paste("threshold", sub("-", "n", thresholds[x]), "to", sub("-", "n", thresholds[x+1]), sep = "_")
    }
  })


  # Create function to calculate degree days
  calc_deg_days <- function(x){
    clim_table <- raster::values(clim_rast)
    if(x == 0){ # For the lowest threshold, create a variable equal to 0 if the
                # value is greater than the threshold, and equal to the
                # threshold minus value otherwise
      clim_table <- ifelse(clim_table > min(thresholds), 0, min(thresholds) - clim_table)

    } else if(x == length(thresholds)){ # For the highest threshold, create a
                                      # variable equal to 0 if value is less
                                      # than threshold and equal to the value
                                      # minus the threshold otherwise
      clim_table <- ifelse(clim_table < max(thresholds), 0, clim_table - max(thresholds))

    } else{ # For all other thresholds, create variable equal to 0 if value is
          # less than threshold, equal to next threshold minus current threshold
          # if the value is greater than the next threshold, and equal to value
          # minus current threshold otherwise
      clim_table <- ifelse(clim_table < thresholds[x], 0,
                           ifelse(clim_table > thresholds[x + 1], thresholds[x + 1] - thresholds[x],
                                  clim_table - thresholds[x]))

    }

    clim_rast_new <- clim_rast
    raster::values(clim_rast_new) <- clim_table

    return(clim_rast_new)
  }

  message(crayon::green("Executing degree days transformation"))

  # For each bin, create new brick of binary values, including edge bins which go from -inf to min, max to inf
  r <- lapply(0:(length(thresholds)), FUN = calc_deg_days)


  create_dt <- function(x){

    # Should output raster cells x/y with 365 days as column names
    dt <- as.data.table.raster(r[[x]], xy=TRUE)

    # Set column names with months
    new_names <- c('x', 'y', layer_names)
    data.table::setnames(dt, new_names)

    # Change from wide to long format
    dt = data.table::melt(dt, id.vars = c("x", "y"))

    # Update variable names
    var_names <- c('date', list_names[x])
    data.table::setnames(dt, old=c('variable', 'value'), new=var_names)
  }

  # Make each raster layer a data.table
  list_dt <- lapply(1:(length(thresholds) + 1), create_dt)

  # Merge all data.tables together if there are multiple
  clim_dt <- list_dt[[1]]
  if(length(thresholds) > 1 ){
    for(i in 2:(length(thresholds) + 1)){
      dt_m <- list_dt[[i]]
      clim_dt <- merge(clim_dt, dt_m, by=c('x', 'y', 'date'))
    }
  }



  # Aggregate by polygon
  sum_by_poly <- polygon_aggregation(clim_dt, overlay_weights, list_names, time_agg)

  return(sum_by_poly)

}
