# Function to convert raster to data.table from https://gist.github.com/etiennebr/9515738
as_data_table_terra <- function(x, row.names = NULL, optional = FALSE, xy=FALSE, inmem = terra::inMemory(x), ...) {
  if(inmem) {
    v <- data.table::as.data.table(terra::as.data.frame(x, row.names=row.names, optional=optional, xy=xy, ...))
    coln <- names(x)
    if(xy) coln <- c("x", "y", coln)
    data.table::setnames(v, coln)
  } else {
    tr <- terra::blocks(x)
    l <- lapply(1:tr$n, function(i) {
      DT <- data.table::as.data.table(as.data.frame(terra::values(x, row = tr$row[i], nrows = tr$nrows[i]), ...))
      if(xy == TRUE) {
        cells <- terra::cellFromRowCol(x, c(tr$row[i], tr$row[i] + tr$nrows[i] - 1), c(1, ncol(x)))
        coords <- terra::xyFromCell(x, cell = cells[1]:cells[2])
        DT[, c("x", "y") := data.frame(terra::xyFromCell(x, cell = cells[1]:cells[2]))]
      }
      DT
    })
    v <- data.table::rbindlist(l)
    coln <- names(x)
    if(xy) {
      coln <- c("x", "y", coln)
      data.table::setcolorder(v, coln)
    }
  }
  v
}

# Function to convert raster to data.table and aggregate to daily values before transformation
daily_aggregation <- function(data, overlay_weights, daily_agg, time_interval='1 hour'){

  # Input validation and read in data
  # ----------------------------------------------------------------------------

  if(!daily_agg %in% c('average', 'sum', 'none')){

    stop(crayon::red("daily_agg must be 'average', 'sum', or 'none'"))
  }

  # Data.table of weights
  weights_dt <- data.table::as.data.table(overlay_weights)

  # read in climate data and coerce if not a spat raster
  if(inherits(data, "SpatRaster")){

    clim_stack <- data

  } else{

    clim_stack <- terra::rast(data)

  }


  # Determine whether the climate data is in climate (0 to 360) or standard (-180 to 180) coordinates
  if(terra::ext(clim_stack)$xmax > 180 + terra::xres(clim_stack)){
    coord_alignment <- "climate"
  } else{
    coord_alignment <- "standard"
  }

  # Now, determine if the polygons were likely in a different coordinate system
  if(coord_alignment == "standard"){
    polygons_split <- (
      max(weights_dt[,x]) > 180 - terra::xres(clim_stack) & # has cell next to 180
      min(weights_dt[,x]) < -180 + terra::xres(clim_stack) # has cell next to -180
    )
  }

  if(coord_alignment == "climate"){
    polygons_split <- (
      max(weights_dt[,x]) > 360 - terra::xres(clim_stack) & # has cell next to 360
      min(weights_dt[,x]) < terra::xres(clim_stack) # has cell next to 0
    )
  }

  # Check for the final case in which the polygons just span the entire globe
  if(length(unique(weights_dt[,x])) >= terra::ncol(clim_stack)){
    polygons_split <- FALSE
  }


  # If we've split the polygons, try to find a more memory efficient spot to do our cropping than across 360 / -180 line,
  # unless we have peculiar cell widths that would make this a bad idea
  if(polygons_split & 360 %% terra::xres(clim_stack) == 0){

    # Find the widest xgap that exists to be the cropping location For instance,
    # if x values are 0, 60, 120, 300, 360, then crop such that left is [0 - 120] and
    # right is [300 - 360] (saves reading in 180 degrees of data)
    x_vector <- sort(unique(weights_dt[,x]))

    # Get cropping locations by finding largest gap and making left side's xmax
    # the x value to the left and right side's xmin the x value to the right
    crop_locs <- data.frame(x_vector) |>
      dplyr::mutate(
        diff = x_vector - dplyr::lag(x_vector),
        is_right_xmin = diff == max(diff, na.rm = TRUE),
        is_left_xmax = dplyr::lead(diff) == max(diff, na.rm = TRUE)
        )

    right_xmin <- crop_locs |>
      dplyr::filter(is_right_xmin) |>
      dplyr::slice(1) |>
      dplyr::pull(x_vector)

    left_xmax <- crop_locs |>
      dplyr::filter(is_left_xmax) |>
      dplyr::slice(1) |>
      dplyr::pull(x_vector)

    # Split into 2 and then merge back together
    left_xmin <- 0 - 2*terra::xres(clim_stack)
    left_xmax <- left_xmax + 2*terra::xres(clim_stack)

    right_xmin <- right_xmin - 2*terra::xres(clim_stack)
    right_xmax <- 360 + 2*terra::xres(clim_stack)

    ymin <- min(weights_dt$y) - 2*terra::yres(clim_stack)
    ymax <- max(weights_dt$y) + 2*terra::yres(clim_stack)

    weights_ext_left <- terra::ext(left_xmin, left_xmax, ymin, ymax)
    weights_ext_right <- terra::ext(right_xmin, right_xmax, ymin, ymax)

    clim_stack_left <- terra::crop(clim_stack, weights_ext_left, snap = 'out')
    clim_stack_right <- terra::crop(clim_stack, weights_ext_right, snap = 'out')


    clim_stack <- terra::merge(clim_stack_left, clim_stack_right)

    # Get layer names (dates) from clim_stack_left
    all_layers <- names(clim_stack_left)


  } else{ # If raster overlay_weights does not span prime meridian, crop as usual

    # Extent of area weights with 2 cell buffer to make sure all cells are included
    xmin <- min(weights_dt$x) - 2*terra::xres(clim_stack)
    xmax <- max(weights_dt$x) + 2*terra::xres(clim_stack)
    ymin <- min(weights_dt$y) - 2*terra::yres(clim_stack)
    ymax <- max(weights_dt$y) + 2*terra::yres(clim_stack)

    weights_ext <- terra::ext(xmin, xmax, ymin, ymax)

    clim_stack <- terra::crop(clim_stack, weights_ext)

    all_layers <- names(clim_stack)

  }


  ## Load climate data
  ## -----------------------------------------------

  # Pass all layers through if not aggregating to daily level
  if(daily_agg == "none"){
    message(crayon::yellow("Skipping pre-transformation aggregation to daily level"))
    all_names <- names(clim_stack)
    clim_hourly <- clim_stack

    return(list(clim_hourly, all_names))
  }

  # Turn the "time_interval" argument into a number of timesteps per day.
  # Throw errors if:
  # (1) the time interval is 1 day or longer or
  # (2) the number of timesteps in a day is not a whole number

  interval_duration <- lubridate::duration(time_interval)
  day_duration <- lubridate::duration("1 day")

  if(interval_duration >= day_duration) {
    stop(crayon::red("The time interval must be less than 1 day in order to perform a daily aggregation. Please set `daily_agg` to `none` to avoid attempting daily aggregation."))
  }

  timesteps_per_day <- as.numeric(day_duration / interval_duration)

  # Check if the number of timesteps in a day is a whole number
  if(timesteps_per_day != as.integer(timesteps_per_day)) {
    stop(crayon::red("The number of timesteps in a day is not a whole number. Please change the `time_interval` argument to a number of hours that can evenly divide 24 hours."))
  }

  # Check that you have a dataset with a number of layers that is divisible by the number of timesteps in a day
  if(!(terra::nlyr(clim_stack)%%timesteps_per_day == 0)){
    stop(crayon::red(sprintf("The data does not contain a number of layers that is a multiple of %d (the number of timesteps in a day calculated using the `time_interval` argument, currently set to %s). Please use complete data with all timesteps available for each day.", timesteps_per_day, time_interval)))
  }

  layer_names <- all_layers[seq(1, length(all_layers), timesteps_per_day)] # Keep one layer name per day

  ## Aggregate to grid-day level
  ## -----------------------------------------------

  ## Average
  if(daily_agg == 'average'){

    message(crayon::green(sprintf("Averaging over %d layers per day to get daily values", timesteps_per_day)))

    # Average over each set of layers representing one day
    indices<-rep(1:(terra::nlyr(clim_stack)/timesteps_per_day),each=timesteps_per_day)
    clim_daily <- terra::tapp(clim_stack, indices, fun=mean)
  }

  ## Sum
  if(daily_agg == 'sum'){

    message(crayon::green(sprintf("Summing over %d layers per day to get daily values", timesteps_per_day)))

    # Sum over each set of layers representing one day
    indices<-rep(1:(terra::nlyr(clim_stack)/timesteps_per_day),each=timesteps_per_day)
    clim_daily <- terra::tapp(clim_stack, indices, fun=sum)
  }



  # Return a list containing, in order, daily aggregated climate data as a spatRaster stack, and the layer_names created.
  return(list(clim_daily, layer_names))
}

# Function to infer date-times for spatRaster layers based on a time interval
infer_layer_datetimes <- function(data, start_date, time_interval) {

  # Running rast() on a stack that's already a spat raster removes the values for some reason
  # read in climate data and coerce if not a spat raster
  if(methods::is(data, "SpatRaster")){

    clim_stack <- data

  } else{

    clim_stack <- terra::rast(data)

  }

  # Number of layers in the spatRaster stack
  num_layers <- terra::nlyr(clim_stack)

  # Convert start date to POSIXct
  start_date <- lubridate::as_datetime(start_date)

  # Generate the sequence of date-times for each layer
  layer_dates <- seq(start_date, by = time_interval, length.out = num_layers)

  # Make sure the full date shows up in the string every time
  formatted_dates <- format(layer_dates, "X%Y.%m.%d.%H.%M.%S")

  # Assign the inferred date-times to the spatRaster layers
  names(clim_stack) <- as.character(formatted_dates)

  return(clim_stack)
}

# Function to merge with geoweights and aggregate by polygon
polygon_aggregation <- function(clim_dt, weights_dt, list_names, time_agg, weights_join_tolerance_x, weights_join_tolerance_y, na_rm){

  ## Merge weights with climate spatRaster
  ## -----------------------------------------------

  # Set key column in the climate data.table
  keycols = c("x", "y")
  data.table::setkeyv(clim_dt, keycols)


  # Convert layer names to dates
  clim_dt[, date := stringr::str_replace(date, "^[^0-9]+", "")] # Remove any non-digit characters from the start of the string
  clim_dt[, date := lubridate::as_datetime(date)]

  # Join overlay_weights and climate data table
  if(weights_join_tolerance_x == 0 & weights_join_tolerance_y == 0){

    # Keyed merge on the x/y column
    merged_dt <- clim_dt[weights_dt, allow.cartesian = TRUE] # cols: x, y, date, value cols 1:k, poly_id, w_area, weight (if weights = T)

  } else{

    # Create a copy of the weights data table so we don't change the user's data
    # without their knowledge. (data.table creates "shallow copies" unless
    # explicitly told to do otherwise).
    copied_weights_dt <- data.table::copy(weights_dt)

    # Since data.table doesn't allow arithmatic in nonequi-joins, create tolerance columns now
    copied_weights_dt[, x_low := x - weights_join_tolerance_x]
    copied_weights_dt[, x_high := x + weights_join_tolerance_x]
    copied_weights_dt[, y_low := y - weights_join_tolerance_y]
    copied_weights_dt[, y_high := y + weights_join_tolerance_y]

    #  Remove x and y columns to avoid confusion in join
    copied_weights_dt[, x := NULL]
    copied_weights_dt[, y := NULL]


    # Determine which columns to keep in join,
    # format as arguments that can directly be supplied to j
    # through eval and parse
    cols_to_keep <- c(

      # Add x and y in manually (referring to clim_dt$x by "x.x" prevents the
      # column from being overwritten in the nonequi-join below)
      'x.x',
      'x.y',

      # Keep all column names except the tolerance columns created above
      setdiff(colnames(clim_dt), c('x', 'y')),
      setdiff(colnames(copied_weights_dt), c('x_low', 'x_high', 'y_low', 'y_high')))


    # Merge based on tolerance columns
    merged_dt <- clim_dt[copied_weights_dt, # Right join
                         allow.cartesian = TRUE,

                         # Specify which columns to keep
                         j = ..cols_to_keep,
                         on = .(x >= x_low,
                                x <= x_high,
                                y >= y_low,
                                y <= y_high)]
    data.table::setnames(merged_dt, c('x.x', 'x.y'), c('x', 'y'))
  }

  ## Multiply weights x climate value (all 1:k values); aggregate by month and polygon
  ## -----------------------------------------------

  # Multiply by secondary weights if weights = TRUE (already normalized by polygon area)
  # Otherwise multiply by just area weights
  if("weight" %in% names(merged_dt)){
    merged_dt[, (list_names) := lapply(list_names, function(x) {get(x) * weight})]
    weight_col <- "weight"
  } else {
    merged_dt[, (list_names) := lapply(list_names, function(x) {get(x) * w_area})]
    weight_col <- "w_area"
  }

  # Separate year, month, day, and time columns
  merged_dt[, ':=' (year = lubridate::year(date),
                    month = lubridate::month(date),
                    day = lubridate::day(date),
                    hour = lubridate::hour(date),
                    minute = lubridate::minute(date))]


  # Determine if we need to sum weights to "ignore" NAs
  if(na_rm){
    sum_cols <- c(weight_col, list_names)
  }else{
    sum_cols <- list_names
  }

  # Temporal Aggregation
  if(time_agg == "year"){
    message(crayon::green("Aggregating by polygon and year"))

    sum_by_poly <- merged_dt[,  lapply(.SD, sum), by = .(poly_id, year),
                             .SDcols = sum_cols]

    ## Order columns
    data.table::setcolorder(sum_by_poly, neworder = c('year', 'poly_id', sum_cols))

  }else if(time_agg == "month"){
    message(crayon::green("Aggregating by polygon and month"))

    sum_by_poly <- merged_dt[,  lapply(.SD, sum), by = .(poly_id, year, month),
                             .SDcols = sum_cols]

    ## Order columns
    data.table::setcolorder(sum_by_poly, neworder = c('year', 'month', 'poly_id', sum_cols))

  }else if(time_agg == "day"){
    message(crayon::green("Aggregating by polygon and date"))

    sum_by_poly <- merged_dt[,  lapply(.SD, sum), by = .(poly_id, year, month, day),
                             .SDcols = sum_cols]

    ## Order columns
    data.table::setcolorder(sum_by_poly, neworder = c('year', 'month', 'day', 'poly_id', sum_cols))

  }else if(time_agg == "hour"){
    message(crayon::green("Aggregating by polygon and hour"))
    sum_by_poly <- merged_dt[,  lapply(.SD, sum), by = .(poly_id, year, month, day, hour),
                             .SDcols = sum_cols]

    ## Order columns
    data.table::setcolorder(sum_by_poly, neworder = c('year', 'month', 'day', 'hour', 'poly_id', sum_cols))

  }else{

    message(crayon::green("Aggregating by polygon and time"))
    sum_by_poly <- merged_dt[,  lapply(.SD, sum), by = .(poly_id, year, month, day, hour, minute),
                             .SDcols = sum_cols]

    ## Order columns
    data.table::setcolorder(sum_by_poly, neworder = c('year', 'month', 'day', 'hour', 'minute', 'poly_id', sum_cols))
  }


  # If we need to re-weight to remove influence of NAs on means, divide by the
  # now summed weights. For instance, if .25 of the cells (by weight) was NA,
  if(na_rm){
    if(weight_col == "weight"){
      sum_by_poly[, (list_names) := lapply(list_names, function(x) {get(x) / weight})]
    } else{
      sum_by_poly[, (list_names) := lapply(list_names, function(x) {get(x) / w_area})]
    }

    # Remove the NA data
    stats::na.omit(sum_by_poly, "year")
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
#' polynomial_output <- staggregate_polynomial(
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
#'   degree = 4 # Highest order
#'   )
#'
#' head(polynomial_output)
#'
#' @export
staggregate_polynomial <- function(data, overlay_weights, daily_agg, time_agg = "month", start_date = NA, time_interval = "1 hour", weights_join_tolerance = 0, na_rm = FALSE, degree){

  # If the start date is supplied, overwrite the spatRaster's layer names to reflect the specified temporal metadata
  if(!is.na(start_date)){
    message(crayon::green(sprintf("Rewriting the data's temporal metadata (layer names) to reflect a dataset starting on the supplied start date and with a temporal interval of %s", time_interval)))
    data <- infer_layer_datetimes(data, start_date, time_interval)
  }

  # Change daily_agg to "none" if time_agg is "hour"
  if(time_agg == "hour" & daily_agg != "none"){
    message(crayon::yellow("Hourly output requested. Automatically setting daily_agg to \'none\'"))
    daily_agg = "none"
  }

  # If tolerance supplied is one number, apply to both x and y.
  # If two numbers, apply first to x and second to y
  if(length(weights_join_tolerance) == 1){
    weights_join_tolerance_x <- weights_join_tolerance
    weights_join_tolerance_y <- weights_join_tolerance
  } else if(length(weights_join_tolerance == 2)){
    weights_join_tolerance_x <- weights_join_tolerance[1]
    weights_join_tolerance_y <- weights_join_tolerance[2]
  } else{
    stop(crayon::red("Please provide one digit or a vector of only two digits for weights_join_tolerance."))
  }

  # Aggregate climate data to daily values
  setup_list <- daily_aggregation(data, overlay_weights, daily_agg, time_interval)

  clim_daily <- setup_list[[1]] # Pulls the daily aggregated spatRaster stack
  layer_names <- setup_list[[2]] # Pulls the saved layer names

  # Polynomial transformation
  poly_orders <- 1:degree # Compute values from 1 to degree
  list_length <- length(poly_orders) # How many lists are in the final object
  list_names <- sapply(1:list_length, FUN=function(x){paste("order", poly_orders[x], sep="_")})

  message(crayon::green("Executing polynomial transformation"))

  # For each daily layer, raise the value to degree, degree-1, degree-2 etc. until 1
  r <- lapply(poly_orders, FUN=function(x){clim_daily ^ x})


  ## Function: Set names of data.table by month, change from wide to long format, rename based on polynomial orders
  create_dt <- function(x){

    # Should output spatRaster cells x/y with 365 days as column names
    dt <- as_data_table_terra(r[[x]], xy=TRUE)

    # Set column names with months
    new_names <- c('x', 'y', layer_names)
    data.table::setnames(dt, new_names)

    # Change from wide to long format
    dt = data.table::melt(dt, id.vars = c("x", "y"))

    # Update variable names
    var_names <- c('date', list_names[x])
    data.table::setnames(dt, old=c('variable', 'value'), new=var_names)
  }

  # Make each layer a data.table
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
  sum_by_poly <- polygon_aggregation(
    clim_dt,
    overlay_weights,
    list_names,
    time_agg,
    weights_join_tolerance_x = weights_join_tolerance_x,
    weights_join_tolerance_y = weights_join_tolerance_y,
    na_rm = na_rm)

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
#' @param data The spatRaster stack with the data to be transformed and aggregated
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
staggregate_spline <- function(data, overlay_weights, daily_agg, time_agg = "month", start_date = NA, time_interval = "1 hour", weights_join_tolerance = 0, na_rm = FALSE, knot_locs){

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

  # If tolerance supplied is one number, apply to both x and y.
  # If two numbers, apply first to x and second to y
  if(length(weights_join_tolerance) == 1){
    weights_join_tolerance_x <- weights_join_tolerance
    weights_join_tolerance_y <- weights_join_tolerance
  } else if(length(weights_join_tolerance == 2)){
    weights_join_tolerance_x <- weights_join_tolerance[1]
    weights_join_tolerance_y <- weights_join_tolerance[2]
  } else{
    stop(crayon::red("Please provide one digit or a vector of only two digits for weights_join_tolerance."))
  }

  # Aggregated climate data to daily values
  setup_list <- daily_aggregation(data, overlay_weights, daily_agg, time_interval)

  clim_daily <- setup_list[[1]] # Pulls the daily aggregated spatRaster stack
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
      clim_daily_table <- terra::values(clim_daily)

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
      terra::values(clim_daily_new) <- clim_daily_table

      return(clim_daily_new)



    }
  }


  message(crayon::green("Executing spline transformation"))

  # For each layer, create new spline variables
  r <- lapply(0:list_length, get_spline)


  ## Function: Set names of data.table by month, change from wide to long format, rename based on polynomial orders
  create_dt <- function(x){

    # Should output raster cells x/y with 365 days as column names
    dt <- as_data_table_terra(r[[x]], xy=TRUE)

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
  sum_by_poly <- polygon_aggregation(
    clim_dt,
    overlay_weights,
    list_names,
    time_agg,
    weights_join_tolerance_x = weights_join_tolerance_x,
    weights_join_tolerance_y = weights_join_tolerance_y,
    na_rm = na_rm)

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
staggregate_bin <- function(data, overlay_weights, daily_agg, time_agg = "month", start_date = NA, time_interval = '1 hour', weights_join_tolerance = 0, na_rm = FALSE, bin_breaks){

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


  # If tolerance supplied is one number, apply to both x and y.
  # If two numbers, apply first to x and second to y
  if(length(weights_join_tolerance) == 1){
    weights_join_tolerance_x <- weights_join_tolerance
    weights_join_tolerance_y <- weights_join_tolerance
  } else if(length(weights_join_tolerance == 2)){
    weights_join_tolerance_x <- weights_join_tolerance[1]
    weights_join_tolerance_y <- weights_join_tolerance[2]
  } else{
    stop(crayon::red("Please provide one digit or a vector of only two digits for weights_join_tolerance."))
  }

  # Aggregate climate data to daily values
  setup_list <- daily_aggregation(data, overlay_weights, daily_agg, time_interval)
  clim_daily <- setup_list[[1]] # Pulls the daily aggregated raster brick
  layer_names <- setup_list[[2]] # Pulls the saved layer names


  clim_daily_table <- terra::values(clim_daily)

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
    clim_daily_table <- terra::values(clim_daily)

    if(x == 0){
      clim_daily_table <- ifelse(min(bin_breaks) > clim_daily_table, 1, 0)
    }
    else if(x == length(bin_breaks)){
      clim_daily_table <- ifelse(max(bin_breaks) <= clim_daily_table, 1, 0)
    }
    else{
      clim_daily_table <- ifelse(bin_breaks[x] <= clim_daily_table &
                                   bin_breaks[x + 1] > clim_daily_table, 1, 0)
    }


    clim_daily_new <- clim_daily
    terra::values(clim_daily_new) <- clim_daily_table

    return(clim_daily_new)
  }

  message(crayon::green("Executing binning transformation"))

  # For each bin, create new brick of binary values, including edge bins which go from -inf to min, max to inf
  r <- lapply(0:(length(bin_breaks)), FUN = check_bins)


  create_dt <- function(x){

    # Should output raster cells x/y with 365 days as column names
    dt <- as_data_table_terra(r[[x]], xy=TRUE)

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
  sum_by_poly <- polygon_aggregation(
    clim_dt,
    overlay_weights,
    list_names,
    time_agg,
    weights_join_tolerance_x = weights_join_tolerance_x,
    weights_join_tolerance_y = weights_join_tolerance_y,
    na_rm = na_rm)

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
#'  required if the `start_date` argument is not `NA`.
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
staggregate_degree_days <- function(data, overlay_weights, time_agg = "month", start_date = NA, time_interval = '1 hour', weights_join_tolerance = 0, na_rm = FALSE, thresholds){

  # If the start date is supplied, overwrite the raster's layer names to reflect the specified temporal metadata
  if(!is.na(start_date)){
    message(crayon::green(sprintf("Rewriting the data's temporal metadata (layer names) to reflect a dataset starting on the supplied start date and with a temporal interval of %s", time_interval)))
    data <- infer_layer_datetimes(data, start_date, time_interval)
  }

  # If tolerance supplied is one number, apply to both x and y.
  # If two numbers, apply first to x and second to y
  if(length(weights_join_tolerance) == 1){
    weights_join_tolerance_x <- weights_join_tolerance
    weights_join_tolerance_y <- weights_join_tolerance
  } else if(length(weights_join_tolerance == 2)){
    weights_join_tolerance_x <- weights_join_tolerance[1]
    weights_join_tolerance_y <- weights_join_tolerance[2]
  } else{
    stop(crayon::red("Please provide one digit or a vector of only two digits for weights_join_tolerance."))
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
    clim_table <- terra::values(clim_rast)
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
    terra::values(clim_rast_new) <- clim_table

    return(clim_rast_new)
  }

  message(crayon::green("Executing degree days transformation"))

  # For each bin, create new brick of binary values, including edge bins which go from -inf to min, max to inf
  r <- lapply(0:(length(thresholds)), FUN = calc_deg_days)


  create_dt <- function(x){

    # Should output raster cells x/y with 365 days as column names
    dt <- as_data_table_terra(r[[x]], xy=TRUE)

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
  sum_by_poly <- polygon_aggregation(
    clim_dt,
    overlay_weights,
    list_names,
    time_agg,
    weights_join_tolerance_x = weights_join_tolerance_x,
    weights_join_tolerance_y = weights_join_tolerance_y,
    na_rm = na_rm)

  return(sum_by_poly)

}
