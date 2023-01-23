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
daily_aggregation <- function(data, overlay_weights, daily_agg){


  # Data.table of weights
  weights_dt <- data.table::as.data.table(overlay_weights)

  # Extent of area weights with slight buffer to make sure all cells are included
  min_x <- min(weights_dt$x) - 0.5
  max_x <- max(weights_dt$x) + 0.5
  min_y <- min(weights_dt$y) - 0.5
  max_y <- max(weights_dt$y) + 0.5

  weights_ext <- raster::extent(min_x, max_x, min_y, max_y)

  ## Load climate data
  ## -----------------------------------------------


  # Immediately crop to weights extent
  clim_raster <- raster::crop(raster::stack(data), weights_ext)


  if(!(raster::nlayers(clim_raster)%%24 == 0)){
    stop(crayon::red("Data does not contain a number of layers that is a multiple of 24. Please use hourly data representing a whole number of days."))
  }

  # Get layer names (dates)
  all_layers <- names(clim_raster)
  layer_names <- all_layers[seq(1, length(all_layers), 24)] # Keep every 24th layer name (1 per day)


  ## Aggregate to grid-day level
  ## -----------------------------------------------

  if(!daily_agg %in% c('average', 'sum')){

    stop(crayon::red("Daily aggregation must be 'average' or 'sum'"))
  }

  ## Average
  if(daily_agg == 'average'){

    message(crayon::green("Averaging hourly values to get daily values"))

    # Average over each set of 24 layers
    indices<-rep(1:(raster::nlayers(clim_raster)/24),each=24)
    clim_daily <- raster::stackApply(clim_raster, indices = indices, fun=mean) #Stack of 365/366 layers
  }

  ## Sum
  if(daily_agg == 'sum'){

    message(crayon::green("Summing hourly values to daily values"))

    # Sum over each set of 24 layers
    indices<-rep(1:(raster::nlayers(clim_raster)/24),each=24)
    clim_daily <- raster::stackApply(clim_raster, indices = indices, fun=sum) #Stack of 365/366 layers
  }



  # Return a list containing, in order, daily aggregated climate data as a raster brick, and the layer_names created.
  return(list(clim_daily, layer_names))
}


# Function to merge with geoweights and aggregate by polygon
polygon_aggregation <- function(clim_dt, weights_dt, list_names, time_agg){

  ## Merge weights with climate raster
  ## -----------------------------------------------

  # Set key column in the climate data.table
  keycols = c("x", "y")
  data.table::setkeyv(clim_dt, keycols)


  # Convert layer names to dates
  message(crayon::yellow("Assuming layer name format which after removal of the first character is compatible with lubridate::as_date()"))
  clim_dt[, date := stringr::str_sub(date, 2, -1)]
  clim_dt[, date := lubridate::as_date(date)]

  # Keyed merge on the x/y column
  merged_dt <- clim_dt[weights_dt, allow.cartesian = TRUE] #cols: x, y, date, value cols 1:k, poly_id, w_area, weight (if weights = T)

  ## Multiply weights x climate value (all 1:k values); aggregate by month and polygon
  ## -----------------------------------------------

  # Multiply by secondary weights if weights = TRUE (already normalized by polygon area)
  # Otherwise multiply by just area weights
  if("weight" %in% names(merged_dt)){
    merged_dt[, (list_names) := lapply(list_names, function(x) {get(x) * weight})]
  } else {
    merged_dt[, (list_names) := lapply(list_names, function(x) {get(x) * w_area})]
  }

  # Separate year, month, and day columns
  merged_dt[, ':=' (year = lubridate::year(date),
                    month = lubridate::month(date),
                    day = lubridate::day(date))]

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



  ## Return the sums by polygon
  return(sum_by_poly)



}



#=====================================================================================================================================================

#' Polynomial Transformation and Aggregation of Climate Data
#'
#' @param data The raster brick with the data to be transformed and aggregated
#' @param overlay_weights A table of weights which can be generated using
#'   the function calc_geoweights()
#' @param daily_agg How to aggregate daily values ('sum' and 'average' currently
#'   supported)
#' @param time_agg the temporal scale to aggregate data to ('day', 'month', and
#'   'year' currently supported)
#' @param degree the highest exponent to raise the data to
#'
#' @examples
#' polynomial_output <- staggregate_polynomial(
#'
#'   data = prcp_kansas_dec2011_era5, # Climate data to transform and aggregate
#'
#'   overlay_weights = overlay_weights_kansas, # Output from overlay_weights()
#'
#'   daily_agg = "sum", # Sum hourly values to produce daily values before transformation
#'
#'   degree = 4 # Highest order
#'   )
#'
#'
#' head(polynomial_output)
#'
#' @export
staggregate_polynomial <- function(data, overlay_weights, daily_agg, time_agg = "month", degree){

  # Get climate data as a data.table and aggregate to daily values
  setup_list <- daily_aggregation(data, overlay_weights, daily_agg)

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

  # Merge all data.tables together
  clim_dt <- list_dt[[1]]
  for(i in 2:list_length){
    dt_m <- list_dt[[i]]
    clim_dt <- merge(clim_dt, dt_m, by=c('x', 'y', 'date'))
  }

  # Aggregate by polygon
  sum_by_poly <- polygon_aggregation(clim_dt, overlay_weights, list_names, time_agg)

  return(sum_by_poly)

}









#==================================================================================================================================================

#' Restricted Cubic Spline Transformation and Aggregation of Climate Data
#'
#' @param data The raster brick with the data to be transformed and aggregated
#' @param overlay_weights A table of weights which can be generated using
#'   the function calc_geoweights()
#' @param daily_agg How to aggregate daily values ('sum' and 'average' currently
#'   supported)
#' @param time_agg the temporal scale to aggregate data to ('day', 'month', and
#'   'year' currently supported)
#'   'polynomial', this is an integer indicating the degree.
#' @param knot_locs where to place the knots
#'
#' @examples

#' spline_output <- staggregate_spline(
#'
#' data = prcp_kansas_dec2011_era5, # Climate data to transform and aggregate
#'
#' overlay_weights = overlay_weights_kansas, # Output from overlay_weights()
#'
#' daily_agg = "sum", # Sum hourly values to produce daily values before transformation
#'
#' knot_locs = c(-1.665335e-16,  1.129775e-06,  1.594356e-02) # Where to place knots
#' )
#'
#'
#' head(spline_output)
#'
#' @export
staggregate_spline <- function(data, overlay_weights, daily_agg, time_agg = "month", knot_locs){

  # Get climate data as a data.table and aggregate to daily values
  setup_list <- daily_aggregation(data, overlay_weights, daily_agg)

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

  i <- 2
  while(i <= list_length + 1){
    dt_m <- list_dt[[i]]
    clim_dt <- merge(clim_dt, dt_m, by=c('x', 'y', 'date'))
    i <- i + 1
  }


  # Aggregate by polygon
  sum_by_poly <- polygon_aggregation(clim_dt, overlay_weights, list_names, time_agg)

  return(sum_by_poly)
}

#==================================================================================================================================================

#' Bin Transformation and Aggregation of Climate Data
#'
#' @param data The raster brick with the data to be transformed and aggregated
#' @param overlay_weights A table of weights which can be generated using the
#'   function calc_geoweights()
#' @param daily_agg How to aggregate daily values ('sum' and 'average' currently
#'   supported) 'polynomial', this is an integer indicating the degree.
#' @param time_agg the temporal scale to aggregate data to ('day', 'month', and
#'   'year' currently supported)
#' @param num_bins number of bins to group the data into
#' @param binwidth width of bins, overrides num_bins
#' @param min the smallest value that must be captured by a non-edge bin,
#'   default is the data minimum, set manually if you are breaking up your data
#' @param max the largest value that must be captured by a non-edge bin, default
#'   is the data maximum, set manually if you are breaking up your data
#' @param start_on where to place the left edge of one of the bins
#' @param center_on where to place the center of one of the bins
#' @param end_on where to place the right edge of one of the bins
#'
#' @examples
#' bin_output <- staggregate_bin(
#'
#'   data = prcp_kansas_dec2011_era5, # Climate data to transform and aggregate
#'
#'   overlay_weights = overlay_weights_kansas, # Output from overlay_weights()
#'
#'   daily_agg = "sum", # Sum hourly values to produce daily values before transformation
#'
#'   bin_breaks = c(0, 2, 4, 6, 8, 10) # Draw 5 bins from 0 to 10, a bin from -inf to 0, and one from 10 to inf
#'   )
#'
#'
#' head(bin_output)
#'
#' @export
staggregate_bin <- function(data, overlay_weights, daily_agg, time_agg = "month", bin_breaks){

  # Get climate data as a data.table and aggregate to daily values
  setup_list <- daily_aggregation(data, overlay_weights, daily_agg)
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
