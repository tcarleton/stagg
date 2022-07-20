#' Aggregate Climate Data
#'
#' @param year The year to aggregate
#' @param data_raster The raster containing climate data
#' @param climate_var The  climate variable of interest ('prcp, 'temp', and 'uv'
#'   currently supported)
#' @param daily_agg How to aggregate daily values ('sum' and 'average' currently
#'   supported)
#' @param trans A statistical transformation ('polynomial' currently supported)
#' @param trans_specs Specifications pertaining to selected transformation. For
#'   'polynomial', this is an integer indicating the degree.
#' @param geoweights_table A table of weights which can be generated using
#'   the function calc_geoweights()
#' @param second_weights A boolean indicating whether secondary weights should be used
#'
#' @examples
#' no_cores <- parallel::detectCores() - 1 # Calculate the number of cores. Leave one in case something else needs to be done on the same computer at the same time.
#' cl <- parallel::makeCluster(no_cores, type="FORK") # Initiate cluster. "FORK" means bring everything in your current environment with you.
#' sum_by_poly_multiyear <- parallel::parLapply(cl, years, agg_climate_data, data_raster, climate_var, daily_agg, trans, trans_specs, geoweights_table)
#' parallel::stopCluster(cl)
#'
#' bin_output <- staggregate_bin(demo_prcp, "prcp", daily_agg = "sum", num_bins = 30, min = 5, year = 2011, geoweights_table = demo_output_geoweights, second_weights = TRUE)
#'
#'

## Sara Orofino
## February 22, 2022
## Aggregate Climate Data: Pipeline Step 02

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

# Function to convert raster to data table and aggregate to daily values before transformation
before_trans <- function(data_raster, climate_var, daily_agg, geoweights_table, second_weights = FALSE){


  # Data.table of weights
  weights_dt <- geoweights_table

  # Extent of area weights with slight buffer to make sure all cells are included
  min_x <- min(weights_dt$x) - 0.5
  max_x <- max(weights_dt$x) + 0.5
  min_y <- min(weights_dt$y) - 0.5
  max_y <- max(weights_dt$y) + 0.5

  weights_ext <- raster::extent(min_x, max_x, min_y, max_y)

  ## Load climate data
  ## -----------------------------------------------

  # Check if the climate variable is one we have data for
  if(!climate_var %in% c('prcp', 'temp', 'uv')){
    stop(crayon::red('No ERA5 data available. Supported variables are: prcp, temp, or uv'))
  }


  # Immediately crop to weights extent
  clim_raster <- raster::crop(raster::stack(data_raster), weights_ext)

  # Get layer names (dates)
  all_layers <- names(clim_raster)
  layer_names <- all_layers[seq(1, length(all_layers), 24)] # Keep every 24th layer name (1 per day)
  layer_names <- paste0('month_', substring(layer_names, 7,8), '_day_', substring(layer_names, 10,11)) # Extract month/day for each layer

  ## Aggregate to grid-day level
  ## -----------------------------------------------

  if(!daily_agg %in% c('average', 'sum')){

    stop(crayon::red("Daily aggregation must be 'average' or 'sum'"))
  }

  ## Average
  if(daily_agg == 'average'){

    # Check if there are 24*365 or 24*366 layers
    if(!raster::nlayers(clim_raster) %in% c(8760, 8784)){

      message(crayon::red("Warning: Incomplete year of data; raster has", length(all_layers),
                          "layers, but a complete year should have 8760 layers or 8784 layers on a leap year"))
    }

    # Average over each set of 24 layers
    indices<-rep(1:(raster::nlayers(clim_raster)/24),each=24)
    clim_daily <- raster::stackApply(clim_raster, indices = indices, fun=mean) #Stack of 365/366 layers
  }

  ## Sum
  if(daily_agg == 'sum'){

    # Check if there are 24*365 or 24*366 layers
    if(!raster::nlayers(clim_raster) %in% c(8760, 8784)){

      message(crayon::red("Warning: Incomplete year of data; raster has", length(all_layers),
                          "layers, but a complete year should have 8760 layers or 8784 layers on a leap year"))
    }

    # Sum over each set of 24 layers
    indices<-rep(1:(raster::nlayers(clim_raster)/24),each=24)
    clim_daily <- raster::stackApply(clim_raster, indices = indices, fun=sum) #Stack of 365/366 layers
  }

  # For temperature convert values in Kelvin to Celsius C = K - 273.15
  if(climate_var == 'temp'){

    clim_daily <- clim_daily - 273.15
  }

  ## convert prcp m to mm
  if(climate_var == 'prcp'){

    clim_daily <- clim_daily * 1000
  }

  # Return a list of the necessary objects
  return(list(clim_daily, weights_dt, layer_names))
}


# Function to merge with geoweights and aggregate by polygon
after_trans <- function(clim_dt, weights_dt, list_names, second_weights, year){

  ## Merge weights with climate raster
  ## -----------------------------------------------

  # Set key column in the climate data table
  keycols = c("x", "y")
  data.table::setkeyv(clim_dt, keycols)

  # Keyed merge on the x/y column
  merged_dt <- clim_dt[weights_dt, allow.cartesian = TRUE] #cols: x, y, date, value cols 1:k, poly_id, w_area, weight (if weights = T)

  ## Multiply weights x climate value (all 1:k values); aggregate by month and polygon
  ## -----------------------------------------------

  # Multiply by secondary weights if weights = TRUE (already normalized by polygon area)
  # Otherwise multiply by just area weights
  if(second_weights){
    merged_dt[, (list_names) := lapply(list_names, function(x) {get(x) * weight})]
  } else {
    merged_dt[, (list_names) := lapply(list_names, function(x) {get(x) * w_area})]
  }

  # Separate month and day columns
  merged_dt[, ':=' (month = substring(date, first=1, last=8),
                    day = substring(date, first=10))]

  # Can customize this in the future to aggregate by day & month
  # Right now just sum by month
  sum_by_poly <- merged_dt[,  lapply(.SD, sum), by = .(poly_id, month),
                           .SDcols = list_names]

  ## Add year column
  sum_by_poly[, year := year]

  ## Order columns
  data.table::setcolorder(sum_by_poly, neworder = c('year', 'month', 'poly_id', list_names))

  ## Return the sums by polygon
  return(sum_by_poly)



}

#' @export
staggregate_polynomial <- function(data_raster, climate_var, daily_agg, k, geoweights_table, second_weights, year){

  # Get climate data as a data table and aggregate to daily values
  setup_list <- before_trans(data_raster, climate_var, daily_agg, geoweights_table, second_weights)

  clim_daily <- setup_list[[1]]
  weights_dt <- setup_list[[2]]
  layer_names <- setup_list[[3]]

  # Polynomial transformation
  poly_orders <- seq(1:k) # Compute values from 1 to K
  list_length <- length(poly_orders) # How many lists are in the final object
  list_names <- sapply(1:list_length, FUN=function(x){paste("order", poly_orders[x], sep="_")})

  # For each daily layer, raise the value to k, k-1, k-2 etc. until 1
  r <- lapply(poly_orders, FUN=function(x){clim_daily ^ x})


    ## Function: Set names of data table by month, change from wide to long format, rename based on polynomial orders
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

  # Merge all data tables together
  clim_dt <- list_dt[[1]]
  for(i in 2:list_length){
    dt_m <- list_dt[[i]]
    clim_dt <- merge(clim_dt, dt_m, by=c('x', 'y', 'date'))
  }

  # Aggregate by polygon
  sum_by_poly <- after_trans(clim_dt, weights_dt, list_names, second_weights, year)

  return(sum_by_poly)

}

#' @export
staggregate_spline <- function(data_raster, climate_var, daily_agg, knot_locs, geoweights_table, second_weights, year){

  # Get climate data as a data table and aggregate to daily values
  setup_list <- before_trans(data_raster, climate_var, daily_agg, geoweights_table, second_weights)

  clim_daily <- setup_list[[1]]
  weights_dt <- setup_list[[2]]
  layer_names <- setup_list[[3]]


  # Spline transformation
  knot_locs <- sort(knot_locs)
  num_knots <- length(knot_locs)
  list_length <- num_knots - 2
  list_names <- sapply(0:list_length, FUN=function(x){if(x == 0){paste(climate_var)}else{paste("term", x, sep="_")}})


  # Define spline function
  get_spline <- function(x){

    # Make first raster returned just the climate variable to preserve it's column in the resulting data table
    if(x == 0){
      return(clim_daily)
    }
    # Add in spline terms, all of which are 0 if negative
    else{
      clim_daily_table <- raster::values(clim_daily)

      clim_daily_table <-ifelse((clim_daily_table - knot_locs[x])^3 > 0,
                                (clim_daily_table - knot_locs[x])^3, 0)
      - ((ifelse((clim_daily_table - knot_locs[num_knots - 1])^3 > 0,
                 (clim_daily_table - knot_locs[num_knots - 1])^3, 0))
         *((knot_locs[num_knots] - knot_locs[x]))
           / (knot_locs[num_knots] - knot_locs[num_knots - 1]))
      + ((ifelse((clim_daily_table - knot_locs[num_knots])^3 > 0,
                 (clim_daily_table - knot_locs[num_knots])^3, 0))
         *((knot_locs[num_knots - 1] - knot_locs[x])
           / (knot_locs[num_knots] - knot_locs[num_knots - 1])))

      clim_daily_new <- clim_daily
      raster::values(clim_daily_new) <- clim_daily_table

      return(clim_daily_new)



    }
  }


  # For each layer, create new spline variables
  r <- lapply(0:list_length, get_spline)


  ## Function: Set names of data table by month, change from wide to long format, rename based on polynomial orders
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

  # Merge all data tables together
  clim_dt <- list_dt[[1]]

  i <- 2
  while(i <= list_length + 1){
    dt_m <- list_dt[[i]]
    clim_dt <- merge(clim_dt, dt_m, by=c('x', 'y', 'date'))
    i <- i + 1
  }



  # Aggregate by polygon
  sum_by_poly <- after_trans(clim_dt, weights_dt, list_names, second_weights, year)

  return(sum_by_poly)
}



#' @export
staggregate_bin <- function(data_raster, climate_var, daily_agg, num_bins = 5, binwidth = NULL, center_on = NULL, start_on = NULL, end_on = NULL, max = NULL, min = NULL, geoweights_table, second_weights, year){
  # Get climate data as a data table and aggregate to daily values
  setup_list <- before_trans(data_raster, climate_var, daily_agg, geoweights_table, second_weights)
  clim_daily <- setup_list[[1]]
  weights_dt <- setup_list[[2]]
  layer_names <- setup_list[[3]]


  clim_daily_table <- values(clim_daily)

  if(is.null(max)){
    max <- max(clim_daily_table) + .01
  }
  if(is.null(min)){
    min <- min(clim_daily_table) - .01
  }

  # Write message that binwidth overrides num_bins
  if(!is.null(binwidth)){
    message(crayon::yellow("Binwidth argument supplied will override num_bins. Bins at the boundaries may extend beyond the maximum and minimum values of the data supplied"))
  }

  # If value not supplied to binwidth, calculate from num_bins
  if(is.null(binwidth)){
    binwidth <- (max - min) / num_bins
  }

  # Stop with message that only one of center_on, start_on, and end_on can be chosen if otherwise
  if((!is.null(center_on) & !is.null(start_on)) |
     (!is.null(center_on) & !is.null(end_on)) |
     (!is.null(start_on) & !is.null(end_on))){
    stop(crayon::red("Too many bin placement arguments specified. Please specify a value for only one variable out of center_on, start_on, or end_on"))
  }

  if(is.null(center_on) & is.null(start_on) & is.null(end_on)){
    message(crayon::yellow("No bin placement argument specified. Drawing bins from min value."))
    start_on <- min
  }

  # Put start_on and end_on in terms of center_on
  if(!is.null(start_on)){
    center_on <- start_on + (binwidth / 2)
  }

  if(!is.null(end_on)){
    center_on <- end_on - (binwidth / 2)
  }



  # Check to see that center_on falls within the range of data, move it if not
  if(max - min < binwidth){
    stop(crayon::red("The binwidth is larger than the range of data. Please specify a different binwidth and rerun"))
  }

  if(center_on > max){
    message(crayon::yellow("Bin center is greater than maximum value in data. Adjusting it by a multiple of the binwidth so that the bin center falls within the range of data"))

    while(center_on > max){
      center_on <- center_on - binwidth
    }
  }
  if(center_on < min){
    message(crayon::yellow("Bin center is less than minimum value in data. Adjusting it by a multiple of the binwidth so that the bin center falls within the range of data"))

    while(center_on < min){
      center_on <- center_on + binwidth
    }
  }


  # Create table of bins
  center <- c(center_on)
  while(max(center) + (binwidth / 2) < max){
    center <- c(center, (max(center) + binwidth))
  }
  while(min(center) - (binwidth / 2) > min){
    center <- c(center, (min(center)) - binwidth)
  }

  center <- sort(center)



  bins_table <- data.table::data.table(center)
  bins_table[, ':=' (start = center - (binwidth / 2), end = center + (binwidth / 2))]

  #To deal with potential overlap, give warning
  if(nrow(bins_table) > num_bins){
    warning(crayon::red("Bins drawn extend beyond data, number of bins may be greater than previously specified by one"))
    num_bins <- nrow(bins_table)
  }


  # bins_table lists center, start, and end of all bins in order

  # Function check_bins to determine which bins data points fall into
  check_bins <- function(x){
    clim_daily_table <- raster::values(clim_daily)

    if(x == 0){
      clim_daily_table <- ifelse(min > clim_daily_table, 1, 0)
    }
    else if(x == num_bins + 1){
      clim_daily_table <- ifelse(max < clim_daily_table, 1, 0)
    }
    else{
      clim_daily_table <- ifelse(bins_table[x, start] <= clim_daily_table &
                                   bins_table[x, end] >= clim_daily_table, 1, 0)
    }


    clim_daily_new <- clim_daily
    raster::values(clim_daily_new) <- clim_daily_table

    return(clim_daily_new)
  }


  # Create names for new columns
  list_names <- sapply(0:(num_bins + 1), FUN=function(x){
    if(x == 0){
      paste("-inf", "to", min, sep = "_")
    }
    else if(x == num_bins + 1){
      paste(max, "to", "inf", sep = "_")
    }
    else{
      paste(bins_table[x, start], "to", bins_table[x, end], sep = "_")
      }
    })

  # For each bin, create new brick of binary values, including end bins which go from -inf to min, max to inf
  r <- lapply(0:(num_bins + 1), FUN = check_bins)


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
  list_dt <- lapply(1:(num_bins + 2), create_dt)

  # Merge all data tables together
  clim_dt <- list_dt[[1]]
  for(i in 2:(num_bins + 2)){
    dt_m <- list_dt[[i]]
    clim_dt <- merge(clim_dt, dt_m, by=c('x', 'y', 'date'))
  }



  # Aggregate by polygon
  sum_by_poly <- after_trans(clim_dt, weights_dt, list_names, second_weights, year)

  return(sum_by_poly)
}



# Create end bins from -inf to min and max to inf
