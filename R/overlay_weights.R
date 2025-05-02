#' Find the spatial overlap between a raster and a set of polygons
#'
#' The `overlay_weights()` function generates a table of weights mapping
#' each grid cell to its respective polygon(s) for use in the `staggregate_*()`
#' family of functions.
#'
#' @param polygons a simple features polygon or multipolygon object
#' @param polygon_id_col the name of a column in the `polygons` object with a
#'   unique identifier for each polygon
#' @param grid a raster layer with the same spatial resolution as the data
#'   (must use geographic coordinates)
#' @param secondary_weights an optional table of secondary weights, output from
#'   the `secondary_weights()` function
#'
#' @return a data.table of area weights and possibly secondary weights for each
#'   cell within each polygon
#'
#' @examples
#' overlay_output_with_secondary_weights <- overlay_weights(
#'   polygons = tigris::counties("nj"), # Polygons outlining the 21 counties of New Jersey
#'   polygon_id_col = "COUNTYFP", # The name of the column with the unique
#'                                # county identifiers
#'   grid = era5_grid, # The grid to use when extracting area weights (era5_grid is the
#'                     # default)
#'   secondary_weights = cropland_world_2015_era5 # Output from
#'                                                # secondary_weights
#'                                                # (cropland_world_2015_era5 is
#'                                                # available to the# user)
#'   )
#'
#' head(overlay_output_with_secondary_weights)
#'
#'
#'
#' overlay_output_without_secondary_weights <- overlay_weights(
#'   polygons = tigris::counties("nj"), # Polygons outlining the 21 counties of New Jersey
#'   polygon_id_col = "COUNTYFP" # The name of the column with the unique county
#'                               # identifiers
#'   )
#'
#' head(overlay_output_without_secondary_weights)
#'
#'
#' @export
overlay_weights <- function(polygons, polygon_id_col, grid = era5_grid, secondary_weights = NULL){

  ## check to make sure climate raster is a spatraster, change if not
  if (!inherits(grid, "SpatRaster")) {
    clim_raster <- terra::rast(grid)
  } else {

    clim_raster <- grid

  }

  ## Raster cell area
  ## -----------------------------------------------
  clim_area_raster <- terra::cellSize(clim_raster, unit = "km")

  ## Raster/polygon alignment
  ## -----------------------------------------------

  message(crayon::yellow('Checking for raster/polygon alignment'))

  ## polygon and raster xmin and xmax values
  poly_xmax <- terra::ext(polygons)$xmax
  rast_xmax <- terra::ext(clim_area_raster)$xmax
  rast_res <-  terra::xres(clim_area_raster)

  ## check if SpatRaster is in geographic coodrinates
  if(!terra::is.lonlat(clim_raster)) {
    stop(crayon::red('Grid does not have geographic coordinates.'))

  }

 ## stop if polygons are not in standard coordinate system
 if(poly_xmax > 180) {

   stop(crayon::red('Polygons must be in standard coordinate system (longitude -180 to 180).'))

 }

  ## check if coordinate systems match, if no shift raster to -180 to 180
  if(rast_xmax > 180 + rast_res) {

    # Make sure the cell widths aren't peculiar otherwise the rotate function will
    # mess things up
    if(360 %% terra::xres(clim_raster) != 0){
      stop(crayon::red('Grid is in climate coordinate system (longitude 0 to 360) and grid cell width does not divide 360 evenly, making accurate alignment impossible.'))
    }

    message(crayon::yellow('Aligning longitudes to standard coordinates.'))

    ## xmin for climate raster
    rast_xmin <- terra::ext(clim_area_raster)$xmin

    ## check if raster needs to be padded, extend if needed
    if(!dplyr::near(rast_xmin, 0, tol = rast_res) | !dplyr::near(rast_xmax, 360, tol = rast_res)) {

      ## create global extent for padding so rotate function can be used
      global_extent <- terra::ext(0, 360, -90, 90)

      ## pad
      clim_area_raster <- terra::extend(clim_area_raster, global_extent)

    }

    ## rotate
    clim_area_raster <- terra::rotate(clim_area_raster)

  }

  # Extend the grid to cover all polygons consistent with the extended rotate
  # above. exact_extract already does this in the background to a certain
  # degree, so this just allows us to be explicit about how we handle NAs later
  # on.
  clim_area_raster <- terra::extend(clim_area_raster, terra::ext(polygons), snap = 'out')

  ## Match raster and polygon crs
  crs_raster <- terra::crs(clim_area_raster)
  polygons_reproj <- sf::st_transform(polygons, crs = crs_raster)

  ## Raster / Polygon overlap (using data.table)
  ## -----------------------------------------------
  message(crayon::green('Extracting raster polygon overlap'))

  overlap <- data.table::rbindlist(exactextractr::exact_extract(clim_area_raster, polygons_reproj, progress = F, include_xy = T), idcol = "poly_id")
  overlap[, ':=' (poly_id = polygons_reproj[[polygon_id_col]][poly_id], cell_area_km2 = value)] # Add the unique id for each polygon based on the input col name


  ## Calculate weights
  ## -----------------------------------------------

  # Calculate area weight per grid cell
  area_weight <- overlap[, .(x, y, poly_id, w_area = coverage_fraction * cell_area_km2)] # area weight = area km2 * coverage fraction

  # IF weights = TRUE, merge secondary weights with area weights
  if(!is.null(secondary_weights)){

    # data.table of secondary weights
    weights_dt <- data.table::as.data.table(secondary_weights)

    ## make sure secondary_weights is not in climate 0-360 coordinates
    s_weight_max <- max(weights_dt$x)

    ## if secondary_weights is in 0-360, adjust x val
    if(s_weight_max > 180 + rast_res / 2) {

      message(crayon::yellow('Adjusting secondary weights longitude to standard coordinates.'))

      weights_dt[, x := data.table::fifelse(x > 180 + rast_res / 2, x - 360, x)]

    }

    ## check if secondary weights fully overlaps with area_weight df
    covers <- min(weights_dt$x) <= min(area_weight$x) &&
              max(weights_dt$x) >= max(area_weight$x) &&
              min(weights_dt$y) <= min(area_weight$y) &&
              max(weights_dt$y) >= max(area_weight$y)

    if (covers) {
      message(crayon::green('Secondary weights fully overlap with the administrative regions.'))
    } else {
      warning(crayon::red('Warning: secondary weights do not fully overlap with the administrative regions. Resulting weights will contain NAs.'))
    }


  ## check if secondary weights table contains NA values
  if(isTRUE(any(is.na(weights_dt[["weight"]])))) {

    ## print warning if there are NAs in the secondary weights
    warning(crayon::red("Warning: secondary weight values contain one or more NAs. The resulting weights for x,y coordinates with NA secondary weight values will be NAs."))

  }

    # Set key column in the merged dt table
    keycols = c("x", "y")
    data.table::setkeyv(area_weight, keycols)

    # Merge with secondary weights, NA for missing values
    w_merged <- merge(area_weight, weights_dt,
                      by = c('x', 'y'),
                      all.x = T)


    # Weight in pixel = w_area * weight
    w_merged[, weight := weight * w_area]

    # Create column that determines if entire polygon has a weight == 0
    zero_polys <- w_merged[, .(sum_weight = sum(weight)),
                                by = .(poly_id)]

    zero_polys <- unique(zero_polys[sum_weight == 0, .(poly_id)])

    if(nrow(zero_polys) > 0) {

      warning(crayon::red("Warning: weight = 0 for all pixels in some of your polygons; NAs will be returned for weights"))

    }

    # List any polygons with NA values in 1 or more grid cells
    na_polys <- unique(w_merged[is.na(weight), .(poly_id)])

    # # Warning if there are polygons with NA weight values
    # if(nrow(na_polys > 0)) {
    #
    #   warning(crayon::red("Warning: some of the secondary weights are NA, meaning weights cannot be calculated. NAs will be returned for weights."))
    #
    # }

    # Update the weight to NA for all grid cells in na_polys
    w_merged <- w_merged[, weight := ifelse(poly_id %in% c(na_polys$poly_id,
                                                           zero_polys$poly_id), NA, weight)]

  }

  # Normalize weights by polygon
  if(!is.null(secondary_weights)){

    w_norm <- data.table::copy(w_merged)

    w_norm <- w_norm[, ':=' (w_area = w_area / sum(w_area), weight = weight / sum(weight)), by = poly_id]


  } else {

    w_norm <- data.table::copy(area_weight)

    w_norm <- w_norm[, w_area := w_area / sum(w_area), by = poly_id]
  }


  message(crayon::yellow('Checking sum of weights within polygons'))
  if(!is.null(secondary_weights)){

    check_weights <- w_norm[, lapply(.SD, sum, na.rm = T), by = poly_id,
                            .SDcols = c('w_area', 'weight')]

  } else {
    check_weights <- w_norm[, w_sum := sum(w_area), by=poly_id]
  }

  # Check that polygon weights sum to 1 or 0 if all weights are NA
  if (!is.null(secondary_weights)){

    for(i in seq_len(nrow(check_weights))){

      if(!dplyr::near(check_weights$w_area[i], 1, tol=0.001)) {

        stop(crayon::red('Area weights for polygon', check_weights$poly_id[i], 'do not sum to 1')) }

        if(!check_weights$poly_id[i] %in% c(na_polys$poly_id, zero_polys$poly_id) & !dplyr::near(check_weights$weight[i], 1, tol=0.001)) {

          stop(crayon::red('Weights for polygon', check_weights$poly_id[i], 'do not sum to 1')) }

     }


    } else {

    for(i in seq_len(nrow(check_weights))){

      if(!dplyr::near(check_weights$w_sum[i], 1, tol=0.001)){

        stop(crayon::red('Area weights for polygon', check_weights$poly_id, 'do not sum to 1'))

      }

    }
  }



  # If it doesn't error out then all weight sums = 1
  message(crayon::green('All weights sum to 1.'))

  ## Return table in coordinate system that matches that of the climate data
  ## ------------------------------------------------------------------------

  if(rast_xmax > 180 + rast_res) {

    w_norm[, x := data.table::fifelse(x < 0 + rast_res / 2, x + 360, x)]

  }


  return(w_norm)

}

