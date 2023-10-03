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
#' @param secondary_weights an optional table of secondary weights, output from
#'   the `secondary_weights()` function
#'
#' @return a data.table of area weights and possibly secondary weights for each
#'   cell within each polygon
#'
#' @examples
#' kansas_counties <- tigris::counties("Kansas")
#'
#' overlay_output_with_secondary_weights <- overlay_weights(
#'   polygons = kansas_counties, # Polygons outlining the 105 counties of Kansas
#'   polygon_id_col = "COUNTYFP", # The name of the column with the unique
#'                                # county identifiers
#'   era5_grid, # The grid to use when extracting area weights (era5_grid is the
#'              # default)
#'   cropland_world_2015_era5 # Output from secondary_weights
#'                            # (cropland_world_2015_era5 is available to the
#'                            # user)
#'   )
#'
#' head(overlay_output_with_secondary_weights)
#'
#'
#'
#' overlay_output_without_secondary_weights <- overlay_weights(
#'   polygons = kansas_counties, # Polygons outlining the 105 counties of Kansas
#'   polygon_id_col = "COUNTYFP" # The name of the column with the unique county
#'                               # identifiers
#'   )
#'
#' head(overlay_output_without_secondary_weights)
#'
#'
#' @export
overlay_weights <- function(polygons, polygon_id_col, grid = era5_grid, secondary_weights = NULL){

  # Create raster
  clim_raster <- raster::raster(grid) # only reads the first band

  ## make sure climate raster is in 0-360
  c_rast_xmin <- raster::extent(clim_raster)@xmin
  c_rast_res <- raster::xres(clim_raster)

  ## if climate raster is not in 0-360, shift
  clim_raster <- if(c_rast_xmin >= 0 - c_rast_res / 2) {
    clim_raster
  } else {
    raster::shift(clim_raster, dx = 360)}

  ## Raster cell area
  ## -----------------------------------------------
  clim_area_raster <- suppressWarnings(raster::area(clim_raster)) # Suppressing warning that is
                                                                  # thrown because era5_grid has
                                                                  # lat extent <-90.1 and >90.1

  ## Raster/polygon alignment
  ## -----------------------------------------------

  message(crayon::yellow('Checking for raster/polygon alignment'))

  ## polygon and raster xmin and xmax values
  poly_xmin <- raster::extent(polygons)@xmin
  # poly_xmax <- raster::extent(polygons)@xmax
  # rast_xmin <- raster::extent(clim_area_raster)@xmin
  # rast_xmax <- raster::extent(clim_area_raster)@xmax
  rast_res <-  raster::xres(clim_area_raster)


 ## do the polygons need to be in standard coord system? possibly delete
  # ## stop if polygons are not in standard coordinate system
  # if(poly_xmin < 0 - c_rast_res / 2) {
  #
  #   stop(crayon::red('Polygons must be in standard coordinate system (longitude -180 to 180).'))
  #
  # }


  ## check if coordinate systems match, if no shift raster in 0-360 format
  if(poly_xmin < 0 - rast_res / 2) {

    message(crayon::yellow('Aligning longitudes to 0-360 coordinates.'))

    polygons <- sf::st_shift_longitude(polygons)

     } else {

    message(crayon::green('Coordinate systems match.'))

  }


  # Match raster and polygon crs
  crs_raster <- raster::crs(clim_area_raster)
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

    # Data.table of secondary weights
    weights_dt <- data.table::as.data.table(secondary_weights)

    # xmin of secondary weights
    weights_xmin <- min(weights_dt$x)

    # If weights in standard coordinates (-180 to 180), convert
    if(weights_xmin < 0) {

      message(crayon::yellow('Adjusting secondary weights longitude 0-360'))

      weights_dt[, x := data.table::fifelse(x < 0, x + 360, x)]

    } else {
      message(crayon::green('No need to adjust secondary weights'))}

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
    zero_polys <- data.frame(w_merged) %>%
      dplyr::group_by(poly_id) %>%
      dplyr::summarise(sum_weight = sum(weight)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(sum_weight == 0) %>%
      dplyr::select(poly_id) %>%
      dplyr::distinct()

    if(nrow(zero_polys > 0)) {

      warning(crayon::red("Warning: weight = 0 for all pixels in some of your polygons; NAs will be returned for weights"))

    }

    # List any polygons with NA values in 1 or more grid cells
    na_polys <- data.frame(w_merged) %>%
      dplyr::filter(is.na(weight)) %>%
      dplyr::select(poly_id) %>%
      dplyr::distinct()


    # Warning if there are polygons with NA weight values
    if(nrow(na_polys > 0)) {

      warning(crayon::red("Warning: some of the secondary weights are NA, meaning weights cannot be calculated. NAs will be returned for weights."))

    }

    # Update the weight to equal w_area for all grid cells in na_polys
    w_merged <- w_merged %>%
      dplyr::mutate(weight = ifelse(poly_id %in% c(na_polys$poly_id, zero_polys$poly_id), NA, weight)) %>%
      data.table::as.data.table()

    ## this doesn't work with dt... figure out or delete and use above
    # w_merged[, weight := data.table::fifelse(poly_id %in% c(na_polys$poly_id, zero_polys$poly_id), NA, weight)]

  }

  # Normalize weights by polygon
  if(!is.null(secondary_weights)){

    w_norm <- data.table::copy(w_merged)

    w_norm <- w_norm[, ':=' (w_area = w_area / sum(w_area), weight = weight / sum(weight)), by = poly_id]


  } else {

    w_norm <- data.table::copy(area_weight)

    w_norm <- area_weight[, w_area := w_area / sum(w_area), by = poly_id]
  }


  message(crayon::yellow('Checking sum of weights within polygons'))
  if(!is.null(secondary_weights)){

    check_weights <- w_norm[, lapply(.SD, sum), by = poly_id,
                            .SDcols = c('w_area', 'weight')]

  } else {
    check_weights <- w_norm[, w_sum := sum(w_area), by=poly_id]
  }

  # Check that polygon weights sum to 1
  if (!is.null(secondary_weights)){

    for(i in nrow(check_weights)){

      if(!dplyr::near(check_weights$w_area[i], 1, tol=0.001)) {

        stop(crayon::red('Weights for polygon', check_weights$poly_id, 'do not sum to 1')) }

        if(!is.na(check_weights$weight[i]) & !dplyr::near(check_weights$weight[i], 1, tol=0.001)) {

          stop(crayon::red('Weights for polygon', check_weights$poly_id, 'do not sum to 1 and are not NA'))

        }

     }


    } else {

    for(i in nrow(check_weights)){

      if(!dplyr::near(check_weights$w_sum[i], 1, tol=0.001)){

        stop(crayon::red('Area weights for polygon', check_weights$poly_id, 'do not sum to 1'))

      }

    }
  }



  # If it doesn't error out then all weight sums = 1
  message(crayon::green('All weights sum to 1.'))

  ## delete
  # ## Convert back to 0-360
  # ## -----------------------------------------------
  #
  # w_norm[, x := data.table::fifelse(x < 0, x + 360, x)]

  return(w_norm)

}

