## Sara Orofino
## January 31, 2022
## Spatial Overlap: Pipeline Step 01


#' Function to find spatial overlap between a raster and a set of polygons
#'
#' @param data_source the source of climate data (default is era5)
#' @param input_polygons a simple features polygon or multipolygon object
#' @param polygon_id the name of a column in the sf object representing a unique
#'   identifier for each polygon
#' @param weights_table an optional data table of secondary weights
#'
#' @return a data.table of geoweights (area weighted raster/polygon overlap)
#'
#' @export
calc_geoweights <- function(data_source,  input_polygons, polygon_id, weights_table = NULL){

  # Create raster
  clim_raster <- raster::raster(data_source) # only reads the first band

  ## Raster cell area
  ## -----------------------------------------------

  clim_area_raster <- raster::area(clim_raster)

  ## Raster/polygon alignment
  ## -----------------------------------------------

  message(crayon::yellow('Checking for raster/polygon alignment'))

  poly_xmin <- raster::extent(input_polygons)@xmin
  poly_xmax <- raster::extent(input_polygons)@xmax
  rast_xmin <- raster::extent(clim_area_raster)@xmin
  rast_xmax <- raster::extent(clim_area_raster)@xmax

  # Shift polygons if initial longitudes don't align
  if(!dplyr::near(poly_xmax, rast_xmax, tol=1.01)) {

    message(crayon::yellow('Adjusting polygon longitude from',
                           round(poly_xmin,0), '-', round(poly_xmax,0),
                           'to', round(rast_xmin,0), '-', round(rast_xmax,0)))


    input_polygons <- sf::st_shift_longitude(input_polygons)


  }

  # Match raster and polygon crs
  crs_raster <- raster::crs(clim_area_raster)
  polygons_reproj <- sf::st_transform(input_polygons, crs = crs_raster)

  ## Raster / Polygon overlap (using data.table)
  ## -----------------------------------------------
  message(crayon::yellow('Extracting raster polygon overlap'))

  overlap <- data.table::rbindlist(exactextractr::exact_extract(clim_area_raster, polygons_reproj, progress = T, include_xy = T), idcol = "poly_id")
  overlap[, ':=' (poly_id = polygons_reproj[[polygon_id]][poly_id], cell_area_km2 = value)] # Add the unique id for each polygon based on the input col name


  ## Calculate weights
  ## -----------------------------------------------

  # Calculate area weight per grid cell
  area_weight <- overlap[, .(x, y, poly_id, w_area = coverage_fraction * cell_area_km2)] # area weight = area km2 * coverage fraction

  # IF weights = TRUE, merge secondary weights with area weights
  if(!is.null(weights_table)){

    # Data.table of secondary weights
    weights_dt <- data.table::as.data.table(weights_table)

    # Min/Max of secondary weights
    weights_xmin <- min(weights_dt$x)
    weights_xmax <- max(weights_dt$x)

    # If weights don't match raster convert them
    if(!dplyr::near(weights_xmax, rast_xmax, tol=1.01)) {

      message(crayon::yellow('Adjusting secondary weights longitude from',
                             round(weights_xmin,0), '-', round(weights_xmax,0),
                             'to', round(rast_xmin,0), '-', round(rast_xmax,0)))

      weights_dt[, x := ifelse(x < 0, x + 360, x)]

    }

    # Set key column in the merged dt table
    keycols = c("x", "y")
    data.table::setkeyv(area_weight, keycols)

    # Merge with secondary weights
    w_merged <- area_weight[weights_dt, nomatch = 0]

    # Weight in pixel = w_area * weight
    w_merged[, weight := weight * w_area]

  }

  # Normalize weights by polygon
  if(!is.null(weights_table)){

    w_norm <- w_merged[, ':=' (w_area = w_area / sum(w_area), weight = weight / sum(weight)), by = poly_id]


  } else {
    w_norm <- area_weight[, w_area := w_area / sum(w_area), by = poly_id]
  }


  message(crayon::yellow('Checking sum of weights within polygons'))
  if(!is.null(weights_table)){

    check_weights <- w_norm[, lapply(.SD, sum), by = poly_id,
                            .SDcols = c('w_area', 'weight')]

  } else{
    check_weights <- w_norm[, w_sum := sum(w_area), by=poly_id]
  }

  # Check that polygon weights sum to 1
  if (!is.null(weights_table)){
    for(i in nrow(check_weights)){

      if(!dplyr::near(check_weights$w_area[i], 1, tol=0.001) | !dplyr::near(check_weights$weight[i], 1, tol=0.001)){

        stop(crayon::red('Weights for polygon', check_weights$poly_id, 'do not sum to 1'))

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

  return(w_norm)

}

