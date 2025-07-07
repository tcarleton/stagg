# Function to generate SpatRaster stack for testing
gen_test_rast <- function(){
  terra::rast(

  )
}


# Function to generate polygons
gen_two_triangles <- function(edge_case = 'none'){

  # Change latitudes based on which edge case we're looking for
  if(edge_case == 'span_dl'){

    # Relocate the polygons so they span the international date_line (-180/180 in standard coords, 180 in standard coords)
    shift = 179
  } else if(edge_case == 'span_pm'){

    # Relocate the polygons so they span the prime meridian (0 in standard coords, 0/360 in climate coords)
    shift = -1
  }else{

    # Shift by 10 to pull it off the 0 line
    shift = 10
  }


  # Make two triangles of width and height 2, each with a point on the origin
  # and together making a square (diagonal like this: [b/a])
  output <- data.frame(
    poly = c(rep('a', 3), rep('b', 3)),
    lon = c(0, 2, 2, 0, 0, 2),
    lat = c(0, 0, 2, 0, 2, 2)
  ) |>

    # Shift by desired amount
    dplyr::mutate(lon = lon + shift) |>

    # Convert coords to points
    sf::st_as_sf(coords = c('lon', 'lat')) |>

    # summarize points to polygons
    dplyr::group_by(poly) |>
    dplyr::summarize(geometry = sf::st_combine(geometry)) |>
    sf::st_cast('POLYGON') |>
    sf::st_set_crs(4326)

  if(edge_case == 'span_dl'){
    output <- output |>
      sf::st_wrap_dateline()
  }

  return(output)

}



# ==============================================================================
# Functions not run in tests, but useful for interactive maintenance
# ==============================================================================

# Function to plot overlays (can add other ggplot2 features w/ "+ geom...")
plot_overlay <- function(raster, polygons, layer = 1, var_name = 'value'){

  # Coerce raster to SpatRaster, get single layer, and rename
  raster <- validate_data(raster)[[layer]]
  names(raster) <- 'value'

  # Coerce to sf
  polygons <- sf::st_as_sf(polygons)

  # Plot
  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = as.data.frame(raster, xy = TRUE),
      mapping = ggplot2::aes(x, y, fill = value)
    ) +
    ggplot2::geom_sf(
      data = polygons,
      fill = NA,
      color = 'red'
    ) +
    ggplot2::labs(
      fill = var_name
    ) +
    ggplot2::theme(
      legend.position = 'top',
      panel.background = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_viridis_c()

}
