
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stagg

<!-- badges: start -->
<!-- badges: end -->

The R Package `stagg` aims to harmonize the preparation of
spatiotemporal data obtained from the [ERA5
project](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5)
for use in statistical analyses. This allows for more efficient useof
researchers’ time, and avoids common missteps. The end goal is a greater
quantity of quality research into the study of coupled human-geophysical
systems.

Further documentation, including the working paper, AGU conference
poster, and cheat sheet are available
[here](https://www.tammacarleton.com/projects-6)

## Installation

Although `stagg` is not yet on CRAN, you can install the development
version of `stagg` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tcarleton/stagg")
```

## Abstract

The increasing availability of high-resolution climate data has greatly
expanded the study of how the climate impacts humans and society.
However, the processing of these multi-dimensional datasets poses
significant challenges for researchers in this growing field, most of
whom are social scientists. This paper introduces `stagg` or “space-time
aggregator”, a new R package that streamlines the three critical
components of climate data processing for impact analysis: nonlinear
transformation, spatial and temporal aggregation, and spatial weighting
by additional social or economic variables. The package consolidates the
entire data processing pipeline into just a few lines of code, lowering
barriers to entry for researchers in the interdisciplinary field of
climate impacts analysis and facilitating a larger and more diverse
research community. While `stagg` is designed to be used with ERA5
reanalysis climate data, it can be easily modified for other input data.
The paper provides an overview of `stagg`’s functions and data
processing pipeline, followed by an applied example demonstrating the
package’s utility in climate impacts research. `stagg` has the potential
to be a valuable tool in generating evidence-based estimates of the
likely impacts of future climate change and quantifying the social cost
of carbon.

## Workflow

Below is example code and commentary aimed at demonstrating expected
typical usage. The order of the steps is important, as output from each
step is used in the one that follows it.

**Important**

The `stagg` package currently does not included functions to download
climate data. There are many packages that can be used to download
climate data. For example,
[`climateR`](https://github.com/mikejohnson51/climateR) provides access
to climate datasets from over 2000 different data providers. ERA5
climate data can be download through
[`ecmwfr`](https://github.com/bluegreen-labs/ecmwfr), which provides an
interface to two European Centre for Medium-Range Weather Forecast APIs.
The [`KrigR`](https://github.com/ErikKusch/KrigR) package provides a
helpful wrapper function for the `ecmwfr` package to enable direct
downloading of ERA5 climate data in R.

While the package is designed to be used with ERA5 climate data, other
types of climate data can be used. For compatibility with `stagg`
climate data should be:  
- a raster or raster stack (netCDF, .tif or other format compatible with
`raster::raster()`) - hourly or daily temporal resolution  
- if hourly, number of layers should be a multiple of 24 (i.e.,
represent whole days)  
- layer names use a character-date-time or character-date format (e.g.,
x2021.01.01.00.00.00; x2021.01.01)

``` r
library(stagg)
```

``` r
 # Using polygons outlining counties of New Jersey as administrative regions
nj_counties <- tigris::counties("NJ")
```

### Step 1 (Optional): Resample a secondary data input and generate secondary weights for Step 2

It is common when studying interactions between human and natural
systems to spatially aggregate climate data using weights derived from
another dataset of interest, such as population or cropland. This allows
the user to retrieve the climate experienced by humans or crops within a
given administrative region. To account for this, `stagg` allows for the
conversion of a raster into a data.table of weights via the
`secondary_weights()` function. These weights can then be used to
compute a weighted average of climate data over each administrative
region.

The following example shows how one would go about generating cropland
weights for the state of New Jersey.

``` r
cropland_weights <- secondary_weights(
  
  secondary_raster = cropland_nj_2015,     # A raster layer of the secondary 
                                           # variable to generate weights from
  
  grid = era5_grid,                        # A raster layer with the same 
                                           # coordinate system and spatial 
                                           # resolution as the climate data 
                                           # (defaults to the era5_grid). 
                                           # You can also pass in your climate 
                                           # data and the grid will be taken 
                                           # from its first layer
  
  extent = "full"                          # The extent to crop the  
                                           # secondary_raster to, use whenever  
                                           # possible to save time (default is 
                                           # "full"). Format must be compatible
                                           # with raster::crop()
)
```

    #> Checking for raster alignment

    #> Longitude coordinates do not match. Aligning longitudes to standard coordinates.

    #> Resampling secondary_raster

    #> Creating a table of weights

``` r
#Display resulting table
cropland_weights
```

    #>           x     y        weight
    #>       <num> <num>         <num>
    #>   1: -75.75 41.50  0.0944147354
    #>   2: -75.50 41.50  0.0796618272
    #>   3: -75.25 41.50  0.0346200160
    #>   4: -75.00 41.50 -0.0006449262
    #>   5: -74.75 41.50 -0.0137744415
    #>  ---                           
    #> 116: -74.50 38.75  0.0000000000
    #> 117: -74.25 38.75  0.0000000000
    #> 118: -74.00 38.75  0.0000000000
    #> 119: -73.75 38.75  0.0000000000
    #> 120: -73.50 38.75  0.0000000000

As you can see from the output, `secondary_weights()` checks for
alignment, and rotates the `secondary_raster` coordinates if necessary.
It also resamples the data to the spatial resolution of the grid, before
outputting a data.table with latitudes, longitudes, and cropland
weights.

### Step 2: Overlay administrative regions onto the data’s grid

A core part of `stagg`’s functionality is to aggregate gridded data to
the level of administrative regions. In order to do this, it first
calculates the portion of each region that is covered by a particular
cell.

These weights may also be scaled by the secondary weights calculated in
Step 1. This is accomplished using the `overlay_weights()` function.

``` r
county_weights <- overlay_weights(
  
  polygons = nj_counties,              # A simple features object with the 
                                       # desired polygons
  
  polygon_id_col = "COUNTYFP",         # The name of the column containing 
                                       # polygons' identifiers
  
  grid = era5_grid,                    # A raster layer with the same coordinate
                                       # system and spatial resolution as the 
                                       # climate_data (defaults to the 
                                       # era5_grid). You can also pass in your 
                                       # climate data and a grid will be taken 
                                       # from its first layer
  
  secondary_weights = cropland_weights # Optional output from Step 1, a table of
                                       # weights, determined here by cropland in
                                       # each grid cell
)
```

    #> Checking for raster/polygon alignment

    #> Aligning longitudes to standard coordinates.

    #> Extracting raster polygon overlap

    #> Checking sum of weights within polygons

    #> All weights sum to 1.

``` r
# Display results
print(county_weights, digits = 4)
```

    #>          x     y poly_id    w_area    weight
    #>      <num> <num>  <char>     <num>     <num>
    #>   1: 284.5 39.25     011 0.0172214 0.0313069
    #>   2: 284.5 39.25     033 0.0058167 0.0055363
    #>   3: 284.5 39.50     011 0.0170326 0.0332240
    #>   4: 284.5 39.50     033 0.3550772 0.3626242
    #>   5: 284.5 39.75     015 0.0219315 0.0245433
    #>  ---                                        
    #> 137: 286.0 40.75     013 0.0095999 0.0065357
    #> 138: 286.0 40.75     003 0.1828317 0.1316480
    #> 139: 286.0 40.75     031 0.0065826 0.0027441
    #> 140: 286.0 41.00     003 0.5410404 0.5931713
    #> 141: 286.0 41.00     031 0.0007093 0.0004502

You can see that the function outputs multiple rows for each polygon,
one for every grid cell it overlaps. The column `w_area` represents the
proportion of that polygon that falls within the grid cell corresponding
to the x and y coordinates. If you included secondary_weights from Step
1, as we have here, `overlay_weights()` also creates a column named
`weight`, which is determined by normalizing the secondary_weights by
`w_area`. This is what will be used in the aggregation of values. If you
wish only to use `w_area` in aggregating by polygon, you need not run
`secondary_weights()` and can omit the argument `secondary_weights` from
your call to `overlay_weights()`.

Given all of this information, we can interpret the top row in the
output as follows: About 1.7% of the area in the county represented by
COUNTYFP 011 falls within the grid cell at 284.50 degrees longitude
(0-360 range), 39.25 degrees latitude. It appears that this particular
pixel has slightly more cropland than other pixels in this polygon
though, since the cropland weight for this cell is 3.1%.

### Step 3: Transform and aggregate data using the `staggregate_*` family of functions

After completing Step 2, you are ready to transform and aggregate your
data. This is the final step before the data is ready for use in
downstream statistical analyses. The `stagg` package provides a family
of functions to perform this final step, each offering a different type
of non-linear transformation. Regardless of the specific
function,`staggregate_*()`’s workflow is to aggregate gridded values to
the daily level (this step can be skipped by specifying daily_agg =
“none”, though is recommended for faster computation), perform a
transformation on the daily values, and aggregate these values to the
administrative regions and desired temporal scale based on the
`overlay_weights()` output from Step 2.

#### Polynomial Transformation

One common transformation is to create new variables by raising a value
to an exponent. Treating these new values as independent variables
within a linear regression allows researchers to identify non-linear
trends within the data. `stagg` prepares the data for this type of
statistical analyses while preserving its high spatiotemporal resolution
by calculating the new variables from daily gridded values prior to
aggregating to the polygon and monthly/yearly level through with
`staggregate_polynomial()`.

``` r
polynomial_output <- staggregate_polynomial(
  
  data = temp_nj_jun_2024_era5 - 273.15, # A raster brick of our primary data, 
                                         # typically but not necessarily climate 
                                         # data. For now, data must start at 
                                         # midnight and be hourly. We're 
                                         # converting from Kelvin to Celsius 
                                         # here.
  
  overlay_weights = county_weights, # Output from Step 2, determined here by 
                                    # area-normalized cropland weights for grid 
                                    # cells within each county in New Jersey
  
  daily_agg = "average",            # How to aggregate hourly values to the 
                                    # daily level (options are "sum", "average",
                                    # and "none"). Here we want average daily 
                                    # temperature 
  
  time_agg = "month",               # The temporal level to aggregate daily 
                                    # transformed values to. Current options are 
                                    # "hour", day", "month", and "year". Note 
                                    # that "hour" is only available if daily_agg
                                    # is set to "none"
  
  degree = 3                        # The highest order of the polynomial. Here 
                                    # this will create variable 3 columns: 
                                    # order_1, order_2, and order_3
  )
```

    #> Averaging hourly values to get daily values

    #> Executing polynomial transformation

    #> Assuming layer name format which after removal of the first character is compatible with lubridate::as_datetime()

    #> Aggregating by polygon and month

``` r
# Display results
polynomial_output
```

    #>      year month poly_id  order_1  order_2  order_3
    #>     <num> <num>  <char>    <num>    <num>    <num>
    #>  1:  2024     6     011 727.4672 17794.91 439122.3
    #>  2:  2024     6     033 734.1446 18138.95 452446.7
    #>  3:  2024     6     015 729.0135 17892.22 443428.8
    #>  4:  2024     6     009 686.5504 15818.35 366967.2
    #>  5:  2024     6     007 730.9461 17986.23 446917.4
    #>  6:  2024     6     041 683.8770 15894.08 376303.5
    #>  7:  2024     6     019 706.1414 16895.95 410720.8
    #>  8:  2024     6     001 724.3873 17655.14 434329.6
    #>  9:  2024     6     005 726.9416 17800.11 440367.0
    #> 10:  2024     6     021 721.1470 17550.20 432291.7
    #> 11:  2024     6     027 686.6788 16005.44 379638.7
    #> 12:  2024     6     037 663.4873 14991.05 345698.7
    #> 13:  2024     6     023 718.8042 17434.85 428000.8
    #> 14:  2024     6     035 709.9632 17060.85 416117.7
    #> 15:  2024     6     029 720.6260 17498.11 429445.2
    #> 16:  2024     6     025 709.4887 16969.68 410464.1
    #> 17:  2024     6     039 702.8470 16685.24 401235.7
    #> 18:  2024     6     013 698.1990 16454.86 392588.5
    #> 19:  2024     6     031 672.3034 15344.78 356438.4
    #> 20:  2024     6     017 705.0655 16740.65 401508.9
    #> 21:  2024     6     003 688.3550 16002.56 376794.2
    #>      year month poly_id  order_1  order_2  order_3

You can see that 3 variables are created. `order_1` represents the
original values, linearly aggregated to the county, monthly level.
`order_2` and `order_3` represent the original values squared and cubed,
respectively, prior to being aggregated to the county and monthly level.
In this case, our example is only 30 days of temperature data and so
each polygon only has one row corresponding to the only month present,
June Were this a full year of data, each polygon would appear 12 times.
Note also that passing `time_agg = "day"` would create a data.table 30
times longer, with another column to the right of `month` called `day`.

#### Restricted Cubic Spline Transformation

Another type of transformation `stagg` supports is a restricted cubic
spline. This, essentially, is a piecewise function where 3rd degree
polynomials intersect at knots such that the function’s first and second
derivatives are continuous from negative infinity to positive infinity,
and the function is linear before the first knot and after the last one.
A more detailed explanation, as well as the formula used to transform
the data, can be found
[here](https://support.sas.com/resources/papers/proceedings16/5621-2016.pdf).
`staggregate_spline()` executes this formula to create K-2 new
variables, where K is the number of knots, in addition to preserving the
original untransformed value of the variable.

Computing knot locations can be very memory intensive therefore the
`stagg` package does not compute any default knot locations. For larger
data sets, users might want to load in a representative subset of data
and calculate different quantiles to help with choosing the knot
locations.

``` r
spline_output <- staggregate_spline(
  
  data = temp_nj_jun_2024_era5 - 273.15, # A raster brick of our primary data, 
                                         # typically but not necessarily climate 
                                         # data. For now, data must start at 
                                         # midnight and be hourly. We're
                                         # converting from Kelvin to Celsius 
                                         # here.
  
  overlay_weights = county_weights, # Output from Step 2, determined here by 
                                    # area-normalized cropland weights for grid
                                    # cells within each county in New Jersey.
  
  daily_agg = "average",            # How to aggregate hourly values to the 
                                    # daily level, "sum" and "average" are the 
                                    # only options. Here we want average daily 
                                    # temperature
  
  time_agg = "month",               # The temporal level to aggregate daily 
                                    # transformed values to. Current options are 
                                    # "day", "month", and "year" 
  
  knot_locs = c(0, 7.5, 12.5, 20)   # Where to place the knots

)
```

    #> Averaging hourly values to get daily values

    #> Executing spline transformation

    #> Assuming layer name format which after removal of the first character is compatible with lubridate::as_datetime()

    #> Aggregating by polygon and month

``` r
# Display output
spline_output
```

    #>      year month poly_id    value   term_1   term_2
    #>     <num> <num>  <char>    <num>    <num>    <num>
    #>  1:  2024     6     011 727.4672 301850.4 61400.11
    #>  2:  2024     6     033 734.1446 306858.5 62652.14
    #>  3:  2024     6     015 729.0135 303010.3 61690.09
    #>  4:  2024     6     009 686.5504 271163.1 53728.30
    #>  5:  2024     6     007 730.9461 304459.7 62052.45
    #>  6:  2024     6     041 683.8770 269228.9 53255.37
    #>  7:  2024     6     019 706.1414 285873.0 57408.31
    #>  8:  2024     6     001 724.3873 299540.5 60822.64
    #>  9:  2024     6     005 726.9416 301456.6 61301.70
    #> 10:  2024     6     021 721.1470 297111.6 60215.60
    #> 11:  2024     6     027 686.6788 271308.0 53771.83
    #> 12:  2024     6     037 663.4873 254033.6 49471.12
    #> 13:  2024     6     023 718.8042 295354.0 59776.14
    #> 14:  2024     6     035 709.9632 288731.9 58121.89
    #> 15:  2024     6     029 720.6260 296720.0 60117.58
    #> 16:  2024     6     025 709.4887 288367.4 58029.49
    #> 17:  2024     6     039 702.8470 283391.0 56786.11
    #> 18:  2024     6     013 698.1990 279904.4 55914.35
    #> 19:  2024     6     031 672.3034 260561.3 51090.41
    #> 20:  2024     6     017 705.0655 285049.9 57200.08
    #> 21:  2024     6     003 688.3550 272528.3 54071.38
    #>      year month poly_id    value   term_1   term_2

You can see that your output looks very similar to the table from the
polynomial transformation. The only difference here is that 4 - 2
(number of knots minus two) new variables are being created. These data
are now ready for use in a regression.

#### Binning Transformation

`stagg` can also divide the daily values into different bins specified
by the user. This can be useful in identifying outliers and
nonlinearities within the data, and accomplished by calling
`staggregate_bin()`.

``` r
bin_output <- staggregate_bin(
  
  data = temp_nj_jun_2024_era5 - 273.15, # A raster brick of our primary data, 
                                         # typically but not necessarily climate 
                                         # data. For now, data must start at 
                                         # midnight and be hourly. We're 
                                         # converting from Kelvin to Celsius
                                         # here.
  
  overlay_weights = county_weights,  # Output from Step 2, determined here by 
                                     # area-normalized cropland weights for grid 
                                     # cells within each county in New Jersey
  
  
  daily_agg = "average",             # How to aggregate hourly values to the 
                                     # daily level, "sum" and "average" are the  
                                     # only options. Here we want average daily 
                                     # temperature. 
  
  time_agg = "month",                # The temporal level to aggregate daily  
                                     # transformed values to. Current options are
                                     # "day", "month", and "year" 
  
  bin_breaks = c(0, 2.5, 5, 7.5, 10) # The values to split the data by
)
```

    #> Averaging hourly values to get daily values

    #> Executing binning transformation

    #> Assuming layer name format which after removal of the first character is compatible with lubridate::as_datetime()

    #> Aggregating by polygon and month

``` r
# Display output
bin_output
```

    #>      year month poly_id bin_ninf_to_0 bin_0_to_2.5 bin_2.5_to_5 bin_5_to_7.5
    #>     <num> <num>  <char>         <num>        <num>        <num>        <num>
    #>  1:  2024     6     011             0            0            0            0
    #>  2:  2024     6     033             0            0            0            0
    #>  3:  2024     6     015             0            0            0            0
    #>  4:  2024     6     009             0            0            0            0
    #>  5:  2024     6     007             0            0            0            0
    #>  6:  2024     6     041             0            0            0            0
    #>  7:  2024     6     019             0            0            0            0
    #>  8:  2024     6     001             0            0            0            0
    #>  9:  2024     6     005             0            0            0            0
    #> 10:  2024     6     021             0            0            0            0
    #> 11:  2024     6     027             0            0            0            0
    #> 12:  2024     6     037             0            0            0            0
    #> 13:  2024     6     023             0            0            0            0
    #> 14:  2024     6     035             0            0            0            0
    #> 15:  2024     6     029             0            0            0            0
    #> 16:  2024     6     025             0            0            0            0
    #> 17:  2024     6     039             0            0            0            0
    #> 18:  2024     6     013             0            0            0            0
    #> 19:  2024     6     031             0            0            0            0
    #> 20:  2024     6     017             0            0            0            0
    #> 21:  2024     6     003             0            0            0            0
    #>      year month poly_id bin_ninf_to_0 bin_0_to_2.5 bin_2.5_to_5 bin_5_to_7.5
    #>     bin_7.5_to_10 bin_10_to_inf
    #>             <num>         <num>
    #>  1:             0            30
    #>  2:             0            30
    #>  3:             0            30
    #>  4:             0            30
    #>  5:             0            30
    #>  6:             0            30
    #>  7:             0            30
    #>  8:             0            30
    #>  9:             0            30
    #> 10:             0            30
    #> 11:             0            30
    #> 12:             0            30
    #> 13:             0            30
    #> 14:             0            30
    #> 15:             0            30
    #> 16:             0            30
    #> 17:             0            30
    #> 18:             0            30
    #> 19:             0            30
    #> 20:             0            30
    #> 21:             0            30
    #>     bin_7.5_to_10 bin_10_to_inf

Like before, the output table features one row for every county for
every time period specified by the `time_agg` argument. What has changed
is that there is a new column for each bin created, representing the
number of days a polygon had a value that fell within that bin during
the timespan specified by the `time_agg` argument. These outputs are not
necessarily integers since the polygon is made up of pixels that are
sorted into bins and then weighted by the `overlay_weights` provided and
aggregated, here, to the county level. Here we specify bins, in degrees
Celsius, from negative infinity to 0, 0 to 2.5, 2.5 to 5, 5 to 7.5, 7.5
to 10, and 10 to infinity by passing `c(0, 2.5, 5, 7.5, 10)` to
`bin_break`. `staggregate_bin()` draws a bin between each pair of
breaks, and adds edge bins that encompass all values below the minimum
break and above the maximum break.

#### Degree Days Transformation

The final transformation offered is degree days, which measures the
degrees over a certain temperature threshold experienced (often) by
crops. This is used to generate estimates for piecewise functions.

``` r
staggregate_degree_days(
  data = temp_nj_jun_2024_era5 - 273.15, # A raster brick of our primary data, 
                                         # typically but not necessarily climate 
                                         # data. For now, data must start at 
                                         # midnight and be hourly. We're
                                         # converting from Kelvin to Celsius
                                         # here. 
  
  overlay_weights = county_weights, # Output from Step 2, determined here by 
                                    # area-normalized cropland weights for grid 
                                    # cells within each county in New Jersey
  
  # Note degree_days() does not take a daily_agg as it uses hourly values
  
  time_agg = "month",               # The temporal level to aggregate daily  
                                    # transformed values to. Current options are
                                    # "day", "month", and "year" 
  
  thresholds = c(0, 10, 20)        # Temperature thresholds between which 
                                   # separate regression coefficients can be 
                                   # estimated
)
```

    #> Skipping pre-transformation aggregation to daily level

    #> Executing degree days transformation

    #> Assuming layer name format which after removal of the first character is compatible with lubridate::as_datetime()

    #> Aggregating by polygon and month

    #>      year month poly_id threshold_ninf_to_0 threshold_0_to_10
    #>     <num> <num>  <char>               <num>             <num>
    #>  1:  2024     6     011                   0              7200
    #>  2:  2024     6     033                   0              7200
    #>  3:  2024     6     015                   0              7200
    #>  4:  2024     6     009                   0              7200
    #>  5:  2024     6     007                   0              7200
    #>  6:  2024     6     041                   0              7200
    #>  7:  2024     6     019                   0              7200
    #>  8:  2024     6     001                   0              7200
    #>  9:  2024     6     005                   0              7200
    #> 10:  2024     6     021                   0              7200
    #> 11:  2024     6     027                   0              7200
    #> 12:  2024     6     037                   0              7200
    #> 13:  2024     6     023                   0              7200
    #> 14:  2024     6     035                   0              7200
    #> 15:  2024     6     029                   0              7200
    #> 16:  2024     6     025                   0              7200
    #> 17:  2024     6     039                   0              7200
    #> 18:  2024     6     013                   0              7200
    #> 19:  2024     6     031                   0              7200
    #> 20:  2024     6     017                   0              7200
    #> 21:  2024     6     003                   0              7200
    #>      year month poly_id threshold_ninf_to_0 threshold_0_to_10
    #>     threshold_10_to_20 threshold_20_to_inf
    #>                  <num>               <num>
    #>  1:           7052.015            3207.199
    #>  2:           6998.201            3421.270
    #>  3:           6929.917            3366.407
    #>  4:           7133.817            2143.393
    #>  5:           6935.154            3407.553
    #>  6:           6482.774            2730.275
    #>  7:           6664.668            3082.725
    #>  8:           6964.307            3220.987
    #>  9:           6901.680            3344.918
    #> 10:           6821.689            3285.840
    #> 11:           6539.259            2741.031
    #> 12:           6339.128            2384.568
    #> 13:           6833.017            3218.283
    #> 14:           6709.898            3129.220
    #> 15:           6893.267            3201.758
    #> 16:           6889.351            2938.378
    #> 17:           6793.669            2874.659
    #> 18:           6806.230            2750.546
    #> 19:           6458.801            2476.480
    #> 20:           6951.579            2769.992
    #> 21:           6738.039            2582.482
    #>     threshold_10_to_20 threshold_20_to_inf

`staggregate_degree_days()` operates directly on the hourly values.
Passing a vector of length `n` to `thresholds` creates `n + 1` columns,
similar to how `staggregate_bin()` operates. For each value in the
climate raster brick (or stack), the function determines which
thresholds the value falls between. For example, a value of 15 falls
between 10 and 20. All the variables corresponding to threshold pairs
below these receive the difference between the two thresholds
(`threshold_0_to_10` gets 10) and all variables above
(`threshold_20_to_inf`) get 0. The variable in which the value falls
gets the difference between the lower threshold and the value
(`threshold_10_to20` gets 5). The low edge variable
(`threshold_ninf_to_0`) is unique in that it measures how far below the
smallest threshold a value falls. A value of -3 would get 3 for this
variable, while any value above 0 would receive 0 here. Once all values
are transformed in this way the hourly values are then aggregated to the
polygon level and temporal scale desired as with all other
`staggregate_*()` functions.
