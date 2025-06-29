
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stagg

<!-- badges: start -->

<!-- badges: end -->

The R Package `stagg` aims to harmonize the preparation of
spatiotemporal raster data, in particular climate observations, for use
in statistical analyses. This allows for more efficient use of
researchers’ time, and avoids common missteps. The end goal is a greater
quantity of quality research into the study of coupled human-geophysical
systems.

Further documentation, including the working paper, AGU conference
poster, and cheat sheet are available
[here](https://www.tammacarleton.com/projects-6).

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

The `stagg` package currently does not include functions to download
climate data. There are many packages that can be used to download
climate data. For example,
[`climateR`](https://github.com/mikejohnson51/climateR) provides access
to climate datasets from over 2000 different data providers. ERA5
climate data can be download through
[`ecmwfr`](https://github.com/bluegreen-labs/ecmwfr), which provides an
interface to two European Centre for Medium-Range Weather Forecast APIs.
The [`KrigR`](https://github.com/ErikKusch/KrigR) package provides a
helpful wrapper function for the `ecmwfr` package to enable direct
downloading of ERA5 climate data in R. ECMWF also provides a general
[guide to downloading ERA5
data](https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5).

In this example we use ERA5 climate data, and the package includes a
global raster layer of ERA5 climate data which can be used as a default
in the `secondary_weights()` and `overlay_weights()` functions. However,
other climate datasets can easily be accommodated.

For compatibility with `stagg`, input climate data should be a raster or
layered raster (read in using `raster::raster()`/`terra::rast()` or
`raster::stack()`/`terra::c()`; data can originally be stored in formats
such as netCDF or .tif) with at least an annual temporal resolution,
since the most coarse temporal aggregation offered by stagg is annual.
If the temporal resolution is sub-daily, the number of layers per day
should be consistent. By default, hourly time steps are assumed (i.e.,
24 time steps per day), though a different number can be specified using
the time_interval argument in the `staggregate_*()` functions.

In order for the temporal aggregation to happen properly, stagg must
have a way of knowing the temporal information associated with the
raster stack being aggregated. This can be achieved one of two ways: 1)
the layer names use a string-date-time or string-date format; or 2) the
start date-time and temporal time steps of the raster stack must be
specified in `staggregate_*()`. For example, layer names for ERA5 and
data retrieved using `climateR` uses the format x2021.01.01.00.00.00 or
x2021.01.01 and variable_2021-01-01 respectively, both of which are
compatible with `stagg`. Alternatively, the user can specify the start
datetime and temporal interval of the data in the `staggregate_*()`
functions using the `start_date` and `time_interval` arguments.

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
  
  extent = "full"                          # The default is "full", which 
                                           # generates weights for entire 
                                           # secondary data raster. User has the 
                                           # option to define an extent for 
                                           # cropping the secondary_raster, 
                                           # which saves compute time. Input 
                                           # must be in a compatible format, such as:
                                           #    nj_extent <- terra::ext(nj_counties)
                                           #    nj_extent <- sf::st_bbox(nj_counties)
                                           #    nj_extent <- raster::extent(nj_counties)
                                           #    nj_extent <- c(-76, -73, 38, 42)
)
```

    #> Loading required namespace: raster

    #> Warning: [is.lonlat] coordinates are out of range for lon/lat

    #> Checking for raster alignment

    #> Resampling secondary_raster

    #> Creating a table of weights

``` r
#Display resulting table
cropland_weights
```

    #>           x     y    weight
    #>       <num> <num>     <num>
    #>   1: -76.25 41.75 0.2294976
    #>   2: -76.00 41.75 0.2198315
    #>   3: -75.75 41.75 0.1777267
    #>   4: -75.50 41.75 0.1362904
    #>   5: -75.25 41.75 0.1061073
    #>  ---                       
    #> 178: -74.25 38.50 0.0000000
    #> 179: -74.00 38.50 0.0000000
    #> 180: -73.75 38.50 0.0000000
    #> 181: -73.50 38.50 0.0000000
    #> 182: -73.25 38.50 0.0000000

As you can see from the output, `secondary_weights()` checks for
alignment, and rotates the `secondary_raster` coordinates if necessary.
It also resamples the data to the spatial resolution of the climate grid
using bilinear interpolation, and returns a `data.table` with latitudes,
longitudes, and cropland weights.

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

    #> Warning: [is.lonlat] coordinates are out of range for lon/lat

    #> Extracting raster polygon overlap

    #> Secondary weights fully overlap with the administrative regions.

    #> Checking sum of weights within polygons

    #> All weights sum to 1.

``` r
# Display results
print(county_weights, digits = 4)
```

    #> Key: <x, y>
    #>          x     y poly_id    w_area    weight
    #>      <num> <num>  <char>     <num>     <num>
    #>   1: -75.5 39.25     011 0.0172211 0.0360401
    #>   2: -75.5 39.25     033 0.0058164 0.0060266
    #>   3: -75.5 39.50     011 0.0170331 0.0307720
    #>   4: -75.5 39.50     033 0.3550727 0.3175946
    #>   5: -75.5 39.75     015 0.0219342 0.0193573
    #>  ---                                        
    #> 137: -74.0 40.75     013 0.0095999 0.0057512
    #> 138: -74.0 40.75     003 0.1828252 0.1307450
    #> 139: -74.0 40.75     031 0.0065823 0.0028849
    #> 140: -74.0 41.00     003 0.5410449 0.6773374
    #> 141: -74.0 41.00     031 0.0007093 0.0005442

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
                                         # data. We're converting from Kelvin to 
                                         # Celsius here.
  
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

    #> Averaging over 24 layers per day to get daily values

    #> Executing polynomial transformation

    #> Aggregating by polygon and month

``` r
# Display results
polynomial_output
```

    #>      year month poly_id  order_1  order_2  order_3
    #>     <num> <num>  <char>    <num>    <num>    <num>
    #>  1:  2024     6     011 729.4372 17894.66 442929.9
    #>  2:  2024     6     033 733.7413 18118.50 451662.2
    #>  3:  2024     6     015 729.0091 17891.97 443418.4
    #>  4:  2024     6     009 686.3299 15807.50 366563.3
    #>  5:  2024     6     007 730.9656 17987.03 446942.0
    #>  6:  2024     6     041 684.4208 15918.45 377139.6
    #>  7:  2024     6     019 705.7509 16878.38 410117.7
    #>  8:  2024     6     001 724.9572 17682.30 435310.4
    #>  9:  2024     6     005 726.9852 17802.18 440440.8
    #> 10:  2024     6     021 721.3909 17561.50 432688.5
    #> 11:  2024     6     027 687.8299 16056.75 381385.2
    #> 12:  2024     6     037 663.3331 14985.01 345520.2
    #> 13:  2024     6     023 718.5814 17424.04 427602.6
    #> 14:  2024     6     035 710.1244 17067.69 416337.1
    #> 15:  2024     6     029 720.6809 17499.87 429477.8
    #> 16:  2024     6     025 711.0168 17042.99 413123.2
    #> 17:  2024     6     039 702.9568 16690.27 401409.4
    #> 18:  2024     6     013 698.6268 16473.52 393206.4
    #> 19:  2024     6     031 670.6924 15279.35 354421.5
    #> 20:  2024     6     017 704.8866 16733.27 401281.2
    #> 21:  2024     6     003 688.9646 16026.15 377469.0
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
                                         # data. We're converting from Kelvin to 
                                         # Celsius here.
  
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

    #> Averaging over 24 layers per day to get daily values

    #> Executing spline transformation

    #> Aggregating by polygon and month

``` r
# Display output
spline_output
```

    #>      year month poly_id    value   term_1   term_2
    #>     <num> <num>  <char>    <num>    <num>    <num>
    #>  1:  2024     6     011 729.4372 303327.9 61769.47
    #>  2:  2024     6     033 733.7413 306556.0 62576.51
    #>  3:  2024     6     015 729.0091 303007.0 61689.27
    #>  4:  2024     6     009 686.3299 270997.7 53686.96
    #>  5:  2024     6     007 730.9656 304474.4 62056.11
    #>  6:  2024     6     041 684.4208 269634.4 53356.44
    #>  7:  2024     6     019 705.7509 285580.9 57335.39
    #>  8:  2024     6     001 724.9572 299968.0 60929.50
    #>  9:  2024     6     005 726.9852 301489.3 61309.88
    #> 10:  2024     6     021 721.3909 297294.4 60261.29
    #> 11:  2024     6     027 687.8299 272168.2 53986.42
    #> 12:  2024     6     037 663.3331 253919.2 49442.70
    #> 13:  2024     6     023 718.5814 295187.0 59734.39
    #> 14:  2024     6     035 710.1244 288852.4 58151.97
    #> 15:  2024     6     029 720.6809 296761.1 60127.85
    #> 16:  2024     6     025 711.0168 289513.5 58315.99
    #> 17:  2024     6     039 702.9568 283473.5 56806.75
    #> 18:  2024     6     013 698.6268 280224.6 55994.31
    #> 19:  2024     6     031 670.6924 259366.6 50793.75
    #> 20:  2024     6     017 704.8866 284915.8 57166.57
    #> 21:  2024     6     003 688.9646 272983.4 54184.86
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
                                         # data. We're converting from Kelvin to 
                                         # Celsius here.
  
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

    #> Averaging over 24 layers per day to get daily values

    #> Executing binning transformation

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
                                         # data. We're converting from Kelvin to 
                                         # Celsius here. 
  
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
    #>  1:           7033.869            3272.623
    #>  2:           6996.451            3413.340
    #>  3:           6929.858            3366.361
    #>  4:           7129.929            2141.989
    #>  5:           6935.445            3407.731
    #>  6:           6486.645            2739.454
    #>  7:           6661.205            3076.817
    #>  8:           6965.392            3233.581
    #>  9:           6902.066            3345.578
    #> 10:           6822.981            3290.399
    #> 11:           6549.321            2758.597
    #> 12:           6336.398            2383.596
    #> 13:           6834.658            3211.297
    #> 14:           6712.330            3130.655
    #> 15:           6894.821            3201.520
    #> 16:           6886.335            2978.070
    #> 17:           6794.598            2876.366
    #> 18:           6812.126            2754.917
    #> 19:           6437.996            2458.621
    #> 20:           6947.556            2769.722
    #> 21:           6755.069            2580.080
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

### A note on spatiotemporally aggregated climate data without further transformations

If desired, a user can employ `stagg` to return spatiotemporally
aggregated climate data without further transformations using the
`staggregate_polynomial()` function. The column `order_1` represents the
first order polynomial, i.e., merely spatiotemporal means (or sums) of
input climate data. The temporal aggregation is dictated by the
user-defined arguments `daily_agg` (options are “sum” “average”, and
“none”) and `time_agg` (options are “hour”, “day”, “month”, and “year”).
The spatial aggregation occurs for the user-defined administrative areas
of interests (e.g., counties), which is incorporated through the
`overlay_weights` input generated by the `overlay_weights()` function.
The user can set the argument `degree = 1` if this is the only desired
output. However, it is important to note that aggregation bias can
become severe when aggregations are taken over large spatial and/or
temporal periods before nonlinear transformations are performed. Indeed,
a main feature of the `stagg` package is that it performs the
transformations and aggregations in a sequence that avoids aggregation
bias. Thus, users are cautioned against performing transformations after
aggregation. Users who wish to use other transformations are encouraged
to submit them to the package via an issue or pull request on GitHub for
future package integration.

## Additional Resources

Below are some helpful resources on climate data and the analysis of
climate impacts on socioeconomic outcomes:

- [UCAR Climate Data Guide](https://climatedataguide.ucar.edu/)
- [Rising et al. (2020) Practical Guide to Climate
  Econometrics](https://climateestimate.net/content/getting-started.html)
- [Auffhammer et al. (2013) Using Weather Data and Climate Model Output
  in Economic Analyses of Climate
  Change](https://www.journals.uchicago.edu/doi/full/10.1093/reep/ret016)
- [Nissan et al. (2019) On the use and misuse of climate change
  projections in international
  development](https://wires.onlinelibrary.wiley.com/doi/10.1002/wcc.579)
- [Hsiang (2016) Climate
  Econometrics](https://www.annualreviews.org/content/journals/10.1146/annurev-resource-100815-095343)
- [Gortan et al. (2024) A unified dataset for pre-processed climate
  indicators weighted by gridded economic
  activity](https://www.nature.com/articles/s41597-024-03304-1)
  - [Gortan et al. (2024) Weighted Climate Dataset
    Dashboard](https://weightedclimatedata.streamlit.app/)
