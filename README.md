
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

## Installation

Although `stagg` is not yet on CRAN, you can install the development
version of `stagg` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tcarleton/stagg")
```

## Abstract

Spatiotemporal climate data is critical to interdisciplinary research on
interactions between human and climate systems. However, it must be
properly transformed and aggregated before it can be utilized in
statistical research. In particular, most statistical analyses of
coupled human-geophysical systems require careful alignment of gridded,
high-resolution climate information with irregular and often temporally
coarse administrative or survey data. While tools exist to aid in
certain aspects of climate data manipulation, there is no standardized
methodology or tools to support the full range of processing necessary
for statistical analysis. Many researchers invest time and resources
into creating individual solutions which may be closed-source,
inconsistent, or employ inadequate methods for handling nonlinearities
in spatial and temporal aggregations, while others forgo using the data
altogether. Here, we present the R package stagg, which aims to
standardize a pipeline for integrating climate and non-climate data for
statistical analysis by automating the following key facets: (1)
resampling, reprojection, and aggregation of raster data to resolve any
spatial and temporal scale mismatches between climate and non-climate
data; (2) transformations of data at grid scale for use in nonlinear
regression causal analysis; (3) aggregation of climate data over
administrative regions, with the capability to weight gridded climate
values based on non-climate datasets, such as population and cropland.
In automating these critical components, the stagg package enables
computationally efficient, simple, and generalizable preparation of
climate data for downstream use in statistical analyses, with the goal
of facilitating a greater quantity of rigorous climate research.

## Workflow

Below is example code and commentary aimed at demonstrating expected
typical usage. The order of the steps is important, as output from each
step is used in the one that follows it.

``` r
library(stagg)
```

### Step 1 (Optional): Resample a secondary data input and generate secondary weights for Step 2

It is common when studying interactions between human and natural
systems to spatially aggregate climate data using weights derived from
another raster dataset of interest, such as population or cropland. This
allows the user to retrieve the climate experienced by humans or crops
average climate experienced by humans or crops within a given
administrative region. To account for this, `stagg` allows for the
conversion of a raster into a data.table of weights via the
`secondary_weights()` function. These weights can then be used to
compute a weighted average of climate data over each administrative
region.

The following example shows how one would go about generating cropland
weights for the state of Kansas.

``` r
small_extent <- c(-95.75, -95.25, 37.25, 37.75)

cropland_weights <- secondary_weights(
  
  secondary_raster = cropland_kansas_2011, # A raster layer of the social 
                                           # variable to generate weights from
  
  grid = era5_grid,                        # A raster layer with the same 
                                           # coordinate system and spatial 
                                           # resolution as the climate data 
                                           # (defaults to the era5_grid). 
                                           # You can also pass in your climate 
                                           # data and the grid will be taken 
                                           # from its first layer
  
  extent = small_extent                   # The extent to crop the  
                                           # secondary_raster to, use whenever  
                                           # possible to save time (default is 
                                           # "full"). Format is a vector of 4 
                                           # numeric values defining boundaries
                                           # in the following order:minimum 
                                           # longitude, maximum longitude, 
                                           # minimum latitude, maximum latitude.
)
```

    #> Checking for raster alignment

    #> Adjusting raster longitude from 0 - 360 to -180 - 180

    #> Longitude ranges match

    #> Resampling secondary_raster

    #> Creating a table of weights

``` r
#Display resulting table
cropland_weights
```

    #>         x     y    weight
    #> 1: -95.75 37.50 0.3615269
    #> 2: -95.50 37.50 0.3518547
    #> 3: -95.75 37.25 0.2755257
    #> 4: -95.50 37.25 0.3613482

As you can see from the output, `secondary_weights()` checks for
alignment, and rotates the `secondary_raster` coordinates if necessary.
It also resamples the data to the spatial resolution of the grid, before
outputting a data.table with latitudes, longitudes, and cropland
weights.

Due to size constraints, cropland_kansas_2011 is very small. Here we’ll
replace the output from this example with pre-loaded, global cropland
weights but keep the same name to demonstrate how the package normally
would flow. Note that the pre-processed global cropland weights as well
as global population weights are available as part of the package and
can be used for analysis in any part of the globe.

``` r
cropland_weights <- dplyr::filter(cropland_world_2003_era5, 
                                  x >= -103, x <= -94, y >= 37, y <= 41)
```

### Step 2: Overlay administrative regions onto the data’s grid

A core part of `stagg`’s functionality is to aggregate gridded data to
the level of administrative regions. In order to do this, it first
calculates the portion of each region that is covered by a particular
cell. These weights may also be scaled by the secondary weights
calculated in Step 1. This is accomplished using the `overlay_weights()`
function.

``` r
 # Using polygons outlining counties of Kansas as administrative regions
kansas_counties <- tigris::counties("Kansas")
```

``` r
county_weights <- overlay_weights(
  
  polygons = kansas_counties,          # A simple features object with the 
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

    #> Warning in .couldBeLonLat(x, warnings = warnings): raster has a longitude/
    #> latitude crs, but coordinates do not match that

    #> Checking for raster/polygon alignment

    #> Adjusting raster longitude from 0 - 360 to -180 to 180

    #> Extracting raster polygon overlap

    #>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%

    #> Adjusting secondary weights longitude from -103 - -94 to -180 - 180

    #> Checking sum of weights within polygons

    #> All weights sum to 1.

``` r
# Display results
county_weights
```

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
output as follows: About 11% of the area in the Kansas county
represented by COUNTYFP 129 falls within the grid cell at 258 degrees
longitude (0-360 range), 37 degrees latitude. It appears that this
particular pixel has slightly less cropland than other pixels in this
polygon though, since the area-normalized cropland weight for this cell
is only around 10%.

### Step 3: Transform and aggregate data using the `staggregate_*` family of functions

After completing Step 2, you are ready to transform and aggregate your
data. This is the final step before the data is ready for use in
downstream statistical analyses. The `stagg` package provides a family
of functions to perform this final step, each offering a different type
of non-linear transformation. Regardless of the specific
function,`staggregate_*`’s workflow is to aggregate gridded values to
the daily level, perform a transformation on the daily values, and
aggregate these values to the administrative regions and desired
temporal scale based on the `overlay_weights()` output from Step 2.

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
  
  data = prcp_kansas_dec2011_era5,  # A raster brick of our primary data, 
                                    # typically but not necessarily climate 
                                    # data. For now, data must start at midnight
                                    # and be hourly.
  
  overlay_weights = county_weights, # Output from Step 2, determined here by 
                                    # area-normalized cropland weights for grid 
                                    # cells within each county in Kansas
  
  daily_agg = "sum",                # How to aggregate hourly values to the 
                                    # daily level, "sum" and "average" are the 
                                    # only options. Here we want total daily 
                                    # precipitation. 
  
  time_agg = "month",               # The temporal level to aggregate daily 
                                    # transformed values to. Current options are 
                                    # "day", "month", and "year" 
  
  degree = 3                        # The highest order of the polynomial. Here 
                                    # this will create variable 3 columns: 
                                    # order_1, order_2, and order_3
  )
```

    #> Summing hourly values to daily values

    #> Executing polynomial transformation

    #> Assuming layer name format which after removal of the first character is compatible with lubridate::as_date()

    #> Aggregating by polygon and month

``` r
# Display results
polynomial_output
```

    #>      year month poly_id    order_1      order_2      order_3
    #>   1: 2011    12     129 0.05643332 0.0007791418 1.265471e-05
    #>   2: 2011    12     187 0.05680940 0.0007600893 1.232165e-05
    #>   3: 2011    12     075 0.05149176 0.0006556775 9.875916e-06
    #>   4: 2011    12     071 0.04700092 0.0005236172 6.729562e-06
    #>   5: 2011    12     199 0.04279329 0.0004099400 4.543805e-06
    #>  ---                                                        
    #> 101: 2011    12     011 0.07018325 0.0009848000 1.687835e-05
    #> 102: 2011    12     107 0.07532783 0.0010772387 1.819141e-05
    #> 103: 2011    12     121 0.07734484 0.0011458564 2.009548e-05
    #> 104: 2011    12     091 0.07875736 0.0012305217 2.269488e-05
    #> 105: 2011    12     209 0.08044296 0.0013109832 2.494272e-05

You can see that 3 variables are created. `order_1` represents the
original values, linearly aggregated to the county, monthly level.
`order_2` and `order_3` represent the original values squared and cubed,
respectively, prior to being aggregated to the county and monthly level.
In this case, we are working with only 30 days of data to meet size
constraints, and so each polygon only has one row corresponding to the
only month present, December. Were this a full year of data, each
polygon would appear 12 times. Note also that passing `time_agg = "day"`
would create a data.table 30 times longer, with another column to the
right of `month` called `day`.

#### Restricted Cubic Spline Transformation

Another type of transformation `stagg` supports is a restricted cubic
spline. This, essentially, is a piecewise function where 3rd degree
polynomials intersect at knots such that the function’s first and second
derivatives are continuous from negative infinity to positive infinity,
and that the function is linear before the first knot and after the last
one. A more detailed explanation, as well as the formula used to
transform the data, can be found
[here](https://support.sas.com/resources/papers/proceedings16/5621-2016.pdf).
`staggregate_spline()` executes this formula to create K-2 new
variables, where K is the number of knots, in addition to preserving the
original untransformed value of the variable.

``` r
quantiles <- daily_quants(          # Use daily_quants() to calculate quantiles 
                                    # in daily values to help determine knot 
                                    # locations
  
  data = prcp_kansas_dec2011_era5,  # A raster brick of our primary data, to 
                                    # aggregate to daily level and calculate 
                                    # quantiles from
  
  overlay_weights = county_weights, # Output from Step 2, determined here by 
                                    # area-normalized cropland weights for grid
                                    # cells within each county in Kansas
  
  daily_agg = "sum",                # How to aggregate hourly values to the 
                                    # daily level, "sum" and "average" are the 
                                    # only options. Here we want total daily
                                    # precipitation
  
  probs = c(.05, .35, .65, .95)     # Probabilities from 0 to 1 to calculate 
                                    # quantiles in the daily data from
  
)
```

    #> Summing hourly values to daily values

``` r
spline_output <- staggregate_spline(
  
  data = prcp_kansas_dec2011_era5,  # A raster brick of our primary data, 
                                    # typically but not necessarily climate 
                                    # data. For now, data must start at midnight
                                    # and be hourly.
  
  overlay_weights = county_weights, # Output from Step 2, determined here by 
                                    # area-normalized cropland weights for grid
                                    # cells within each county in Kansas
  
  daily_agg = "sum",                # How to aggregate hourly values to the 
                                    # daily level, "sum" and "average" are the 
                                    # only options. Here we want total daily 
                                    # precipitation
  
  time_agg = "month",               # The temporal level to aggregate daily 
                                    # transformed values to. Current options are 
                                    # "day", "month", and "year" 
  
  knot_locs = as.vector(quantiles)  # Where to place the knots
)
```

    #> Summing hourly values to daily values

    #> Executing spline transformation

    #> Assuming layer name format which after removal of the first character is compatible with lubridate::as_date()

    #> Aggregating by polygon and month

``` r
# Display output
spline_output
```

    #>      year month poly_id      value       term_1       term_2
    #>   1: 2011    12     129 0.05643332 1.668146e-07 1.668146e-07
    #>   2: 2011    12     187 0.05680940 1.632206e-07 1.632206e-07
    #>   3: 2011    12     075 0.05149176 1.455279e-07 1.455279e-07
    #>   4: 2011    12     071 0.04700092 1.237896e-07 1.237896e-07
    #>   5: 2011    12     199 0.04279329 1.016673e-07 1.016673e-07
    #>  ---                                                        
    #> 101: 2011    12     011 0.07018325 2.058852e-07 2.058852e-07
    #> 102: 2011    12     107 0.07532783 2.265836e-07 2.265836e-07
    #> 103: 2011    12     121 0.07734484 2.358506e-07 2.358506e-07
    #> 104: 2011    12     091 0.07875736 2.461432e-07 2.461432e-07
    #> 105: 2011    12     209 0.08044296 2.578355e-07 2.578355e-07

You can see that your output looks very similar to the table from the
polynomial transformation. The only difference here is that 4 - 2
(number of knots minus two) new variables are being created. This data
is now ready for use in a regression.

#### Binning Transformation

The last tranformation `stagg` offers is to divide the daily values into
different bins specified by the user. This can be useful in identifying
outliers and nonlinearities within the data, and accomplished by calling
`staggregate_bin()`.

``` r
bin_output <- staggregate_bin(
  
  data = prcp_kansas_dec2011_era5,  # A raster brick of our primary data, 
                                    # typically but not necessarily climate 
                                    # data. For now, data must start at midnight
                                    # and be hourly.
  
  overlay_weights = county_weights, # Output from Step 2, determined here by 
                                    # area-normalized cropland weights for grid 
                                    # cells within each county in Kansas
  
  
  daily_agg = "sum",                # How to aggregate hourly values to the 
                                    # daily level, "sum" and "average" are the  
                                    # only options. Here we want total daily 
                                    # precipitation. 
  
  time_agg = "month",               # The temporal level to aggregate daily  
                                    # transformed values to. Current options are
                                    # "day", "month", and "year" 
  
  bin_breaks = c(0, 2, 4, 6)        # The values to split the data by
)
```

    #> Summing hourly values to daily values

    #> Executing binning transformation

    #> Assuming layer name format which after removal of the first character is compatible with lubridate::as_date()

    #> Aggregating by polygon and month

``` r
# Display output
bin_output
```

    #>      year month poly_id bin_ninf_to_0 bin_0_to_2 bin_2_to_4 bin_4_to_6
    #>   1: 2011    12     129      16.57524   13.42476          0          0
    #>   2: 2011    12     187      17.62299   12.37701          0          0
    #>   3: 2011    12     075      17.50892   12.49108          0          0
    #>   4: 2011    12     071      17.81776   12.18224          0          0
    #>   5: 2011    12     199      18.79249   11.20751          0          0
    #>  ---                                                                  
    #> 101: 2011    12     011      13.98931   16.01069          0          0
    #> 102: 2011    12     107      13.95431   16.04569          0          0
    #> 103: 2011    12     121      13.44070   16.55930          0          0
    #> 104: 2011    12     091      14.71583   15.28417          0          0
    #> 105: 2011    12     209      14.48745   15.51255          0          0
    #>      bin_6_to_inf
    #>   1:            0
    #>   2:            0
    #>   3:            0
    #>   4:            0
    #>   5:            0
    #>  ---             
    #> 101:            0
    #> 102:            0
    #> 103:            0
    #> 104:            0
    #> 105:            0

Like before, the output table features one row for every county for
every time period specified by the `time_agg` argument. What has changed
is that there is a new column for each bin created, representing the
number of days a polygon had a value that fell within that bin during
the timespan specified by the `time_agg` argument. These outputs are not
necessarily integers since the polygon is made up of pixels that are
sorted into bins and then weighted by the `overlay_weights` provided and
aggregated, here, to the county level. Here we specify bins from
negative infinity to 0, 0 to 2, 2 to 4, 4 to 6, and 6 to infinity by
passing c(0,2,4,6) to bin_break. `staggregate_bin` draws a bin between
each break, and adds edge bins that encompass all values below the
minimum break and above the maximum break.
