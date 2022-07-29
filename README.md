
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stagg

<!-- badges: start -->
<!-- badges: end -->

The R Package `stagg` aims to harmonize the preparation spatiotemporal
data for use in statistical analyses. This allows for more efficient use
of researchers’ time, and avoids common missteps. The end goal is a
greater quantity of quality research into the study of coupled
human-geophysical systems.

## Abstract (Work in Progress)

Spatiotemporal climate data must be properly transformed and aggregated
before it can be utilized in statistical research. In particular, most
statistical analyses of coupled human-geophysical systems require a
careful aligning of gridded, high-resolution climate information with
irregular and often temporally coarse administrative or survey data. For
researchers familiar with popular climate datasets and the steps
required to use them, the preparation process is tedious; to the
uninitiated, it can act as a barrier to entry into the study of
interactions between human and natural systems. Here, we present the R
package stagg, which aims to reduce the time and expertise required to
use climate data for statistical analysis by automating the following
key facets: (1) resampling, reprojection, and aggregation of raster data
to resolve any spatial and temporal scale mismatches; (2) nonlinear
transformations of data at grid scale for use in nonlinear regression
causal analysis; (3) aggregation of climate data over administrative
regions, with the capability to weight gridded climate values based on
non-climate datasets, such as population and cropland.

In automating these critical components, the stagg package enables
computationally efficient, simple, and generalizable preparation of
climate data for downstream use in statistical analyses.

## Workflow

Below is example code and commentary aimed at demonstrating expected
typical usage. The order of the steps is important, as output from each
step is used in the one that follows it.

``` r
library(stagg)
```

### Step 1 (Optional): Resample a secondary data input and generate secondary weights for Step 2

It is common when studying interactions between human and natural
systems to weight a climate variable by another variable of interest
such as population or cropland. This allows the user to retrieve the
climate experienced by humans or crops. Thus, `stagg` allows for the
conversion of a raster into a data.table of weights via the
`secondary_weights()` function.

The following example shows how one would go about generating cropland
weights for the state of kansas.

``` r
kansas_extent <- c(-95.75, -95.25, 37.25, 37.75)

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
  
  extent = kansas_extent # The extent to crop the  
                                           # secondary_raster to, use whenever  
                                           # possible to save time (default is 
                                           # "full")
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

### Step 2: Overlay administrative regions onto the your data’s grid

A core part of `stagg`’s functionality is to aggregate gridded data to
the level of administrative regions. In order to do this, it first
calculates the portion of each region is covered by a particular cell,
or the weight each grid cell has for a given administrative region.
These weights may also be scaled by the secondary weights calcualted in
Step 1. This is accomplished using the `overlay_weights()` function.

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

    #> Adjusting polygon longitude from -102 - -95 to 0 - 360

    #> Extracting raster polygon overlap

    #>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%

    #> Adjusting secondary weights longitude from -103 - -94 to 0 - 360

    #> Checking sum of weights within polygons

    #> All weights sum to 1.

``` r
# Display results
county_weights
```

    #>          x     y poly_id     w_area      weight
    #>   1: 258.0 37.00     129 0.11483006 0.102767664
    #>   2: 258.0 37.25     129 0.21733466 0.197589254
    #>   3: 258.0 37.50     187 0.21911435 0.186872708
    #>   4: 258.0 37.50     129 0.01223492 0.013395133
    #>   5: 258.0 37.75     075 0.08685436 0.081127788
    #>  ---                                           
    #> 854: 265.5 38.75     091 0.01778473 0.017353580
    #> 855: 265.5 38.75     121 0.01046101 0.008027165
    #> 856: 265.5 39.00     091 0.02229991 0.012353309
    #> 857: 265.5 39.00     209 0.03417110 0.020933949
    #> 858: 265.5 39.25     209 0.02636328 0.022149179

You can see that the function outputs multiple rows for each polygon,
one for every grid cell it overlaps. The column `w_area` represents the
proportion of that polygon that falls within the grid cell corresponding
to the x and y coordinates. If you included secondary_weights from Step
1, as we have here, `overlay_weights()` also creates a column named
`weight`, which is determined by normalizing the secondary_weights by
`w_area`. This is what will be used in the aggregation of values. If you
wish only to use `w_area` in aggregating by polygon, you need not run
`secondary_weights()` and can omit the argument `secondary_weights` from
your call to `overlay_weights()`. Once again, note that the function is

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
  
  prcp_kansas_dec2011_era5, # A raster brick of our primary data, typically but 
                            # not necessarily climate data. For now, data must 
                            # start at midnight and be hourly.
  
  county_weights,           # Output from Step 2, determined here by 
                            # area-normalized cropland weights for grid cells 
                            # within each county in Kansas
  
  daily_agg = "sum",        # How to aggregate hourly values to the daily level,
                            # "sum" and "average" are the only options. Here we 
                            # want total daily precipitation. 
  
  time_agg = "month",       # The temporal level to aggregate daily transformed 
                            # values to. Current options are "day", "month", and
                            # "year" 
  
  degree = 3                # The highest order of the polynomial. Here this 
                            # will create variable 3 columns: x, x^2, and x^3
  )
```

    #> Summing hourly values to daily values

    #> Executing polynomial transformation

    #> Aggregating by polygon and month

``` r
# Display results
polynomial_output
```

    #>           year    month poly_id    order_1      order_2      order_3
    #>   1: year_2011 month_12     129 0.05643332 0.0007791418 1.265471e-05
    #>   2: year_2011 month_12     187 0.05680940 0.0007600893 1.232165e-05
    #>   3: year_2011 month_12     075 0.05149176 0.0006556775 9.875916e-06
    #>   4: year_2011 month_12     071 0.04700092 0.0005236172 6.729562e-06
    #>   5: year_2011 month_12     199 0.04279329 0.0004099400 4.543805e-06
    #>  ---                                                                
    #> 101: year_2011 month_12     011 0.07018325 0.0009848000 1.687835e-05
    #> 102: year_2011 month_12     107 0.07532783 0.0010772387 1.819141e-05
    #> 103: year_2011 month_12     121 0.07734484 0.0011458564 2.009548e-05
    #> 104: year_2011 month_12     091 0.07875736 0.0012305217 2.269488e-05
    #> 105: year_2011 month_12     209 0.08044296 0.0013109832 2.494272e-05

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

#### Restricted Cubic Spline Transfromation

Another type of transformation `stagg` supports is a restricted cubic
spline. This, essentially, is a piecewise function where 3rd degree
polynomials intersect at knots such that the function’s first and second
derivatives are continuous from negative infinity to positive infinity,
and that the function is linear before the first knot and after the last
one. A more detailed explanation, as well as the formula used to
transform the data, can be found
[here](https://support.sas.com/resources/papers/proceedings16/5621-2016.pdf).
`staggregate_spline()` executes this formula to create number of knots
minus two new variables, in addition to preserving the original
untransformed values.

``` r
spline_output <- staggregate_spline(
  
  prcp_kansas_dec2011_era5, # A raster brick of our primary data, typically but 
                            # not necessarily climate data. For now, data must 
                            # start at midnight and be hourly.
  
  county_weights,           # Output from Step 2, determined here by 
                            # area-normalized cropland weights for grid cells 
                            # within each county in Kansas
  
  daily_agg = "sum",        # How to aggregate hourly values to the daily level,
                            # "sum" and "average" are the only options. Here we 
                            # want total daily precipitation. 
  
  time_agg = "month",       # The temporal level to aggregate daily transformed 
                            # values to. Current options are "day", "month", and
                            # "year" 
  
  knot_locs = c(1, 2, 3, 4) # Where to place the knots. Most likely, you'd want 
                            # to calculate different quantiles in your daily 
                            # gridded values and choose knot_locs based on 
                            # those.
)
```

    #> Summing hourly values to daily values

    #> Executing spline transformation

    #> Aggregating by polygon and month

``` r
# Display output
spline_output
```

    #>           year    month poly_id untransformed_value term_1 term_2
    #>   1: year_2011 month_12     129          0.05643332      0      0
    #>   2: year_2011 month_12     187          0.05680940      0      0
    #>   3: year_2011 month_12     075          0.05149176      0      0
    #>   4: year_2011 month_12     071          0.04700092      0      0
    #>   5: year_2011 month_12     199          0.04279329      0      0
    #>  ---                                                             
    #> 101: year_2011 month_12     011          0.07018325      0      0
    #> 102: year_2011 month_12     107          0.07532783      0      0
    #> 103: year_2011 month_12     121          0.07734484      0      0
    #> 104: year_2011 month_12     091          0.07875736      0      0
    #> 105: year_2011 month_12     209          0.08044296      0      0

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
  
  prcp_kansas_dec2011_era5, # A raster brick of our primary data, typically but 
                            # not necessarily climate data. For now, data must 
                            # start at midnight and be hourly.
  
  county_weights,           # Output from Step 2, determined here by 
                            # area-normalized cropland weights for grid cells 
                            # within each county in Kansas
  
  
  daily_agg = "sum",        # How to aggregate hourly values to the daily level,
                            # "sum" and "average" are the only options. Here we 
                            # want total daily precipitation. 
  
  time_agg = "month",       # The temporal level to aggregate daily transformed 
                            # values to. Current options are "day", "month", and
                            # "year" 
  
  binwidth = 2,             # The width of the bins to draw (overrides any 
                            # num_bin argument)
  
  min = 0,                  # The smallest value that non-edge bins must cover,
                            # use this if you are chunking your data
  
  max = 6,                  # The largest value that non-edge bins must cover,
                            # use this if you are chunking your data
  
  center_on = 4,            # Where to center the first bin drawn on. Bins of 
                            # equal width will be drawn around this one until 
                            # both min and max are contained by non-edge bins
)
```

    #> Summing hourly values to daily values

    #> Binwidth argument supplied will override num_bins

    #> Non-edge bins extend beyond max value

    #> Non-edge bins extend beyond min value

    #> Executing binning transformation

    #> Aggregating by polygon and month

``` r
# Display output
bin_output
```

    #>           year    month poly_id neg_inf_to_-1 -1_to_1 1_to_3 3_to_5 5_to_7
    #>   1: year_2011 month_12     129             0      30      0      0      0
    #>   2: year_2011 month_12     187             0      30      0      0      0
    #>   3: year_2011 month_12     075             0      30      0      0      0
    #>   4: year_2011 month_12     071             0      30      0      0      0
    #>   5: year_2011 month_12     199             0      30      0      0      0
    #>  ---                                                                      
    #> 101: year_2011 month_12     011             0      30      0      0      0
    #> 102: year_2011 month_12     107             0      30      0      0      0
    #> 103: year_2011 month_12     121             0      30      0      0      0
    #> 104: year_2011 month_12     091             0      30      0      0      0
    #> 105: year_2011 month_12     209             0      30      0      0      0
    #>      7_to_inf
    #>   1:        0
    #>   2:        0
    #>   3:        0
    #>   4:        0
    #>   5:        0
    #>  ---         
    #> 101:        0
    #> 102:        0
    #> 103:        0
    #> 104:        0
    #> 105:        0

Like before, the output table features one row for every county for
every time period specified by the `time_agg` argument. What has changed
is that there is a new column for each bin created, representing the
number of days a polygon had a value that fell within that bin during
the timespan specified by the `time_agg` argument. These outputs are not
necessarily integers since the polygon is made up of pixels that are
sorted into bins and then weighted by the `overlay_weights` provided and
aggregated, here, to the county level.

Because there are many ways to specify bins, `staggregate_bin()` has
many optional arguments which can influence the placement, extent, and
width of the bins. First, the non-edge bins are all of equal width. This
is taken from the `binwidth` argument if it is supplied, or otherwise
from the `num_bins` argument by dividing the range (`max` minus `min`)
by the number of bins. The max and min, if not supplied are taken
directly from the maximum and minimum values in the data. Once the width
has been established, a bin is drawn using one of the placement
arguments: `start_on`, `center_on`, or `end_on`, which draw the bin’s
left edge, center, or right edge on that value, respectively. If this
value falls outside the range, it is moved over by a bin-width at a time
until it is within the range. If no placement value is specified, `min`
is passed to `start_on`. Bins are then constructed around that bin until
the full range is covered. Note that if you specify `num_bins` but
choose a placement where the max or min value will be overlapped, you
will get one more non-edge bin than requested. Lastly, edge bins, from
negative infinity to the start of the leftmost bin and from the end of
the rightmost bin to infinity, are constructed to capture any other
data.

## Installation

Although `stagg` is not yet on CRAN, you can install the development
version of `stagg` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tcarleton/stagg")
```
