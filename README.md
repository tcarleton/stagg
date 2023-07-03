
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
  
  secondary_raster = cropland_kansas_2011, # A raster layer of the secondary 
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
cropland_weights <- dplyr::filter(cropland_world_2015_era5, 
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

    #> Warning in .couldBeLonLat(x, warnings = warnings): raster has a
    #> longitude/latitude crs, but coordinates do not match that

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
  
  data = temp_kansas_jan_2020_era5, # A raster brick of our primary data, 
                                    # typically but not necessarily climate 
                                    # data. For now, data must start at midnight
                                    # and be hourly.
  
  overlay_weights = county_weights, # Output from Step 2, determined here by 
                                    # area-normalized cropland weights for grid 
                                    # cells within each county in Kansas
  
  daily_agg = "average",            # How to aggregate hourly values to the 
                                    # daily level (options are "sum", "average",
                                    # and "none"). Here we want average daily 
                                    # temperature 
  
  time_agg = "month",               # The temporal level to aggregate daily 
                                    # transformed values to. Current options are 
                                    # "hour", day", "month", and "year". Note 
                                    # that "hour" is only availabel if daily_agg
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

    #>      year month poly_id   order_1  order_2   order_3
    #>   1: 2020     1     129 112.92002 640.0153 3412.8025
    #>   2: 2020     1     187 102.00318 590.2717 3012.2352
    #>   3: 2020     1     075  80.82003 474.7735 1970.4388
    #>   4: 2020     1     071  69.47838 418.2560 1316.2572
    #>   5: 2020     1     199  66.90692 413.5508 1190.2729
    #>  ---                                                
    #> 101: 2020     1     011  79.68634 947.9243 7415.2695
    #> 102: 2020     1     107  63.81250 878.3050 5521.9911
    #> 103: 2020     1     121  43.60497 800.5253 3118.5685
    #> 104: 2020     1     091  29.09490 764.9933 1514.0289
    #> 105: 2020     1     209  15.17167 756.0888  207.5949

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
spline_output <- staggregate_spline(
  
  data = temp_kansas_jan_2020_era5, # A raster brick of our primary data, 
                                    # typically but not necessarily climate 
                                    # data. For now, data must start at midnight
                                    # and be hourly.
  
  overlay_weights = county_weights, # Output from Step 2, determined here by 
                                    # area-normalized cropland weights for grid
                                    # cells within each county in Kansas
  
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

    #> year month poly_id     value   term_1      term_2
    #>   1: 2020     1     129 112.86916 3510.695   0.9194369
    #>   2: 2020     1     187 101.84892 3161.328   1.8694796
    #>   3: 2020     1     075  80.78173 2236.894   0.3241656
    #>   4: 2020     1     071  69.47118 1700.457   0.0000000
    #>   5: 2020     1     199  66.88787 1631.257   0.0000000
    #>  ---                                                  
    #> 101: 2020     1     011  79.81634 8122.806 413.2174304
    #> 102: 2020     1     107  63.92166 6609.667 268.2332054
    #> 103: 2020     1     121  43.53806 4795.768 148.4238966
    #> 104: 2020     1     091  28.85431 3780.715 122.6036406
    #> 105: 2020     1     209  15.57286 3171.073 107.3190563

You can see that your output looks very similar to the table from the
polynomial transformation. The only difference here is that 4 - 2
(number of knots minus two) new variables are being created. This data
is now ready for use in a regression.

#### Binning Transformation

`stagg` can also divide the daily values into different bins specified
by the user. This can be useful in identifying outliers and
nonlinearities within the data, and accomplished by calling
`staggregate_bin()`.

``` r
bin_output <- staggregate_bin(
  
  data = temp_kansas_jan_2020_era5,  # A raster brick of our primary data, 
                                     # typically but not necessarily climate 
                                     # data. For now, data must start at midnight
                                     # and be hourly.
  
  overlay_weights = county_weights,  # Output from Step 2, determined here by 
                                     # area-normalized cropland weights for grid 
                                     # cells within each county in Kansas
  
  
  daily_agg = "average",             # How to aggregate hourly values to the 
                                     # daily level, "sum" and "average" are the  
                                     # only options. Here we want average daily 
                                     # temperature. 
  
  time_agg = "month",                # The temporal level to aggregate daily  
                                     # transformed values to. Current options are
                                     # "day", "month", and "year" 
  
  bin_breaks = c(0, 2.5, 5, 7.5, 10) # The values to split the data by

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
    #>   1: 2020     1     129      2.549250     5.636038    13.299869     7.514844
    #>   2: 2020     1     187      3.113396     6.702054    12.197221     7.079204
    #>   3: 2020     1     075      4.497534     8.875200    11.645565     5.107956
    #>   4: 2020     1     071      5.102531    10.185469    10.569460     5.142540
    #>   5: 2020     1     199      5.000000    10.673271    10.326729     5.000000
    #>  ---                                                                        
    #> 101: 2020     1     011      7.437066     8.605006     6.957929     3.000000
    #> 102: 2020     1     107      8.069160    10.482130     4.537403     2.964308
    #> 103: 2020     1     121      9.421720    11.553613     4.032575     2.100239
    #> 104: 2020     1     091      9.281905    11.718095     4.018300     4.075213
    #> 105: 2020     1     209     10.279019    10.720981     4.566454     4.414998
    #>      bin_7.5_to_10 bin_10_to_inf
    #>   1:    2.00000000      0.000000
    #>   2:    1.90812541      0.000000
    #>   3:    0.87374504      0.000000
    #>   4:    0.00000000      0.000000
    #>   5:    0.00000000      0.000000
    #>  ---                            
    #> 101:    2.00000000      3.000000
    #> 102:    2.74871671      2.198282
    #> 103:    2.48594898      1.405904
    #> 104:    0.90648694      1.000000
    #> 105:    0.01854763      1.000000

Like before, the output table features one row for every county for
every time period specified by the `time_agg` argument. What has changed
is that there is a new column for each bin created, representing the
number of days a polygon had a value that fell within that bin during
the timespan specified by the `time_agg` argument. These outputs are not
necessarily integers since the polygon is made up of pixels that are
sorted into bins and then weighted by the `overlay_weights` provided and
aggregated, here, to the county level. Here we specify bins, in degrees
Celsius, from negative infinity to 0, 0 to 2.5, 2.5 to 5, 5 to 7.5, 7.5
to 10, and 10 to infinity by passing c(0, 2.5, 5, 7.5, 10) to bin_break.
`staggregate_bin` draws a bin between each pair of breaks, and adds edge
bins that encompass all values below the minimum break and above the
maximum break.

#### Degree Days Transformation

The final transformation offered is degree days, which measures the
degrees over a certain temperature threshold experienced (often) by
crops. This is used to generate estimates for piecewise functions.

``` r
staggregate_degree_days(
  data = temp_kansas_jan_2020_era5, # A raster brick of our primary data, 
                                    # typically but not necessarily climate 
                                    # data. For now, data must start at midnight
                                    # and be hourly.
  
  overlay_weights = county_weights, # Output from Step 2, determined here by 
                                    # area-normalized cropland weights for grid 
                                    # cells within each county in Kansas
  
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
    #>   1: 2020     1     129            539.9509          2994.173
    #>   2: 2020     1     187            656.5378          2871.995
    #>   3: 2020     1     075            820.5155          2598.197
    #>   4: 2020     1     071            874.7042          2422.126
    #>   5: 2020     1     199            890.6823          2370.442
    #>  ---                                                         
    #> 101: 2020     1     011            793.1640          2437.138
    #> 102: 2020     1     107            925.6238          2250.300
    #> 103: 2020     1     121           1100.3957          2017.507
    #> 104: 2020     1     091           1239.0013          1841.484
    #> 105: 2020     1     209           1388.8817          1672.952
    #>      threshold_10_to_20 threshold_20_to_inf
    #>   1:          255.85814                   0
    #>   2:          232.61921                   0
    #>   3:          161.99939                   0
    #>   4:          120.05884                   0
    #>   5:          126.00658                   0
    #>  ---                                       
    #> 101:          268.49831                   0
    #> 102:          206.82323                   0
    #> 103:          129.40779                   0
    #> 104:           95.79451                   0
    #> 105:           80.04994                   0

`staggregate_degree_days()` operates directly on the hourly values.
Passing a vector of length n to `thresholds` creates n + 1 columns,
similar to how `staggregate_bin()` opperates. For each value in the
climate raster brick (or stack), the function determines which
thresholds the value falls between. For example, a value of 15 falls
between 10 and 20. All the variables corresponding to threshold pairs
below these receive the difference between the two thresholds
(threshold_0\_to_10 gets 10) and all variables above
(threshold_20_to_inf) get 0. The variable in which the value falls gets
the difference between the lower threshold and the value
(threshold_10_to20 gets 5). The low edge variable (threshold_ninf_to_0)
is unique in that it measures how far below the smallest threshold a
value falls. A value of -3 would get 3 for this variable, while any
value above 0 would receive 0 here. Once all values are transformed in
this way the hourly values are then aggregated to the polygon level and
temporal scale desired as with all other `staggregate_*` functions.
