% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/staggregate.R
\name{staggregate_degree_days}
\alias{staggregate_degree_days}
\title{Degree days transformation and aggregation of climate data}
\usage{
staggregate_degree_days(data, overlay_weights, time_agg = "month", thresholds)
}
\arguments{
\item{data}{The raster brick with the data to be transformed and aggregated}

\item{overlay_weights}{A table of weights which can be generated using the
function \code{overlay_weights()}}

\item{time_agg}{the temporal scale to aggregate data to. Options are \code{'day'},
\code{'month'}, or \code{'year'}}

\item{thresholds}{A vector of temperature thresholds critical to a crop}
}
\description{
The function \code{staggregate_degree_days()} aggregates climate data to the daily
level, performs a degree days transformation on these daily values, and
aggregates the transformed values to the polygon level and desired temporal
scale
}
\examples{
degree_days_output <- staggregate_degree_days(
  data = temp_nj_jun_2024_era5 - 273.15, # Climate data to transform and aggregate
  overlay_weights = overlay_weights_nj, # Output from overlay_weights()
  time_agg = "month", # Sum the transformed daily values across months
  thresholds = c(0, 10, 20) # Calculate degree days above 0, 10, and 20
                            # degrees Celsius
  )

head(degree_days_output)

}
