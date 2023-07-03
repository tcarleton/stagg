# Without this line, check() throws a note about "undefined global functions or variables".
# This is because we are referencing column names in dplyr / data.table objects which use
# Non-Standard Evaluation.

globalVariables(c("cell_area_km2", "coverage_fraction", "day", "era5_grid",
                  "hour", "month", "poly_id", "sum_weight", "value", "w_area",
                  "w_sum", "weight", "x", "y", "year", "."))
