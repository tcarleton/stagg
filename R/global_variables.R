# Without this line, check() throws a note about "undefined global functions or variables".
# This is because we are referencing column names in dplyr / data.table objects which use
# Non-Standard Evaluation.

globalVariables(
  c(
    'x',
    'y',
    'era5_grid',
    'poly_id',
    'value',
    '.',
    'coverage_fraction',
    'cell_area_km2',
    'weight',
    'w_area',
    'sum_weight',
    'w_sum',
    'is_right_xmin',
    'is_left_xmax',
    'x_low',
    'x_high',
    'y_low',
    'y_high',
    '..cols_to_keep',
    'w_area',
    'weight',
    'year',
    'month',
    'day',
    'hour',
    'minute',
    'weight_sum'

  )
)
