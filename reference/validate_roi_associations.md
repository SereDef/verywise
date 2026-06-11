# Validate ROI association specifications

Checks that `roi_associations` is a named list of named numeric vectors
whose names correspond to known Desikan-Killiany ROI labels.

## Usage

``` r
validate_roi_associations(roi_associations, simulate_other_rois)
```

## Arguments

- roi_associations:

  Named list of named numeric vectors specifying fixed-effect beta
  coefficients per ROI. May be `NULL` or an empty
  [`list()`](https://rdrr.io/r/base/list.html).

- simulate_other_rois:

  Logical. If `FALSE`, an empty `roi_associations` is an error because
  no data would be simulated.

## Value

The validated `roi_associations` list (invisibly).
