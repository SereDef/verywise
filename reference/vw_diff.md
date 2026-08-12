# Compute vertex-wise differences between two maps

Calculates the absolute difference (`a - b`) between two cortical
surface maps for both the left and right hemispheres. Prints statistical
summaries of the difference distribution.

## Usage

``` r
vw_diff(
  lh_a = NULL,
  lh_b = NULL,
  rh_a = NULL,
  rh_b = NULL,
  label_a = "a",
  label_b = "b",
  cutoffs = 0,
  digits = 3
)

vw_diff(
  lh_a = NULL,
  lh_b = NULL,
  rh_a = NULL,
  rh_b = NULL,
  label_a = "a",
  label_b = "b",
  cutoffs = 0,
  digits = 3
)
```

## Arguments

- lh_a, lh_b:

  Left hemisphere data maps (numeric vectors or file paths) to contrast.

- rh_a, rh_b:

  Right hemisphere data maps (numeric vectors or file paths) to
  contrast.

- label_a:

  Character string representing the name of the first map (`a`).

- label_b:

  Character string representing the name of the second map (`b`).

- cutoffs:

  Numeric vector of absolute-difference thresholds to report, in
  addition to the base `a < b` comparison at 0. For each nonzero cutoff
  `c`, prints the count/percentage of vertices where `a` exceeds `b` by
  more than `c`, and where `b` exceeds `a` by more than `c`. `0` is
  always included even if not explicitly passed. Default: `0`.

- digits:

  Integer, number of decimal places for the summary statistics row.
  Default: `3`.

## Value

A named list containing `lh` and `rh` numeric vectors representing the
calculated differences (`a - b`).
