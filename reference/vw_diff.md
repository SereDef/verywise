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
  label_b = "b"
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

## Value

A named list containing `lh` and `rh` numeric vectors representing the
calculated differences (`a - b`).
