# Check if all row elements are 0 in FBM

Check if all row elements in a FBM are 0. This is used to clean up the
super-subject matrix in
[`mask_cortex`](https://seredef.github.io/verywise/reference/mask_cortex.md).

## Usage

``` r
fbm_col_has_0(
  X,
  n_cores = 1,
  row.ind = bigstatsr::rows_along(X),
  col.ind = bigstatsr::cols_along(X),
  row.mask = NULL,
  col.mask = NULL,
  verbose = TRUE
)
```

## Arguments

- X:

  : the file-backed matrix (FBM) object

- n_cores:

  : (default = 1) number of cores for parellalization

- row.ind:

  : indicator for rows

- col.ind:

  : indicator for columns

- row.mask:

  : (default = NULL) specify a subset of rows

- col.mask:

  : (default = NULL) specify a subset of columns

- verbose:

  : (default = TRUE)

## Value

A (large) logical vector for where rows are all 0.

## Author

Serena Defina, 2024.
