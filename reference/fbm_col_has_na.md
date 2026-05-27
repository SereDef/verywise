# Check if any elements are NA in FBM

Check if any row elements in a FBM are NA. This is used to clean up the
effect sizes vectors before running a meta-analysis.

## Usage

``` r
fbm_col_has_na(
  X,
  n_cores = 1,
  row.ind = bigstatsr::rows_along(X),
  col.ind = bigstatsr::cols_along(X),
  row.mask = NULL,
  col.mask = NULL,
  verbose = FALSE
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

A (large) logical vector for where a column contained at least 2
observed values.

## Author

Serena Defina, 2026.
