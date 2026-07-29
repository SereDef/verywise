# Load and threshold coefficient maps

Loads a vertex-wise coefficient map from a specified results directory
and generates a significance mask based on a threshold method (e.g.,
Cluster-Wise Significance, FDR, uncorrected p-values, or an absolute
numeric threshold).

## Usage

``` r
load_and_mask_coef(hemi, measure, stack, res_dir, threshold)
```

## Arguments

- hemi:

  Hemisphere identifier (`"lh"` or `"rh"`).

- measure:

  The morphological measure (e.g., `"thickness"`, `"area"`).

- stack:

  The stack or contrast identifier.

- res_dir:

  Path to the directory containing the `.mgh` results files.

- threshold:

  The threshold criteria. Can be `"cws"` for cluster-wise significance,
  a string expression like `"fdr < 0.05"` or `"p < 0.01"`, a numeric
  absolute value, or `NULL` for no thresholding.

## Value

A named list containing:

- `coef`: The numeric vector of loaded coefficients.

- `mask`: A logical vector indicating significant vertices, or `NULL`.
