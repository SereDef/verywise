# Convert statistical result FBMs to FreeSurfer `.mgh` format

This function takes a list of
[FBM](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
objects storing statistical results and writes them to
FreeSurfer-compatible `.mgh` files. It supports coefficient, standard
error, p-value, and residual maps, with optional on-the-fly
\\-\log\_{10}\\ transformation of p-values.

## Usage

``` r
convert_to_mgh(
  vw_results,
  result_path,
  stacks,
  stat_names = c("coef", "se", "p", "-log10p", "resid"),
  verbose = TRUE
)
```

## Arguments

- vw_results:

  A named list of
  [FBM](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
  objects containing the statistical results.

- result_path:

  Character string indicating the base output path where the `.mgh`
  files will be written.

- stacks:

  Vector of stack identifiers (e.g. hemisphere stack IDs) to be included
  in output filenames (for non-residual stats).

- stat_names:

  Character vector of statistic names to process. Default:
  `c("coef","se","p","-log10p","resid")`. The special name `"-log10p"`
  triggers the on-the-fly p-value transformation.

- verbose:

  Logical. Default: `TRUE`

## Value

Invisibly returns `NULL`. Side effects: `.mgh` files are written to
disk.

## Details

For residuals, all rows are written into a single `.mgh` file. For other
stats, data is written to one file per row (i.e. term). Note, by
default, the p vector is cast down to float (single precision, 32‑bit),
meaning about 7 significant decimal digits are stored (accurately). This
is done to cut disk space needed for the results in half.

When computing \\-\log\_{10}(p)\\ values, the transformation is applied
**in chunks of columns** to avoid loading the full FBM into memory.

## Examples

``` r
if (FALSE) { # \dontrun{
library(bigstatsr)

# Dummy vw_results list with small FBMs
vw_results <- list(
  coef = FBM(5, 10, init = rnorm(50)),
  se   = FBM(5, 10, init = runif(50, 0.1, 1)),
  p    = FBM(5, 10, init = runif(50, 0, 1)),
  resid= FBM(20, 10, init = rnorm(200))
)

convert_to_mgh(
  vw_results = vw_results,
  result_path = "my_results/stat",
  stacks = 1:5
)
} # }
```
