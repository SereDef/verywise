# Convert statistical result FBMs to FreeSurfer `.mgh` format

This function takes a list of
[bigstatsr::FBM](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
objects storing statistical results and writes them to
FreeSurfer-compatible `.mgh` files. It supports coefficient, standard
error, p-value, and residual maps, with optional on-the-fly
\\-\log\_{10}\\ transformation of p-values.

## Usage

``` r
convert_to_mgh(
  vw_results,
  result_path,
  fixed_terms = NULL,
  random_terms = NULL,
  stat_names = c("coef", "se", "p", "-log10p", "resid"),
  verbose = TRUE
)
```

## Arguments

- vw_results:

  A named list of
  [bigstatsr::FBM](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
  objects containing the statistical results.

- result_path:

  Character string indicating the base output path where the `.mgh`
  files will be written.

- fixed_terms:

  Vector of fixed term names to be included in output filenames (as
  "stack1", "stack2"...).

- random_terms:

  Vector of random term names used to label ICC statistics.

- stat_names:

  Character vector of statistic names to process. Default:
  `c("coef", "se","p", "-log10p","resid")`. The special name `"-log10p"`
  triggers the on-the-fly p-value transformation.

- verbose:

  Logical. Default:`TRUE`

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
