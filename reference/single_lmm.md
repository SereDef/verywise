# Run a single linear mixed model and extract statistics

Fits a linear mixed model to a single vertex outcome using
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) and extracts
fixed effects statistics. This function is called repeatedly across all
cortical vertices during vertex-wise analysis.

## Usage

``` r
single_lmm(imp, y, y_name, model_template = NULL, weights = NULL)
```

## Arguments

- imp:

  A data.frame containing the phenotype dataset (in verywise format).

- y:

  A numeric vector of outcome values representing a single vertex from
  the super-subject matrix.

- y_name:

  String with the outcome name (as specified in the formula)

- model_template:

  Pre-compiled model object for faster estimation. `single_lmm` uses an
  "update"-based workflow instead of refitting the model from scratch.
  This minimizes repeated parsing and model construction overhead,
  significantly reducing computation time for large-scale vertex-wise
  analyses.

- weights:

  Optional string or numeric vector of weights for the linear mixed
  model. You can use this argument to specify inverse-probability
  weights. If this is a string, the function look for a column with that
  name in the phenotype data. Note that these are not normalized or
  standardized in any way. Default: `NULL` (no weights).

## Value

A list with two elements:

- `stats` - A data.frame with columns:

  - `term` - Fixed effect term names

  - `qhat` - Parameter estimates

  - `se` - Standard errors

- `resid` - A numeric vector of model residuals

- `warning` - Warning message(s) if any

## See also

[`run_vw_lmm`](https://seredef.github.io/verywise/reference/run_vw_lmm.md)
for the main interface.
