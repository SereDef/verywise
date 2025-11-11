# Run a single linear mixed model and extract statistics

Fits a linear mixed model to a single vertex outcome using
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) and extracts
fixed effects statistics. This function is called repeatedly across all
cortical vertices during vertex-wise analysis.

## Usage

``` r
single_lmm(
  imp,
  y,
  formula,
  model_template = NULL,
  weights = NULL,
  lmm_control = lme4::lmerControl()
)
```

## Arguments

- imp:

  A data.frame containing the phenotype dataset (in verywise format).

- y:

  A numeric vector of outcome values representing a single vertex from
  the super-subject matrix.

- formula:

  An R formula object describing the linear mixed model using `lme4`
  notation.

- model_template:

  Optional pre-compiled model object for faster estimation. When
  provided, `single_lmm` will use an "update"-based workflow instead of
  refitting the model from scratch. This minimizes repeated parsing and
  model construction overhead, significantly reducing computation time
  for large-scale vertex-wise analyses. Default: `NULL`.

- weights:

  Optional string or numeric vector of weights for the linear mixed
  model. You can use this argument to specify inverse-probability
  weights. If this is a string, the function look for a column with that
  name in the phenotype data. Note that these are not normalized or
  standardized in any way. Default: `NULL` (no weights).

- lmm_control:

  Optional list (of correct class, resulting from `lmerControl()`
  containing control parameters to be passed to
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) (e.g.
  optimizer choice, convergence criteria, see the `*lmerControl`
  documentation for details. Default: (uses default settings).

## Value

A list with two elements:

- `stats` - A data.frame with columns:

  - `term` - Fixed effect term names

  - `qhat` - Parameter estimates

  - `se` - Standard errors

- `resid` - A numeric vector of model residuals

- `warning` - Warning message(s) if any

## Details

Additional parameters are currently passed to the
[`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) call using the
`lmm_control` argument.

## See also

[`run_vw_lmm`](https://seredef.github.io/verywise/reference/run_vw_lmm.md)
for the main interface.
