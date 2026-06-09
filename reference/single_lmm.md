# Estimate a single outcome linear mixed model and extract statistics

Fits a linear mixed model to a single vertex outcome using
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) and extracts
fixed effects and model fit statistics. This function is called
repeatedly across all cortical vertices during vertex-wise analysis in
the
[`run_vw_lmm2()`](https://seredef.github.io/verywise/reference/run_vw_lmm2.md)
pipeline.

## Usage

``` r
single_lmm(
  imp,
  y,
  y_name,
  model_formula = NULL,
  REML,
  lmm_control,
  weights = NULL
)
```

## Arguments

- imp:

  A data.frame containing the phenotype dataset (in verywise format).

- y:

  A numeric vector of outcome values representing a single vertex from
  the super-subject matrix.

- y_name:

  String with the outcome name (as specified in the formula)

- model_formula:

  A model formula object. This should specify a linear mixed model
  `lme4` syntax. Example:
  `vw_thickness ~ age * sex + site + (1|participant_id)`.

- REML:

  Logical specifying whether to optimize the REML criterion (as opposed
  to the log-likelihood). Default: TRUE. Use `REML = FALSE` if you
  intend to do model comparison (using AIC output).

- lmm_control:

  Optional list (of correct class, resulting from
  [`lme4::lmerControl()`](https://rdrr.io/pkg/lme4/man/lmerControl.html)
  containing control parameters to be passed to
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) (e.g.
  optimizer choice, convergence criteria, see the `?lmerControl`
  documentation for details. Default: (uses default settings).

- weights:

  Optional string or numeric vector of weights for the linear mixed
  model. You can use this argument to specify inverse-probability
  weights. If this is a string, the function look for a column with that
  name in the phenotype data. Note that these are not normalized or
  standardized in any way. Default: `NULL` (no weights).

## Value

A named `list` with one of two shapes: **On error:**

- `error`:

  `character(1)`. The error message produced by
  [`lme4::refit()`](https://rdrr.io/pkg/lme4/man/refit.html).

**On success:**

- `stats`:

  A `data.frame` with columns `term` (character), `qhat` (numeric,
  fixed-effect estimate), and `se` (numeric, standard error), one row
  per fixed-effect term.

- `resid`:

  Named numeric vector of model residuals from
  [`stats::residuals()`](https://rdrr.io/r/stats/residuals.html), used
  downstream for smoothness estimation.

- `model_fit`:

  Unnamed numeric vector of length 5 containing, in order: `is_singular`
  (0/1), `AIC`, `ICC`, `R2_marginal`, `R2_conditional`.

- `warning`:

  `character` vector of any warning or message strings captured during
  refitting. Zero-length if none occurred. Singular-fit boundary
  warnings are suppressed from this vector when singularity is already
  flagged in `model_fit`.

## See also

[`run_vw_lmm2()`](https://seredef.github.io/verywise/reference/run_vw_lmm2.md)
for the main interface.
