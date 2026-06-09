# Refit a pre-compiled linear mixed model with a new outcome and extract statistics

Refits a linear mixed model with data from a single vertex outcome and
extracts fixed effect and model fit statistics. Errors and warnings are
caught and returned in the output list rather than signalled to the
caller, because the function is designed to be used inside parallel
loops (over vertices).

## Usage

``` r
refit_lmm(model_template_i, y)
```

## Arguments

- model_template_i:

  Pre-compiled model object for (MI) dataset `i`, generated using
  [`precompile_model()`](https://seredef.github.io/verywise/reference/precompile_model.md).
  The random-effects structure is reused as-is; only the response vector
  is replaced.

- y:

  A numeric vector of outcome values representing a single vertex from
  the super-subject matrix.

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

## Details

`verywise` uses an "update" (or rather refit)-based workflow instead of
fitting the model from scratch at each vertex. This minimizes repeated
parsing and model construction overhead, significantly reducing
computation time for large-scale vertex-wise analyses.

Variance components used for ICC and R\\^2\\ are computed directly from
[`lme4::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) with no
additional package dependencies:

- **Var(random)** — sum of diagonal elements across all random-effect
  covariance matrices.

- **Var(residual)** — residual variance \\\hat{\sigma}^2\\.

- **Var(fixed)** — variance of the marginal linear predictor
  \\\mathrm{Var}(\mathbf{X}\hat{\beta})\\.

The three R\\^2\\/ICC quantities follow Nakagawa & Schielzeth (2013):

\$\$ \mathrm{ICC} =
\frac{\sigma^2\_{\mathrm{rand}}}{\sigma^2\_{\mathrm{rand}} +
\sigma^2\_{\varepsilon}} \$\$

\$\$ R^2\_{\mathrm{marginal}} =
\frac{\sigma^2\_{\mathrm{fix}}}{\sigma^2\_{\mathrm{fix}} +
\sigma^2\_{\mathrm{rand}} + \sigma^2\_{\varepsilon}} \$\$

\$\$ R^2\_{\mathrm{conditional}} = \frac{\sigma^2\_{\mathrm{fix}} +
\sigma^2\_{\mathrm{rand}}}{\sigma^2\_{\mathrm{fix}} +
\sigma^2\_{\mathrm{rand}} + \sigma^2\_{\varepsilon}} \$\$

Division-by-zero or other numeric edge cases in ICC and R\\^2\\ are
handled by the internal helper
[`safe_calc()`](https://seredef.github.io/verywise/reference/safe_calc.md).

## References

Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for
obtaining R² from generalized linear mixed-effects models. *Methods in
Ecology and Evolution*, 4(2), 133–142.
[doi:10.1111/j.2041-210x.2012.00261.x](https://doi.org/10.1111/j.2041-210x.2012.00261.x)

## See also

[`precompile_model()`](https://seredef.github.io/verywise/reference/precompile_model.md)
and
[`run_vw_lmm()`](https://seredef.github.io/verywise/reference/run_vw_lmm.md)
for the main interface.
