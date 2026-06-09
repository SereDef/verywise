# Pre-compile a linear mixed model template

Fits a [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) model (on
each imputed dataset) using a dummy outcome variable. It produces a list
of fitted model templates that can be updated with vertex outcome data,
avoiding the overhead of re-specifying the random-effects structure for
every vertex/surface measure.

## Usage

``` r
precompile_model(
  formula,
  data_list,
  tmp_y,
  measure,
  weights = NULL,
  REML = TRUE,
  lmm_control = lme4::lmerControl(calc.derivs = FALSE),
  verbose
)
```

## Arguments

- formula:

  A two-sided `formula` object passed to
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html). The
  left-hand side is replaced by a dummy outcome column named
  `vw_<measure>` before fitting.

- data_list:

  A `list` of `data.frame` objects, typically constructed by
  [`imp2list()`](https://seredef.github.io/verywise/reference/imp2list.md)
  Each element must contain all variables referenced by `formula`.

- tmp_y:

  A numeric vector used as the dummy outcome. Must have the same length
  as the number of rows in each element of `data_list`.

- measure:

  A character string naming the brain surface measure (e.g.
  `"thickness"`, `"area"`). Used to construct the temporary outcome
  column `vw_<measure>` injected into each imputed dataset.

- weights:

  Either `NULL` (no observation weights), a numeric vector of weights of
  the same length as the rows in each imputed dataset, or a character
  string giving the name of a column in each dataset that contains the
  weights.

- REML:

  Logical. Whether to use restricted maximum likelihood estimation.
  Passed directly to
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html). Defaults to
  `TRUE`.

- lmm_control:

  A `lmerControl` object produced by
  [`lme4::lmerControl()`](https://rdrr.io/pkg/lme4/man/lmerControl.html),
  used to fine-tune the optimiser behaviour. Defaults to
  `lme4::lmerControl(calc.derivs = FALSE)`.

- verbose:

  Logical. If `TRUE`, emits a `cli` progress step during model
  construction.

## Value

A `list` of the same length as `data_list`, where each element is a
fitted `lmerMod` object estimated on the corresponding imputed dataset
with the dummy outcome. These objects are intended for use with
[`lme4::refit()`](https://rdrr.io/pkg/lme4/man/refit.html) or equivalent
downstream functions in the `verywise` pipeline.

## Author

Serena Defina, 2026
