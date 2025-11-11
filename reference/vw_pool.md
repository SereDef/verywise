# Pool `lme4::lmer` model output across imputed datasets

This function combines estimates, standard errors and p-values across
imputed datasets, for a single model (i.g., one vertex). The function
was largely taken from
[`mice::pool()`](https://amices.org/mice/reference/pool.html) and
[`mice::summary.mipo()`](https://amices.org/mice/reference/mipo.html)
code. It averages the estimates of the complete data model, and computes
relevant statistics, following Rubin's rules (Rubin, 1987, p. 76).

**P-values estimation** Default: P-values are estimated using the
*t-as-z* approach. This is known to the anti-conservative for small
sample sizes but provides a computationally efficient solution. Type I
error control is addressed more rigorously at the cluster-forming stage.
Wald Chi-square test Satterthwaite Approximation

The residuals of the model are currently simply averaged across imputed
datasets, for lack of a better idea of how to combine them.

## Usage

``` r
vw_pool(out_stats, m, pvalue_method = "t-as-z", min_pvalue = 2^-149)
```

## Arguments

- out_stats:

  Output of
  [`single_lmm`](https://seredef.github.io/verywise/reference/single_lmm.md),
  i.e.: a list containing:

  - `"stats"`: a dataframe with estimates and SEs for each fixed effect
    term)

  - `"resid"`: a vector of residuals for the given model.

  - `"warning"`: a character vector with warning messages (if any)

  Note that if a model failed for that dataset `out_stats` has form
  `list("error"="Error message")`

- m:

  Integer indicating the number of (imputed) dataset

- pvalue_method:

  String indicating which approximation methods should be used to
  compute p-values. Options: `"t-as-z"`, `"wald-chi2"`,
  `"satterthwaite"` (see Details). Default: `"t-as-z"`.

- min_pvalue:

  Float, used to avoid pvalues == 0L for which log10 is Inf. Set this to
  0L if no trimming should be applied. Default: `2^-149` ( ==
  1.401298e-45) the smallest positive subnormal float.

## Value

A list containing the pooled coefficients, SEs, t- and p- values.

## Note

Used inside
[`run_vw_lmm`](https://seredef.github.io/verywise/reference/run_vw_lmm.md).

## References

Rubin, D.B. (1987). *Multiple Imputation for Nonresponse in Surveys*.
New York: John Wiley and Sons.

## Author

Serena Defina, 2024.
