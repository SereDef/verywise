# Pooled degrees of freedom calculation

Edited from `mice::barnard.rubin`. This function computes the residual
degrees of freedom for hypothesis testing, as proposed by Barnard &
Rubin (1999). It is used by the
[`vw_pool`](https://seredef.github.io/verywise/reference/vw_pool.md)
function.

## Usage

``` r
barnard.rubin(lambda, m, dfcom = Inf)
```

## Arguments

- lambda:

  : Proportion of total variance due to missingness

- m:

  : Number of imputed datasets

- dfcom:

  : Residual degrees of freedom, estimated using
  [`stats::df.residual()`](https://rdrr.io/r/stats/df.residual.html)

## Value

Value for the residual degrees of freedom

## References

Barnard, J. and Rubin, D.B. (1999). Small sample degrees of freedom with
multiple imputation. *Biometrika*, 86, 948-955.
