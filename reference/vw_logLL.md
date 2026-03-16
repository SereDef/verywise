# Profile (RE)ML log-likelihood for a single variance ratio

Computes the profile (RE)ML log-likelihood (up to constants) for a given
candidate \\\lambda\\ and a vector of residual sum-of-squares values
\\Q\\ across vertices. This is evaluated inside the \\\lambda\\ grid
search to select the optimal variance ratio per vertex.

## Usage

``` r
vw_logLL(REML = TRUE, df, Q, logdetV, logdetA)
```

## Arguments

- REML:

  Logical. If `TRUE` (default), evaluate the REML objective; otherwise
  the ML objective.

- df:

  Integer. Degrees of freedom: \\N - p\\ (REML) or \\N\\ (ML).

- Q:

  Numeric vector of length \\C\\ (chunk size). Residual sum of squares
  \\Q = Y^\top V^{-1} Y - \hat\beta^\top X^\top V^{-1} Y\\ per vertex.
  Vertices with `Q <= 0` or non-finite values are pre-set to `NA` by the
  caller.

- logdetV:

  Numeric scalar. \\\log \|V(\lambda)\|\\, computed as \\\sum_k \log(1 +
  n_k \lambda)\\.

- logdetA:

  Numeric scalar. \\\log \|A(\lambda)\|\\, computed as \\2 \sum_j \log
  R\_{jj}\\ from the Cholesky factor.

## Value

Numeric vector of length \\C\\; the profile log-likelihood value at this
\\\lambda\\ for each vertex in the chunk. Vertices with `NA` in `Q`
propagate `NA`.

## Details

**REML** profile log-likelihood (constants dropped):
\$\$\ell\_{\text{REML}} = -\frac{1}{2} \left\[ (N-p)
\log\\\frac{Q}{N-p} + \log\|V\| + \log\|A\| \right\]\$\$

**ML** profile log-likelihood (constants dropped): \$\$\ell\_{\text{ML}}
= -\frac{1}{2} \left\[ N \log\\\frac{Q}{N} + \log\|V\| \right\]\$\$

where \\\log\|V\|\\ and \\\log\|A\|\\ are the pre-computed
log-determinants passed as `logdetV` and `logdetA`.
