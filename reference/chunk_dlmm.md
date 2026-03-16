# Vertex-wise distributed LMM: per-chunk estimation

Core computational workhorse called by
[`run_vw_fed_aggr`](https://seredef.github.io/verywise/reference/run_vw_fed_aggr.md)
for each chunk of vertices. Given the pre-computed static lambda profile
`lambda_prof` and the chunk-specific sufficient statistics (\\X^\top
Y\\, \\\mathbf{1}^\top Y\\, \\Y^\top Y\\), this function selects the
optimal \\\lambda\\ per vertex and returns fixed-effect estimates,
standard errors, p-values, and variance components.

## Usage

``` r
chunk_dlmm(
  lambda_prof,
  X1,
  XtY,
  sumY,
  YtY,
  n_terms,
  chunk_size,
  n_sites,
  REML,
  df
)
```

## Arguments

- lambda_prof:

  List. Output of the static lambda profile loop in
  [`run_vw_fed_aggr`](https://seredef.github.io/verywise/reference/run_vw_fed_aggr.md).
  Each element is a named list with fields `lambda`, `w` (site weights),
  `R` (Cholesky factor of \\A\\), `logdetA`, and `logdetV`.

- X1:

  List of length \\K\\. Each element is the column vector of row-sums of
  the design matrix for site \\k\\ (\\X_k^\top \mathbf{1}\_{n_k}\\).

- XtY:

  List of length \\K\\. Each element is a \\p \times C\\ matrix
  \\X_k^\top Y_k\\ for the current chunk of \\C\\ vertices.

- sumY:

  List of length \\K\\. Each element is a length-\\C\\ vector
  \\\mathbf{1}^\top Y_k\\ (column sums of the outcome chunk).

- YtY:

  List of length \\K\\. Each element is a length-\\C\\ vector of
  element-wise squared column sums \\Y_k^\top Y_k\\.

- n_terms:

  Integer. Number of fixed-effect terms \\p\\.

- chunk_size:

  Integer. Number of vertices in the current chunk \\C\\.

- n_sites:

  Integer. Number of sites \\K\\.

- REML:

  Logical. Use REML (`TRUE`) or ML (`FALSE`) likelihood.

- df:

  Integer. Residual degrees of freedom (\\N - p\\ for REML, \\N\\ for
  ML).

## Value

A named list with four matrices, all with \\C\\ columns:

- `coef`:

  \\p \times C\\ fixed-effect estimates.

- `se`:

  \\p \times C\\ standard errors.

- `pval`:

  \\p \times C\\ two-sided p-values.

- `bw_var`:

  \\2 \times C\\ variance components; row 1 =
  \\\hat\sigma^2\_\varepsilon\\, row 2 = \\\hat\tau^2\\.

## Details

For each candidate \\\lambda\\ the dynamic quantities are:
\$\$B(\lambda) = X^\top V^{-1} Y = \sum_k \bigl(X_k^\top Y_k - w_k
\mathbf{x}\_{k1} \mathbf{1}\_k^\top Y_k \bigr)\$\$ \$\$C(\lambda) =
Y^\top V^{-1} Y = \sum_k \bigl(Y_k^\top Y_k - w_k (\mathbf{1}\_k^\top
Y_k)^2 \bigr)\$\$ where \\w_k = \lambda / (1 + n_k \lambda)\\ is the
site weight. The profile (RE)ML log-likelihood is evaluated at each
candidate \\\lambda\\ and the maximiser selected per vertex via
[`max.col`](https://rdrr.io/r/base/maxCol.html).

Standard errors are derived from \\\text{Var}(\hat\beta) = \hat\sigma^2
A^{-1}\\, with \\A^{-1}\\ obtained via
[`chol2inv`](https://rdrr.io/r/base/chol2inv.html) from the pre-computed
Cholesky factor.
