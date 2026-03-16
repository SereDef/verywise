# Robust Cholesky decomposition with positive-definiteness check

A safe wrapper around [`chol`](https://rdrr.io/r/base/chol.html) that
returns `NULL` instead of throwing an error when the input matrix `A` is
not numerically positive-definite (PD). An additional near-singularity
check rejects decompositions whose smallest diagonal pivot falls below
`tol`, preventing numerically unstable downstream solves.

## Usage

``` r
vw_chol(A, tol = 1e-12)
```

## Arguments

- A:

  Numeric square matrix. Should be symmetric PD (e.g. a Hessian \\X^\top
  V^{-1} X\\).

- tol:

  Numeric. Minimum acceptable diagonal pivot. Default: `1e-12`.

## Value

The upper-triangular Cholesky factor \\R\\ such that \\R^\top R = A\\,
or `NULL` if `A` is not sufficiently PD.
