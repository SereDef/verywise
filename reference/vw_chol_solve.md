# Solve a linear system given a Cholesky factor

Solves \\A X = B\\ for \\X\\ given the upper-triangular Cholesky factor
\\R\\ of \\A\\ (i.e. \\R^\top R = A\\), using forward and backward
substitution. This is numerically preferable to forming \\A^{-1}\\
explicitly.

## Usage

``` r
vw_chol_solve(R, B)
```

## Arguments

- R:

  Upper-triangular matrix. Cholesky factor of \\A\\ as returned by
  [`vw_chol`](https://seredef.github.io/verywise/reference/vw_chol.md)
  or [`chol`](https://rdrr.io/r/base/chol.html).

- B:

  Numeric matrix (or vector). Right-hand side; must have the same number
  of rows as `R`.

## Value

Numeric matrix \\X = A^{-1} B\\ of the same dimensions as `B`.
