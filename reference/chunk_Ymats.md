# Compute the sum and sum of squares of Y

Calculates the cross-product of the design matrix X with the vertex
chunk Y (XtY), along with the sum (1'Y) and sum of squares (Y'Y) for
each vertex in the chunk.

## Usage

``` r
chunk_Ymats(X, ss_chunk)
```

## Arguments

- X:

  Numeric matrix. The design matrix of dimension \\n \times p\\.

- ss_chunk:

  Numeric matrix. A chunk of the ss matrix of dimension \\n \times B\\,
  where B is the block size (number of vertices in the chunk).

## Value

A list containing:

- XtY:

  A numeric matrix of dimension \\p \times B\\, representing \\X^T Y\\.

- psumsY:

  A numeric matrix of dimension \\2 \times B\\. The first row contains
  the column sums of Y (\\1^T Y\\), and the second row contains the
  column sums of squared Y (\\Y^T Y\\).

## Author

Serena Defina, 2026.
