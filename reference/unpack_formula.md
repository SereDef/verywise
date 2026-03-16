# Unpack R formula

Get the terms for the (fixed) effect for the given formula. Wors for
`lme4` as well as regular `lm` style formulas.

## Usage

``` r
unpack_formula(formula, df, return_X = FALSE)
```

## Arguments

- formula:

  A formula object, for example specifying a LME model.

- df:

  A data.frame, for example the first element in the list output of
  [`imp2list`](https://seredef.github.io/verywise/reference/imp2list.md).

- return_X:

  Logical. Whether to return the whole design matrix or just the (fixed)
  term names.

## Value

A design matrix or a character vector of (fixed) term names.

## Author

Serena Defina, 2026.
