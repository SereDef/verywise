# Unpack `lme4` formula

Get the terms for the fixed effect for the given formula.

## Usage

``` r
get_terms(formula, dset)
```

## Arguments

- formula:

  : model formula object (this should specify a LME model)

- dset:

  : the data.frame, for example the first element in the list output of
  [`imp2list`](https://seredef.github.io/verywise/reference/imp2list.md).

## Value

A character vector of fixed terms.

## Author

Serena Defina, 2024.
