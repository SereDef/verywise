# Convert imputation object to a list of dataframes

Converts known imputation objects (e.g., `amelia`, `aregImpute`, `mi`,
`mids`, `missForest`) to a `list` of dataframes.

## Usage

``` r
imp2list(obj)
```

## Arguments

- obj:

  The imputation object that contains the imputed datasets

## Value

A list of imputed datasets

## Details

This function attempts to unify imputation object formats in supplying
datasets. This is done by extracting the imputed datasets from the
imputation object and assembling them into a `list`. If an unknown
object is supplied, the function will throw an error with the request to
supply a list of data frames instead.

This is a generic function: methods can be defined for it directly, see
`methods(imp2list)` for a list of all available methods.

## Author

Sander Lamballais, 2018.
