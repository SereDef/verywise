# Assert that suggested packages are installed

Call at the top of any function that uses packages listed in
`Suggests:`. Throws an informative error via
[`vw_error()`](https://seredef.github.io/verywise/reference/vw_error.md)
if any are missing, naming the function that requires them so the user
knows why.

## Usage

``` r
require_packages(..., call_fn = "verywise")
```

## Arguments

- ...:

  Package names as unquoted symbols or character strings.

- call_fn:

  Character string; name of the calling function shown in the error
  message. Defaults to the name of the immediate caller.

## Value

`NULL` invisibly if all packages are available.

## Examples

``` r
my_function <- function() {
  require_packages(reticulate, bigreadr)
  # ... rest of function
}
```
