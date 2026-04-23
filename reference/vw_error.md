# Themed error

Wrapper around
[`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)
that applies `VW_THEME` and reports the call from the correct frame.
Accepts the same `msg` conventions as
[`vw_message()`](https://seredef.github.io/verywise/reference/vw_message.md):
a plain string, a named character vector (for bullet-style bodies), or
inline cli markup.

## Usage

``` r
vw_error(msg, ..., call = rlang::caller_env())
```

## Arguments

- msg:

  A character string or named character vector passed to
  [`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html).
  Use names `"i"`, `"x"`, `"v"`, `">"` for bullet formatting, e.g.
  `c("Something failed.", "i" = "Check {.arg x}.")`.

- ...:

  Additional arguments forwarded to
  [`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)
  (e.g. `class`, `call`, `.internal`).

## Author

Serena Defina, 2026.
