# Print a styled message to the console

Routes each call to the appropriate `cli` alert type based on the
conventional leading-character prefixes used in the existing codebase:

## Usage

``` r
vw_message(..., verbose = TRUE, type = NULL)
```

## Arguments

- ...:

  Concatenated via [`paste0()`](https://rdrr.io/r/base/paste.html) to
  form the message string.

- verbose:

  Logical (default `TRUE`).

- type:

  Optional override: `"info"`, `"step"`, `"bullet"`, `"warning"`,
  `"success"`, `"note"`, `"sub"`.

## Details

|                  |                        |
|------------------|------------------------|
| Prefix in `msg`  | Rendered as            |
| `" * "` / `"* "` | bullet (`cli_bullets`) |
| `" ! "` / `" !"` | warning alert          |
| `"NOTE:"`        | info alert             |
| anything else    | `cli_text`             |

Use the explicit typed variants (`vw_step`, `vw_bullet`, …) for new call
sites rather than relying on auto-detection.
