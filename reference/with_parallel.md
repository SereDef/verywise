# Run a chunked foreach loop sequentially or in parallel

Helper to execute a vertex-chunk loop either sequentially (e.g. for
profiling) or in parallel using `foreach` and `doParallel`. The
expression passed in `expr` is evaluated for each chunk in `chuck_seq`.

## Usage

``` r
with_parallel(n_cores, expr, progress_file = NULL, seed = 3108, verbose = TRUE)
```

## Arguments

- n_cores:

  Integer number of CPU cores to use. Values greater than 1 enable
  parallel execution via a fork cluster; 1 forces sequential execution
  via `registerDoSEQ()`.

- expr:

  An expression to be evaluated in parallel or in sequence.

- progress_file:

  Optional character string specifying a log file path. If `NULL`,
  worker output is discarded.

- seed:

  Optional integer seed for reproducible `%dopar%` loops via
  [`doRNG::registerDoRNG()`](https://rdrr.io/pkg/doRNG/man/registerDoRNG.html).

- verbose:

  Logical; if `TRUE`, prints basic messages about the selected execution
  mode.

## Value

Invisibly returns `NULL`. Side effects are produced by `expr`.

## Details

TMP: On Unix-like systems, a fork cluster is used for parallel execution
(type `"FORK"`). On non-Unix systems this will fail anyay while
FreeSurfer is required.

## Author

Serena Defina, 2026.
