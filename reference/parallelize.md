# Run a chunked foreach loop sequentially or in parallel

Helper to execute a vertex-chunk loop either sequentially (e.g. for
profiling) or in parallel using `foreach` and `doParallel`. The
expression passed in `expr` is evaluated for each chunk in `chuck_seq`.

## Usage

``` r
parallelize(
  n_cores,
  chunk_seq,
  progress_file = NULL,
  packages = NULL,
  exports = NULL,
  seed = NULL,
  verbose = TRUE,
  expr
)
```

## Arguments

- n_cores:

  Integer number of CPU cores to use. Values greater than 1 enable
  parallel execution via a fork cluster; 1 forces sequential execution
  via `registerDoSEQ()`.

- chunk_seq:

  List of integer vectors produced by
  [`make_chunk_sequence`](https://seredef.github.io/verywise/reference/make_chunk_sequence.md).

- progress_file:

  Optional character string specifying a log file path. If `NULL`,
  worker output is discarded.

- packages:

  Character vector of package names to make available on workers (passed
  to `foreach`'s `.packages`).

- exports:

  Character vector of object names to export to workers (passed to
  `foreach`'s `.export`).

- seed:

  Optional integer seed for reproducible `%dopar%` loops via
  [`doRNG::registerDoRNG()`](https://rdrr.io/pkg/doRNG/man/registerDoRNG.html).

- verbose:

  Logical; if `TRUE`, prints basic messages about the selected execution
  mode.

- expr:

  An expression to be evaluated for each chunk. Inside `expr`, the
  symbol `chunk` refers to the current element of `chunk_seq`.

## Value

Invisibly returns `NULL`. Side effects are produced by `expr`.

## Details

TMP: On Unix-like systems, a fork cluster is used for parallel execution
(type `"FORK"`). On non-Unix systems this will fail anyay while
FreeSurfer is required.

## Author

Serena Defina, 2026.
