# Initialize per-chunk progress tracking

Prepare progress-tracking metadata for a given chunk of vertices. When
`verbose = TRUE`, prints a message indicating which chunk is being
processed and returns markers for 25\\ that chunk.

## Usage

``` r
init_progress_tracker(chunk, chunk_seq, verbose)
```

## Arguments

- chunk:

  Integer vector of vertex indices for the current chunk. Expected to be
  an element of `chunk_seq` created by
  [`make_chunk_sequence`](https://seredef.github.io/verywise/reference/make_chunk_sequence.md).

- chunk_seq:

  Full list of chunks (as returned by `make_chunk_sequence`), used to
  compute the total number of chunks.

- verbose:

  Logical; if `TRUE`, prints progress messages.

## Value

If `verbose = TRUE`, a list with elements:

- `chunk_idx_str`: character label of the form `"<i>/<n_chunks>"`.

- `milestone_markers`: named list of integer vertex indices at 25\\

If `verbose = FALSE`, returns `NULL` invisibly.

## Author

Serena Defina, 2026.
