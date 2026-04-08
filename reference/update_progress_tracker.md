# Update within-chunk progress for a (milestone) vertex

Check whether the current vertex index matches any pre-defined
within-chunk milestones (25\\ message if so.

## Usage

``` r
update_progress_tracker(
  v,
  progress_tracker,
  progress_file = "",
  verbose = TRUE
)
```

## Arguments

- v:

  Integer; current vertex index.

- progress_tracker:

  List returned by
  [`init_progress_tracker`](https://seredef.github.io/verywise/reference/init_progress_tracker.md),
  containing `chunk_idx_str` and `milestone_markers`.

- progress_file:

  Optional character string specifying a log file path. If `NULL`,
  worker output is discarded.

- verbose:

  Logical; if `TRUE`, prints progress messages.

## Value

Invisibly returns `NULL`. Used for its side effects (messages).

## Author

Serena Defina, 2026.
