# Update within-chunk progress for a (milestone) vertex

Check whether the current vertex index matches any pre-defined
within-chunk milestones (25\\ message if so.

## Usage

``` r
update_progress_tracker(v, progress_tracker, verbose)
```

## Arguments

- v:

  Integer; current vertex index.

- progress_tracker:

  List returned by
  [`init_progress_tracker`](https://seredef.github.io/verywise/reference/init_progress_tracker.md),
  containing `chunk_idx_str` and `milestone_markers`.

- verbose:

  Logical; if `TRUE`, prints progress messages.

## Value

Invisibly returns `NULL`. Used for its side effects (messages).

## Author

Serena Defina, 2026.
