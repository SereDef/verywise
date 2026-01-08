# Running vertex-wise meta-analyses

`verywise` also implements vertex-wise meta-analysis via the `metafor`
package

``` r
# Run a meta-analysis
run_vw_meta(
  term = "age", # Which "term" / predictor / effect to pool
  hemi = "lh", # (default) or "rh": which hemisphere to run
  measure = "area", # (default) or any available FreeSurfer metric.
  res_dirs = c("/path/to/study1/results", "/path/to/study2/results"),
  study_names = c("Study1", "Study2"),
  n_cores = 4  # parallel processing
)
```

\[TODO: expand with more details\]
