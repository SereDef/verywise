# Vertex-wise Random-Effects Meta-Analysis Across Studies

This function runs a vertex-wise random-effects meta-analysis across
multiple studies (sites or cohorts). Is expects individual study results
as they are outputted but
[`verywise::run_vw_lmm`](https://seredef.github.io/verywise/reference/run_vw_lmm.md)
or by `QDECR`.

## Usage

``` r
run_vw_meta(
  term,
  hemi = c("lh", "rh"),
  measure = "area",
  study_names,
  study_weights = NULL,
  res_dirs,
  outp_dir = NULL,
  fs_template = "fsaverage",
  seed = 3108,
  n_cores = 1,
  chunk_size = 1000,
  verbose = TRUE
)
```

## Arguments

- term:

  Character. Name of the model term to visualize (matches entries in
  `stack_names.txt`).

- hemi:

  Character string specifying which hemisphere to analyze. Options:
  `"lh"` (left hemisphere: default), `"rh"` (right hemisphere).

- measure:

  Character. Surface measure, e.g. `'area'`, `'thickness'`, `'volume'`.
  Defaults to `'area'`.

- study_names:

  Character vector of study names (must match the length of `res_dirs`).

- study_weights:

  Numeric vector of study weights (e.g. sample size)

- res_dirs:

  Character vector. Path to the directories containing vertex-wise
  result files (`*.mgh`) of each study.

- outp_dir:

  Character string specifying the output directory for results. If
  `NULL` (default), creates a "verywise_results" sub-directory in the
  current working directory (not recommended).

- fs_template:

  Character string specifying the FreeSurfer template for vertex
  registration. Options:

  - `"fsaverage"` (default) = 163842 vertices (highest resolution),

  - `"fsaverage6"` = 40962 vertices,

  - `"fsaverage5"` = 10242 vertices,

  - `"fsaverage4"` = 2562 vertices,

  - `"fsaverage3"` = 642 vertices

  Note that lower resolutions should be only used to downsample the
  brain map, for faster model tuning. The final analyses should also run
  using `fs_template = "fsaverage"` to avoid (small) imprecisions in
  vertex registration and smoothing.

- seed:

  Integer specifying the random seed for reproducibility Default: 3108.

- n_cores:

  Integer specifying the number of CPU cores for parallel processing.
  Default: 1.

- chunk_size:

  Integer specifying the number of vertices processed per chunk in
  parallel operations. Larger values use more memory but may be faster.
  Default: 1000.

- verbose:

  Logical. verbose execution (default: TRUE)

## Value

A named list with three `FBM` objects:

- coef:

  Meta-analytic effect size estimates at each vertex.

- se:

  Meta-analytic standard errors at each vertex.

- p:

  Meta-analytic p-values at each vertex.

Additionally, these are exported as MGH files in the working directory.

## Details

For each vertex, the function loads the effect size and standard error
from each study, computes the variance, and runs a random-effects
meta-analysis using
[`rma`](https://wviechtb.github.io/metafor/reference/rma.uni.html).
Results are saved as Filebacked Big Matrices (FBM) and as MGH files for
downstream neuroimaging analysis.

## See also

[`rma`](https://wviechtb.github.io/metafor/reference/rma.uni.html),
[`FBM`](https://privefl.github.io/bigstatsr/reference/FBM-class.html)

## Examples

``` r
if (FALSE) { # dir.exists("study1/results")
run_vw_meta(
  term = "age",
  hemi = "lh",
  measure = "area",
  study_names = c("Study1", "Study2"),
  res_dirs = c("study1/results", "study2/results")
)
}
```
