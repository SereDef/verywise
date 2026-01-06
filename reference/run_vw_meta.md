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
  res_dirs,
  study_names,
  fs_template = "fsaverage",
  n_cores = 1,
  verbose = TRUE
)
```

## Arguments

- term:

  Character. Name of the model term to visualize (matches entries in
  `stack_names.txt`).

- hemi:

  Character. Hemisphere to analyse (`"lh"` or `"rh"`).

- measure:

  Character. Surface measure, e.g. `'area'`, `'thickness'`, `'volume'`.
  Defaults to `'area'`.

- res_dirs:

  Character vector. Path to the directories containing vertex-wise
  result files (`*.mgh`) of each study.

- study_names:

  Character vector of study names (must match the length of `res_dirs`).

- fs_template:

  Character string specifying the FreeSurfer template surface. The
  following values are accepted:

  - fsaverage (default) = 163842 vertices (highest resolution),

  - fsaverage6 = 40962 vertices,

  - fsaverage5 = 10242 vertices,

  - fsaverage4 = 2562 vertices,

  - fsaverage3 = 642 vertices

- n_cores:

  Integer. Number of CPU cores to use for parallel processing (default:
  1).

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
if (FALSE) { # \dontrun{
run_vw_meta(
  term = "age",
  hemi = "lh",
  measure = "area",
  res_dirs = c("study1/results", "study2/results"),
  study_names = c("Study1", "Study2"),
  fs_template = "fsaverage",
  n_cores = 4
)
} # }
```
