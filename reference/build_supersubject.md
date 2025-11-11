# Build "supersubject" by stacking all vertex data in one large file-backed matrix with dimensions n_subjects x n_vertices.

Build "supersubject" by stacking all vertex data in one large
file-backed matrix with dimensions n_subjects x n_vertices.

## Usage

``` r
build_supersubject(
  subj_dir,
  folder_ids,
  supsubj_dir,
  measure,
  hemi,
  n_cores,
  fwhmc = "fwhm10",
  fs_template = "fsaverage",
  backing,
  error_cutoff = 20,
  save_rds = FALSE,
  verbose = TRUE
)
```

## Arguments

- subj_dir:

  : path to the FreeSurfer data, this expects a verywise structure.

- folder_ids:

  : the vector of observations to include. This holds the relative path
  (from `subj_dir`) to the FreeSurfer data folder (e.g.
  "site1/sub-1_ses-01").

- supsubj_dir:

  : output path, where logs, backing files and the matrix itself (if
  `save_rds == TRUE`) will be stored.

- measure:

  : vertex-wise measure, used to identify files.

- hemi:

  : hemisphere, used to identify files.

- n_cores:

  : number of cores to use for parallel processing.

- fwhmc:

  : (default = "fwhm10") full-width half maximum value, used to identify
  files.

- fs_template:

  : (default = "fsaverage") template on which to register vertex-wise
  data. The following values are accepted:

  - fsaverage (default) = 163842 vertices (highest resolution),

  - fsaverage6 = 40962 vertices,

  - fsaverage5 = 10242 vertices,

  - fsaverage4 = 2562 vertices,

  - fsaverage3 = 642 vertices Note that, at the moment, these are only
    used to downsample the brain map, for faster model tuning.
    `verywise` expects the input data to be always registered on the
    "fsaverage" template and the final analyses should also be run using
    `fs_template = "fsaverage"` to avoid (small) imprecisions in vertex
    registration and smoothing.

- backing:

  : (default = `supsubj_dir`) location to save the matrix `backingfile`.

- error_cutoff:

  : (default = 20) how many missing directories or brain surface files
  for the function to stop with an error. If \< `error_cutoff`
  directories/files are not found a warning is thrown and missing files
  are registered in the `issues.log` file.

- save_rds:

  : (default = FALSE) save the supersubject file metadata for re-use in
  other sessions.

- verbose:

  : (default = TRUE)

## Value

A Filebacked Big Matrix with vertex data for all subjects (dimensions:
n_subjects x n_vertices)

## Author

Serena Defina, 2024.
