# Compute significant clusters of vertices

This function is a wrapper for the FreeSurfer command `mri_surfcluster`.
It computes the clusters of significant vertices.

## Usage

``` r
compute_clusters(
  stack_path,
  hemi,
  fwhm,
  FS_HOME,
  cwp_thr = 0.025,
  mcz_thr = 30,
  csd_sign = "abs",
  full_surfcluster_output = FALSE,
  mask = NULL,
  verbose = FALSE
)
```

## Arguments

- stack_path:

  : path where this term's p values are stored and where results are
  generated.

- hemi:

  : hemisphere.

- fwhm:

  : full-width half maximum estimate of data smoothness.

- FS_HOME:

  : FreeSurfer directory, i.e. `$FREESURFER_HOME`

- cwp_thr:

  : (default = 0.025, when both hemispheres are ran, else 0.05) the
  cluster-wise p-value threshold on top of all corrections.

- mcz_thr:

  : (default = 0.001) numeric value for the Monte Carlo simulation
  threshold. Any of the following are accepted (equivalent values
  separate by `/`):

  - 13 / 1.3 / 0.05,

  - 20 / 2.0 / 0.01, . \* 23 / 2.3 / 0.005,

  - 30 / 3.0 / 0.001, \\ default

  - 33 / 3.3 / 0.0005,

  - 40 / 4.0 / 0.0001.

- csd_sign:

  : (default = "abs")

- full_surfcluster_output:

  : (default = FALSE) whether to save additional files: `.cluster.mgh`:
  a map of cluster-wise significances files: `.voxel.mgh`: a map of
  corrected voxel-wise significance files: `.ocn.annot`: output clusters
  as an annotation files: `.masked.mgh`: input with non-clusters set to
  0

- mask:

  : apply a mask (default = cortical mask + excluding problematic
  vertices)

- verbose:

  : (default = FALSE) verbosity.
