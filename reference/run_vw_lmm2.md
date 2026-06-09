# Run vertex-wise linear mixed model using [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)

This is an alternative to the main function for conducting vertex-wise
linear mixed model analyses on brain surface metrics. This is almost
identical to
[`run_vw_lmm()`](https://seredef.github.io/verywise/reference/run_vw_lmm.md)
but it estimates a linear mixed model from scratch at each vertex
(instead of refitting from a template) This is done using the
[`single_lmm()`](https://seredef.github.io/verywise/reference/single_lmm.md)
function.

Use this pipeline when you have reason to believe that the random effect
structure is very different across vertices. In all test I have run so
far,
[`run_vw_lmm()`](https://seredef.github.io/verywise/reference/run_vw_lmm.md)
and `run_vw_lmm2()` give identical fixed term results.

## Usage

``` r
run_vw_lmm2(
  formula,
  pheno,
  subj_dir,
  outp_dir = NULL,
  hemi = c("lh", "rh"),
  fs_template = "fsaverage",
  apply_cortical_mask = TRUE,
  folder_id = "folder_id",
  tolerate_surf_not_found = 20,
  weights = NULL,
  lmm_control = lme4::lmerControl(calc.derivs = FALSE),
  REML = TRUE,
  seed = 3108,
  n_cores = 1,
  chunk_size = 1000,
  FS_HOME = Sys.getenv("FREESURFER_HOME"),
  fwhm = 10,
  mcz_thr = 30,
  cwp_thr = 0.025,
  save_optional_cluster_info = FALSE,
  save_ss = FALSE,
  save_residuals = FALSE,
  verbose = TRUE
)
```

## Arguments

- formula:

  A model formula object. This should specify a linear mixed model
  `lme4` syntax. The outcome variable should be one of the supported
  brain surface metrics (see Details). Example:
  `vw_thickness ~ age * sex + site + (1|participant_id)`.

- pheno:

  Either a `data.frame`/`tibble` containing the "phenotype" data (i.e.,
  already loaded in the global environment), or a string specifying the
  file path to phenotype data. Supported file formats: .rds, .csv, .txt,
  .sav (SPSS). The data should be in **long** format and it should
  contain all the variables specified in the left-hand side of the
  `formula` (i.e., after the `~`) plus the `folder_id` column.

- subj_dir:

  Character string specifying the path to FreeSurfer data directory.
  Must follow the verywise directory structure (see package vignette for
  details).

- outp_dir:

  Character string specifying the output directory for results. If
  `NULL` (default), creates a "verywise_results" sub-directory in the
  current working directory (not recommended).

- hemi:

  Character string specifying which hemisphere to analyze. Options:
  `"lh"` (left hemisphere: default), `"rh"` (right hemisphere).

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

- apply_cortical_mask:

  Logical indicating whether to exclude non-cortical vertices from
  analysis. Default: `TRUE` (recommended).

- folder_id:

  Character string specifying the column name in `pheno` that contains
  subject directory names of the input neuroimaging data (e.g.
  "sub-001_ses-baseline" or "site1/sub-010_ses-F1"). These are expected
  to be nested inside `subj_dir`. Default: `"folder_id"`.

- tolerate_surf_not_found:

  Integer indicating how many brain surface files listed in `folder_id`
  can be missing from `subj_dir`. If the number of missing or corrupted
  files is ` > tolerate_surf_not_found ` execution will stop. Default:
  `20L`.

- weights:

  Optional string or numeric vector of weights for the linear mixed
  model. You can use this argument to specify inverse-probability
  weights. If this is a string, the function look for a column with that
  name in the phenotype data. Note that these are not normalized or
  standardized in any way. Default: `NULL` (no weights).

- lmm_control:

  Optional list (of correct class, resulting from `lmerControl()`
  containing control parameters to be passed to
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) (e.g.
  optimizer choice, convergence criteria, see the `?lmerControl`
  documentation for details. Default: (uses default settings).

- REML:

  Logical specifying whether to optimize the REML criterion (as opposed
  to the log-likelihood). Default: TRUE. Use `REML = FALSE` if you
  intend to do model comparison (using AIC output).

- seed:

  Integer specifying the random seed for reproducibility Default: 3108.

- n_cores:

  Integer specifying the number of CPU cores for parallel processing.
  Default: 1.

- chunk_size:

  Integer specifying the number of vertices processed per chunk in
  parallel operations. Larger values use more memory but may be faster.
  Default: 1000.

- FS_HOME:

  Character string specifying the FreeSurfer home directory. Defaults to
  `FREESURFER_HOME` environment variable.

- fwhm:

  Numeric value specifying the full-width half-maximum for smoothing
  kernel. Default: 10.

- mcz_thr:

  Numeric value for the Monte Carlo simulation threshold. Any of the
  following are accepted (equivalent values are separated by `/`):

  - 13 / 1.3 / 0.05,

  - 20 / 2.0 / 0.01,

  - 23 / 2.3 / 0.005,

  - 30 / 3.0 / 0.001, \\ default

  - 33 / 3.3 / 0.0005,

  - 40 / 4.0 / 0.0001.

- cwp_thr:

  Numeric value for cluster-wise p-value threshold (on top of all
  corrections). Set this should be set to `0.025` when both hemispheres
  are analyzed, and `0.05` for single hemisphere analyses. Default:
  `0.025`.

- save_optional_cluster_info:

  Logical indicating whether to save additional output form
  `mri_surfcluster` call. See
  [`compute_clusters`](https://seredef.github.io/verywise/reference/compute_clusters.md)
  for details. Default: `FALSE`.

- save_ss:

  Logical indicating whether to save the super-subject matrix ("ss") as
  an .rds file that can be then re-used in future analyses. This can
  also be a character string specifying the directory where ss should be
  saved. When `TRUE`, the ss matrix will be saved in `<outp_dir>/ss` by
  default. Default: `FALSE`.

- save_residuals:

  Logical indicating whether to save the residuals.mgh file. Default:
  `FALSE`.

- verbose:

  Logical indicating whether to display progress messages. Default:
  `TRUE`.

## Details

See
[`run_vw_lmm()`](https://seredef.github.io/verywise/reference/run_vw_lmm.md)
for details.

## Author

Serena Defina, 2026.
