# Run vertex-wise linear mixed model using `lme4::lmer()`

This is is the main function for conducting vertex-wise linear mixed
model analyses on brain surface metrics. It will first check use inputs,
prepare the phenotype data(list) and run a linear mixed model at each
vertex of the specified hemisphere using the
[`single_lmm`](https://seredef.github.io/verywise/reference/single_lmm.md)
function.

The function supports analysis of both single and multiple imputed
datasets. It also automatically handles cortical masking, and provides
cluster-wise correction for multiple testing using FreeSurfer's Monte
Carlo simulation approach.

## Usage

``` r
run_vw_lmm(
  formula,
  pheno,
  subj_dir,
  outp_dir = NULL,
  hemi = c("lh", "rh"),
  fs_template = "fsaverage",
  apply_cortical_mask = TRUE,
  folder_id = "folder_id",
  tolerate_surf_not_found = 20,
  use_model_template = TRUE,
  weights = NULL,
  lmm_control = lme4::lmerControl(),
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
  contain all the variables specified in the `formula` plus the
  `folder_id` column.

- subj_dir:

  Character string specifying the path to FreeSurfer data directory.
  Must follow the verywise directory structure (see package vignette for
  details).

- outp_dir:

  Character string specifying the output directory for results. If
  `NULL` (default), creates a "verywise_results" sub-directory in the
  current working directory (not recommended).

- hemi:

  Character string specifying which hemisphere(s) to analyze. Options:
  `"both"` (default), `"lh"` (left only), `"rh"` (right only).

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

- use_model_template:

  Logical indicating whether to pre-compile the model template for
  faster estimation. Default: `TRUE` (recommended).

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
  optimizer choice, convergence criteria, see the `*lmerControl`
  documentation for details. Default: (uses default settings).

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

  Numeric value for Monte Carlo simulation threshold. Any of the
  following are accepted (equivalent values separate by `/`):

  - 13 / 1.3 / 0.05,

  - 20 / 2.0 / 0.01,

  - 23 / 2.3 / 0.005,

  - 30 / 3.0 / 0.001, \\ default

  - 33 / 3.3 / 0.0005,

  - 40 / 4.0 / 0.0001.

- cwp_thr:

  Numeric value for cluster-wise p-value threshold (on top of all
  corrections). Set this to 0.025 when both hemispheres are analyzed,
  0.05 for single hemisphere. Default: 0.025.

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

## Value

A list of file-backed matrices
([`bigstatsr::FBM`](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
objects) containing pooled coefficients, SEs, t- and p- values and
residuals. Results are also automatically saved to disk in .mgh format.

## Details

**Supported Brain Surface Metrics:** The outcome specified in `formula`
should be a brain surface metric among:

- `vw_thickness` - Cortical thickness

- `vw_area` - Cortical surface area (white surface)

- `vw_area.pial` - Cortical surface area (pial surface)

- `vw_curv` - Mean curvature

- `vw_jacobian_white` - Jacobian determinant (white surface)

- `vw_pial` - Pial surface coordinates

- `vw_pial_lgi` - Local gyrification index (pial surface)

- `vw_sulc` - Sulcal depth

- `vw_volume` - Gray matter volume

- `vw_w_g.pct` - White/gray matter intensity ratio

- `vw_white.H` - Mean curvature (white surface)

- `vw_white.K` - Gaussian curvature (white surface)

**Statistical Approach:** The function uses
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) for
mixed-effects modeling, enabling analysis of longitudinal and
hierarchical data. P-values are computed using the t-as-z approximation,
with cluster-wise correction applied using FreeSurfer's Monte Carlo
simulation approach.

**Multiple Imputation:** The function automatically detects and handles
multiple imputed datasets (created with `mice` or similar packages),
pooling results according to Rubin's rules.

**Parallel processing:** The `verywise` package employs a carefully
designed parallelization strategy to maximize computational efficiency
while avoiding the performance penalties associated with nested
parallelization. Left and right cortical hemispheres are processed
sequentially by default. Parallel processing of the two hemispheres
(and/or different metrics, models) should be handled by the user (e.g.,
using SLURM job arrays or similar, see vignette on parallelisation
...COMING UP). Within each hemisphere, vertices are divided into chunks
of size `chunk_size` and processed in parallel across `n_cores` workers
(when `n_cores > 1`). When multiple imputed datasets are present, these
are processed sequentially within each vertex.

Note that, on some systems, implicit parallelism in low-level matrix
algebra libraries (BLAS/LAPACK) can interfere with explicit
parallelization. If you feel like processing is taking too long, I
recommend disabling these implicit threading libraries before starting
R. For example:

    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export VECLIB_MAXIMUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1

Also note that using a very large number of cores (e.g. \>120) may
sometimes cause worker initialization or other issues (e.g. R parallel
processes limits)

**Output Files:** Results are saved in FreeSurfer-compatible .mgh format
for visualization with
[verywiseWIZard](https://github.com/SereDef/verywise-wizard), FreeView
or other neuroimaging software.

## Note

- Ensure FreeSurfer is properly installed and the `FREESURFER_HOME`
  environment variable is set.

- Large datasets may require substantial memory. Consider adjusting
  `chunk_size` and `n_cores` based on your system specifications.

- For reproducibility, always specify a `seed`.

## See also

[`single_lmm`](https://seredef.github.io/verywise/reference/single_lmm.md)
for single-vertex modeling,
`vignette("03-run-vw-lmm", package = "verywise")` for detailed usage
examples.

## Author

Serena Defina, 2024.
