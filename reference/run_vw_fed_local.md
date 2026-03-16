# Run vertex-wise linear model on local site data

This the "local" function for conducting a distributed linear model
analyses on brain surface metrics. It will first check use inputs,
prepare the phenotype data and compute vertwx-wise sufficient statistics
for the specified hemisphere. These can be then compressed, shared and
finally aggregatedusing using the
[`run_vw_fed_aggr`](https://seredef.github.io/verywise/reference/run_vw_fed_aggr.md)
function.

## Usage

``` r
run_vw_fed_local(
  site_name,
  formula,
  pheno,
  subj_dir,
  outp_dir = NULL,
  hemi = c("lh", "rh"),
  fs_template = "fsaverage",
  apply_cortical_mask = TRUE,
  folder_id = "folder_id",
  tolerate_surf_not_found = 20,
  fwhm = 10,
  seed = 3108,
  n_cores = 1,
  chunk_size = 1000,
  save_ss = FALSE,
  verbose = TRUE
)
```

## Arguments

- site_name:

  A character string indicating the site name or ID.

- formula:

  A model formula object. This should specify a linear model. The
  outcome variable should be one of the supported brain surface metrics
  (see Details). Example: `vw_thickness ~ age * sex + ethinicy`.

- pheno:

  Either a `data.frame`/`tibble` containing the "phenotype" data (i.e.,
  already loaded in the global environment), or a string specifying the
  file path to phenotype data. Supported file formats: .rds, .csv, .txt,
  .sav (SPSS). The data contain all the variables specified in the
  left-hand side of the `formula` (i.e., after the `~`) plus the
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

- fwhm:

  Numeric value specifying the full-width half-maximum for smoothing
  kernel. This is used to read FreeSurfer files. Default: 10.

- seed:

  Integer specifying the random seed for reproducibility Default: 3108.

- n_cores:

  Integer specifying the number of CPU cores for parallel processing.
  Default: 1.

- chunk_size:

  Integer specifying the number of vertices processed per chunk in
  parallel operations. Larger values use more memory but may be faster.
  Default: 1000.

- save_ss:

  Logical indicating whether to save the super-subject matrix ("ss") as
  an .rds file that can be then re-used in future analyses. This can
  also be a character string specifying the directory where ss should be
  saved. When `TRUE`, the ss matrix will be saved in `<outp_dir>/ss` by
  default. Default: `FALSE`.

- verbose:

  Logical indicating whether to display progress messages. Default:
  `TRUE`.

## Value

A list of site-speficic information summary matrices, of which some are
([`bigstatsr::FBM`](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
objects). These should be compressed before sending them to the
aggregation center.

## Details

The function does not currently support multiple imputed datasets or IPW
weights (this is for future developement)

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

**Parallel processing:** The `verywise` package employs a carefully
designed parallelization strategy to maximize computational efficiency
while avoiding the performance penalties associated with nested
parallelization. Left and right cortical hemispheres are processed
sequentially by default. Parallel processing of the two hemispheres
(and/or different metrics, models) should be handled by the user (e.g.,
using SLURM job arrays or similar, see vignette on parallelisation).
Within each hemisphere, vertices are divided into chunks of size
`chunk_size` and processed in parallel across `n_cores` workers (when
`n_cores > 1`). When multiple imputed datasets are present, these are
processed sequentially within each vertex.

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

## Note

- Large datasets may require substantial memory. Consider adjusting
  `chunk_size` and `n_cores` based on your system specifications.

- For reproducibility, always specify a `seed`.

## See also

[`chunk_Ymats`](https://seredef.github.io/verywise/reference/chunk_Ymats.md)
for vertex-chunk modeling,

## Author

Serena Defina, 2026.
