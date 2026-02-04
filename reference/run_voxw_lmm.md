# Run voxel-wise linear mixed model using `lme4::lmer()`

This is is the main function for conducting voxel-wise linear mixed
model analyses on brain morphology. It will first check use inputs,
prepare the phenotype data(list) and run a linear mixed model at each
voxel using the
[`single_lmm`](https://seredef.github.io/verywise/reference/single_lmm.md)
function.

The function supports analysis of both single and multiple imputed
datasets. –TODO– It also automatically handles roi masking, and provides
cluster-wise correction for multiple testing using cluster size
permutation.

## Usage

``` r
run_voxw_lmm(
  formula,
  pheno,
  ss_file,
  outp_dir = NULL,
  brain_template = NULL,
  apply_mask = NULL,
  weights = NULL,
  lmm_control = lme4::lmerControl(),
  seed = 3108,
  n_cores = 1,
  chunk_size = 1000,
  save_residuals = FALSE,
  verbose = TRUE
)
```

## Arguments

- formula:

  A model formula object. This should specify a linear mixed model
  `lme4` syntax. The outcome variable should be one of the supported
  brain surface metrics (see Details). Example:
  `vw_value ~ age * sex + site + (1|participant_id)`.

- pheno:

  Either a `data.frame`/`tibble` containing the "phenotype" data (i.e.,
  already loaded in the global environment), or a string specifying the
  file path to phenotype data. Supported file formats: .rds, .csv, .txt,
  .sav (SPSS). The data should be in **long** format and it should
  contain all the variables specified in the left-hand side of the
  `formula` (i.e., after the `~`) –TODO: no checking done so far– plus
  the `obs_id` column.

- ss_file:

  A path to the super-subject matrix file (observations by voxels). The
  rows are assumed in the same order as the phenotype. .csv formats are
  currently supported.

- outp_dir:

  Character string specifying the output directory for results. If
  `NULL` (default), creates a "results" sub-directory in the current
  working directory (not recommended).

- brain_template:

  Character string specifying the brain template for voxel registration.
  TMP: number of voxels. Options: –TODO–

- apply_mask:

  Logical vector for ROI masking –TODO–

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

- seed:

  Integer specifying the random seed for reproducibility Default: 3108.

- n_cores:

  Integer specifying the number of CPU cores for parallel processing.
  Default: 1.

- chunk_size:

  Integer specifying the number of vertices processed per chunk in
  parallel operations. Larger values use more memory but may be faster.
  Default: 1000.

- save_residuals:

  Logical indicating whether to save the residuals.csv file. Default:
  `FALSE`.

- verbose:

  Logical indicating whether to display progress messages. Default:
  `TRUE`.

  **Statistical Approach:** The function uses
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) for
  mixed-effects modeling, enabling analysis of longitudinal and
  hierarchical data. P-values are computed using the t-as-z
  approximation, with cluster-wise correction applied using –TODO–.

  **Multiple Imputation:** The function automatically detects and
  handles multiple imputed datasets (created with `mice` or similar
  packages), pooling results according to Rubin's rules.

  **Parallel processing:** The `verywise` package employs a carefully
  designed parallelization strategy to maximize computational efficiency
  while avoiding the performance penalties associated with nested
  parallelization. Parallel processing of multiple models should be
  handled by the user (e.g., using SLURM job arrays or similar, see
  vignette on parallelisation). Within a single run_voxw_lmm call,
  voxels are divided into chunks of size `chunk_size` and processed in
  parallel across `n_cores` workers (when `n_cores > 1`). When multiple
  imputed datasets are present, these are processed sequentially within
  each voxel.

  Note that, on some systems, implicit parallelism in low-level matrix
  algebra libraries (BLAS/LAPACK) can interfere with explicit
  parallelization. If you feel like processing is taking too long, I
  recommend disabling these implicit threading libraries before
  starting R. For example:

      export OPENBLAS_NUM_THREADS=1
      export OMP_NUM_THREADS=1
      export MKL_NUM_THREADS=1
      export VECLIB_MAXIMUM_THREADS=1
      export NUMEXPR_NUM_THREADS=1

  Also note that using a very large number of cores (e.g. \>120) may
  sometimes cause worker initialization or other issues (e.g. R parallel
  processes limits)

  **Output Files:** Results are saved in –TODO– format for visualization
  with –TODO–

## Value

A list of file-backed matrices
([`bigstatsr::FBM`](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
objects) containing pooled coefficients, SEs, t- and p- values and
residuals. Results are also automatically saved to disk in –TODO–
format.

## Note

- Large datasets may require substantial memory. Consider adjusting
  `chunk_size` and `n_cores` based on your system specifications.

- For reproducibility, always specify a `seed`.

## See also

[`single_lmm`](https://seredef.github.io/verywise/reference/single_lmm.md)
for single-voxel modeling

## Author

Serena Defina, 2026.
