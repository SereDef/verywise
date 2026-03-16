# Run federated aggregation for vertex-wise distributed LMM

Performs the aggregation step of a privacy-preserving federated linear
mixed model (LMM) on vertex-wise neuroimaging data. This function never
sees raw participant data — it operates entirely on pre-computed
sufficient statistics (site payloads) uploaded by each site. It fits a
random-intercept LMM with a single shared variance component ratio
\\\lambda = \sigma_b^2 / \sigma\_\varepsilon^2\\, selected via a 1-D
profile likelihood grid search (rather than iterative site
communication).

## Usage

``` r
run_vw_fed_aggr(
  site_names,
  formula,
  inpt_dir,
  outp_dir = NULL,
  hemi = c("lh", "rh"),
  fs_template = "fsaverage",
  fwhm = 10,
  seed = 3108,
  n_cores = 1,
  chunk_size = 1000,
  REML = TRUE,
  lambda_grid = NULL,
  chol_tol = 1e-12,
  ridge = 0,
  verbose = TRUE
)
```

## Arguments

- site_names:

  Character vector of site identifiers. Each site must have a
  corresponding `<site>.tar.gz` archive in `inpt_dir`.

- formula:

  A two-sided [`formula`](https://rdrr.io/r/stats/formula.html)
  specifying the fixed effects. The LHS names the surface measure (e.g.
  `thickness ~ age + sex`).

- inpt_dir:

  Character. Directory containing the site payload archives.

- outp_dir:

  Character or `NULL`. Directory for output FBM backing files. Defaults
  to `inpt_dir` if `NULL`.

- hemi:

  Character. Hemisphere: `"lh"` (default) or `"rh"`.

- fs_template:

  Character. FreeSurfer template surface used for the cortical mask
  (e.g. `"fsaverage"`, `"fsaverage5"`). Default: `"fsaverage"`.

- fwhm:

  Numeric. Full-width at half-maximum (mm) of the smoothing kernel
  applied at the site level. Used for informational messages only;
  smoothing must be applied before payload generation. Default: `10`.

- seed:

  Integer. Random seed passed to the parallel backend for
  reproducibility. Default: `3108`.

- n_cores:

  Integer. Number of parallel workers. `1` (default) runs sequentially;
  values `> 1` use doParallel.

- chunk_size:

  Integer. Number of vertices processed per parallel job. Larger values
  are faster but use more memory. Default: `1000`.

- REML:

  Logical. If `TRUE` (default), use restricted maximum likelihood
  (REML); otherwise use ML.

- lambda_grid:

  Numeric vector of candidate \\\lambda\\ values for the grid search, or
  `NULL` to use the default log-spaced grid spanning \\\[10^{-6},
  10^{2}\]\\.

- chol_tol:

  Numeric. Minimum diagonal pivot in the Cholesky factor below which a
  \\\lambda\\ value is rejected as numerically singular. Default:
  `1e-12`.

- ridge:

  Numeric. Ridge penalty added to the diagonal of \\A(\lambda)\\ before
  factorisation. Use a small positive value (e.g. `1e-6`) when
  predictors are near-collinear. Default: `0` (no regularisation).

- verbose:

  Logical. Whether to print progress messages. Default: `TRUE`.

## Value

A named list of four
[`FBM`](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
objects, each of dimension described below. Rows correspond to model
terms; columns to vertices.

- `coef`:

  \\p \times V\\ matrix of fixed-effect estimates \\\hat\beta\\.

- `se`:

  \\p \times V\\ matrix of standard errors.

- `pval`:

  \\p \times V\\ matrix of two-sided p-values (t-distribution, `df`
  degrees of freedom).

- `variances`:

  \\2 \times V\\ matrix; row 1 = \\\hat\sigma^2\_\varepsilon\\, row 2 =
  \\\hat\tau^2 = \hat\lambda \hat\sigma^2\_\varepsilon\\.

## Details

**Model structure.** The covariance matrix is assumed to follow \$\$V =
\sigma\_\varepsilon^2 (I + \lambda Z Z^\top)\$\$ where \\Z\\ is a
block-diagonal indicator of site membership. Under this structure, the
Sherman–Morrison–Woodbury identity reduces the \\N \times N\\ inverse to
site-level rank-1 updates, making the algorithm \\O(K p^2)\\ rather than
\\O(N^3)\\.

**Grid search.** For each candidate \\\lambda\\, the function
reconstructs the global Hessian \\A(\lambda) = X^\top V^{-1} X\\ and its
Cholesky factor, then computes the profile (RE)ML log-likelihood per
vertex in chunks. The \\\lambda\\ that maximises the log-likelihood is
selected independently for each vertex.

**Chunked processing.** Vertices are processed in chunks of size
`chunk_size` to keep memory usage bounded regardless of surface
resolution. Results are written directly to memory-mapped bigstatsr FBM
files on disk.

## See also

[`chunk_dlmm`](https://seredef.github.io/verywise/reference/chunk_dlmm.md),
[`decompress_site_payloads`](https://seredef.github.io/verywise/reference/decompress_site_payloads.md),
[`vw_logLL`](https://seredef.github.io/verywise/reference/vw_logLL.md)

## Author

Serena Defina, 2026.
