# Simulate longitudinal FreeSurfer vertex-wise data

Writes one `.mgh` file per (subject, session) pair following the
generative model documented in
[`simulate_longit_dataset2`](https://seredef.github.io/verywise/reference/simulate_longit_dataset2.md).
This function is exported so advanced users can generate brain surface
data from a pre-existing phenotype `data.frame` (e.g., after multiple
imputation).

Random effects drawn per call (all seeded via `seed`):

- Vertex-specific intercepts \\\mu_v\\ — drawn once, shared across sites
  and subjects.

- Site random intercepts \\u_k\\ — one per site.

- Subject random intercepts \\b_i\\ and slopes \\b^{\mathrm{time}}\_i\\
  — one pair per subject per site.

## Usage

``` r
simulate_freesurfer_data2(
  path,
  data_structure = list(cohort1 = list(sessions = c("01", "02"), n_subjects = 100),
    cohort2 = list(sessions = c("01", "02"), n_subjects = 150)),
  pheno = NULL,
  betas = c(age = -0.1, sex = 0.05, time = -0.2),
  sigma2 = 0.25,
  subj_intercept_sd = 0.2,
  subj_slope_sd = 0.1,
  site_sd = 0.1,
  vw_mean = 2.5,
  vw_sd = 0.5,
  fs_template = "fsaverage",
  roi_subset = c("temporalpole", "frontalpole", "entorhinal"),
  location_association = NULL,
  measure = "thickness",
  hemi = "lh",
  fwhmc = "fwhm10",
  seed = 3108,
  verbose = TRUE
)
```

## Arguments

- path:

  Character string. Root directory for the dataset. Created recursively
  if it does not exist.

- data_structure:

  Named list defining cohorts/sites. Each element is a list with:

  `"sessions"`

  :   Character vector of session labels.

  `"n_subjects"`

  :   Integer number of subjects.

- pheno:

  `data.frame` as returned by
  [`simulate_long_pheno_data2`](https://seredef.github.io/verywise/reference/simulate_long_pheno_data2.md),
  containing at minimum the columns `site`, `id`, `time` (0-based
  integer), `folder_id`, plus any covariate columns named in `betas`.

- betas:

  Named numeric vector of fixed-effect coefficients. Names identify the
  covariates; currently supported names are `"age"`, `"sex"`, `"time"`,
  and `"wisdom"`. Set a coefficient to `0` to include the covariate
  without injecting a signal (useful for null-effect benchmarks).
  Example: `c(age = -0.1, sex = 0.05, time = -0.2)`.

- sigma2:

  Numeric scalar (\> 0). Within-subject residual variance
  \\\varepsilon\_{itkv}\\. Default: `0.25`.

- subj_intercept_sd:

  Numeric scalar (\\\geq 0\\). Standard deviation of the subject-level
  random intercept \\b\_{iv}\\. Default: `0.2`.

- subj_slope_sd:

  Numeric scalar (\\\geq 0\\). Standard deviation of the subject-level
  random slope for `time`. Set to `0` to remove random slopes entirely.
  Default: `0.1`.

- site_sd:

  Numeric scalar (\\\geq 0\\). Standard deviation of the site-level
  random intercept \\u\_{kv}\\. Default: `0.1`.

- vw_mean:

  Numeric scalar. Mean of the vertex-specific intercept \\\mu_v\\. For
  cortical thickness a realistic value is around `2.5` mm. Default:
  `2.5`.

- vw_sd:

  Numeric scalar. Standard deviation of \\\mu_v\\. Default: `0.5`.

- fs_template:

  Character string. FreeSurfer template; determines the number of
  vertices per surface file. Options:

  - `"fsaverage"` — 163 842 vertices (default)

  - `"fsaverage6"` — 40 962 vertices

  - `"fsaverage5"` — 10 242 vertices

  - `"fsaverage4"` — 2 562 vertices

  - `"fsaverage3"` — 642 vertices

- roi_subset:

  Character vector of ROI names that define the active vertices. The
  remaining vertices are set to `0` and excluded from analysis. Default:
  `c("temporalpole", "frontalpole", "entorhinal")`.

- location_association:

  Optional character vector. If supplied, the fixed-effect signal from
  `betas` is injected *only* within these ROIs; the remaining active
  vertices in `roi_subset` are generated under a null model (\\\beta =
  0\\). Useful for testing spatial localisation of discovered effects.

- measure:

  Character string. Surface measure label; used for file naming only.
  Default: `"thickness"`.

- hemi:

  Character string. Hemisphere: `"lh"` or `"rh"`. Default: `"lh"`.

- fwhmc:

  Character string. Smoothing label; used for file naming only. Default:
  `"fwhm10"`.

- seed:

  Integer. Random seed for reproducibility. Default: `3108`.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

## Value

Invisibly returns a named list of ground-truth simulation parameters:

- `beta0`:

  Length-\\V\_{\mathrm{roi}}\\ vector of vertex-specific intercepts
  \\\mu_v\\.

- `betas`:

  The `betas` argument as supplied.

- `sigma2`, `subj_intercept_sd`, `subj_slope_sd`, `site_sd`:

  Variance parameters as supplied.

- `u_site`:

  Length-\\K\\ vector of realised site random intercepts.

- `signal_vertices`:

  Logical vector (length = total vertices) indicating which vertices
  carry a non-zero fixed effect.

## See also

[`simulate_longit_dataset2`](https://seredef.github.io/verywise/reference/simulate_longit_dataset2.md)

## Author

Serena Defina, 2024.
