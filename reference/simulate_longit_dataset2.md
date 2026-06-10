# Simulate a longitudinal brain surface dataset with associated phenotype data

Generates a synthetic longitudinal dataset for multiple sites/cohorts,
each with multiple timepoints/sessions per subject. The generative model
for vertex \\v\\, subject \\i\\ nested in site \\k\\, at session \\t\\
is:

\$\$ y\_{itkv} = \mu_v + \sum_j \beta_j \cdot x\_{ijtk} + b\_{iv} +
b^{\mathrm{time}}\_{iv} \cdot t + u\_{kv} + \varepsilon\_{itkv} \$\$

where \\\mu_v \sim \mathcal{N}(\code{vw\\mean},\\ \code{vw\\sd}^2)\\ is
a vertex-specific intercept drawn once and shared across all subjects
and sites, \\\beta_j\\ are fixed effects supplied via `betas`, \\b\_{iv}
\sim \mathcal{N}(0,\\ \code{subj\\intercept\\sd}^2)\\ is a subject-level
random intercept, \\b^{\mathrm{time}}\_{iv} \sim \mathcal{N}(0,\\
\code{subj\\slope\\sd}^2)\\ is a subject-level random slope for time
(set `subj\_slope\_sd = 0` to suppress it), \\u\_{kv} \sim
\mathcal{N}(0,\\ \code{site\\sd}^2)\\ is a site-level random intercept,
and \\\varepsilon\_{itkv} \sim \mathcal{N}(0,\\ \code{sigma2})\\ is
residual noise.

The function writes:

- Brain surface data in FreeSurfer `.mgh` format organised in a verywise
  folder structure (see vignettes for details).

- A phenotype `data.frame` saved as `phenotype.csv` in `path`.

## Usage

``` r
simulate_longit_dataset2(
  path,
  data_structure = list(cohort1 = list(sessions = c("01", "02", "03"), n_subjects = 10),
    cohort2 = list(sessions = c("01", "02"), n_subjects = 20)),
  betas = c(age = -0.1, sex = 0.05, time = -0.2),
  sigma2 = 0.25,
  subj_intercept_sd = 0.2,
  subj_slope_sd = 0.1,
  site_sd = 0.1,
  vw_mean = 2.5,
  vw_sd = 0.5,
  mean_age = 20,
  sd_age = 2,
  interval_years = 2,
  fs_template = "fsaverage",
  roi_subset = c("temporalpole", "frontalpole", "entorhinal"),
  location_association = NULL,
  measure = "thickness",
  hemi = "lh",
  fwhmc = "fwhm10",
  overwrite = TRUE,
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

- mean_age:

  Numeric scalar. Population mean of baseline age. Default: `20`.

- sd_age:

  Numeric scalar. Population standard deviation of baseline age.
  Default: `2`.

- interval_years:

  Numeric scalar. Mean number of years between consecutive sessions.
  Default: `2`.

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

- overwrite:

  Logical. If `FALSE` and `phenotype.csv` already exists in `path`, it
  is reused and no new brain files are written. Default: `TRUE`.

- seed:

  Integer. Random seed for reproducibility. Default: `3108`.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

## Value

Invisibly returns a named list of ground-truth parameters:

- `beta0`:

  Length-\\V\_{\mathrm{roi}}\\ vector of simulated vertex-specific
  intercepts \\\mu_v\\.

- `betas`:

  The `betas` argument as supplied.

- `sigma2`:

  Scalar residual variance as supplied.

- `subj_intercept_sd`, `subj_slope_sd`, `site_sd`:

  Random-effect standard deviations as supplied.

- `u_site`:

  Length-\\K\\ vector of realised site random intercepts.

- `signal_vertices`:

  Logical vector of length \\n\_{\mathrm{verts}}\\ indicating which
  global vertex indices carry a non-zero fixed effect.

Data and phenotype files are written to `path` as a side-effect.

## See also

[`simulate_long_pheno_data2`](https://seredef.github.io/verywise/reference/simulate_long_pheno_data2.md),
[`simulate_freesurfer_data2`](https://seredef.github.io/verywise/reference/simulate_freesurfer_data2.md),
[`simulate_distrib_dataset`](https://seredef.github.io/verywise/reference/simulate_distrib_dataset.md)

## Author

Serena Defina, 2024.

## Examples

``` r
truth <- simulate_longit_dataset2(
  path              = tempfile("longit_sim"),
  betas             = c(age = -0.1, sex = 0.05, time = -0.2),
  sigma2            = 0.25,
  subj_intercept_sd = 0.3,
  subj_slope_sd     = 0.05,
  site_sd           = 0.1,
  fs_template       = "fsaverage3"
)
#> ── Simulating longitudinal dataset ────────────────────────── verywise v1.3.6 ──
#> ! `path` specified (/tmp/RtmpeqW9nk/longit_sim1d49264765d) does not exist.
#>   I'll try to create it.
#> • creating phenotype file...
#> * creating FreeSurfer dataset (
#> left
#> hemisphere)...
#> ℹ Selected 11 vertices (1.7%) from 3 regions: entorhinal, temporalpole, and
#>   frontalpole
#> *
#> cohort1
#> :
#> 10
#> (subjects) x
#> 3
#> (sessions) cortical
#> thickness
#> files
#> *
#> cohort2
#> :
#> 20
#> (subjects) x
#> 2
#> (sessions) cortical
#> thickness
#> files
#> 
#> ── Done! :) ────────────────────────────────────────────────────────────────────
truth$betas
#>   age   sex  time 
#> -0.10  0.05 -0.20 
sum(truth$signal_vertices)
#> [1] 11
```
