# Simulate longitudinal phenotype data

Generates synthetic phenotype data in long format for multiple cohorts
and multiple sessions per subject. The set of columns produced is driven
by `betas`: only the covariates whose names appear in `betas` are
generated, plus mandatory bookkeeping columns.

Covariate generation rules:

- `age` — Baseline age \\\sim \mathcal{N}(\code{mean\\age},\\
  \code{sd\\age}^2)\\, z-scored within site. At each follow-up session
  the value is incremented by `interval_years` plus \\\mathcal{N}(0,\\
  0.5^2)\\ jitter.

- `sex` — Binary 0 / 1 \\\sim \mathrm{Bernoulli}(0.5)\\; constant across
  sessions.

- `time` — 0-based integer session index (`0, 1, 2, ...`); enters the
  model as a within-person linear time trend.

- `wisdom` — Continuous nuisance covariate \\\sim \mathcal{N}(10,\\
  3^2)\\, z-scored within site; constant across sessions.

## Usage

``` r
simulate_long_pheno_data2(
  data_structure = list(cohort1 = list(sessions = c("01", "02"), n_subjects = 100),
    cohort2 = list(sessions = c("01", "02"), n_subjects = 150)),
  betas = c(age = -0.1, sex = 0.05, time = -0.2),
  mean_age = 20,
  sd_age = 2,
  interval_years = 2,
  seed = 3108,
  verbose = TRUE
)
```

## Arguments

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

- mean_age:

  Numeric scalar. Population mean of baseline age. Default: `20`.

- sd_age:

  Numeric scalar. Population standard deviation of baseline age.
  Default: `2`.

- interval_years:

  Numeric scalar. Mean number of years between consecutive sessions.
  Default: `2`.

- seed:

  Integer. Random seed for reproducibility. Default: `3108`.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

## Value

A `data.frame` in long format with mandatory columns `site`, `id`,
`session` (original character label), `time` (0-based integer index),
`folder_id`, followed by any covariate columns implied by `betas`.

## See also

[`simulate_longit_dataset2`](https://seredef.github.io/verywise/reference/simulate_longit_dataset2.md),
[`simulate_freesurfer_data2`](https://seredef.github.io/verywise/reference/simulate_freesurfer_data2.md)

## Author

Serena Defina, 2024.
