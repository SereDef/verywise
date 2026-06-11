# Simulate longitudinal phenotype data

Generates a long-format phenotype `data.frame` for one or more
cohorts/sites across multiple sessions / time points.

## Usage

``` r
simulate_long_pheno_data(
  data_structure = list(cohort1 = list(sessions = c("01", "02"), n_subjects = 100),
    cohort2 = list(sessions = c("01", "02"), n_subjects = 150)),
  baseline = list(age = c(mean = 10, sd = 0.5), sex = c(levels = c("Male", "Female")),
    wisdom = c(mean = 0, sd = 1)),
  change = list(age = c(mean = 4, sd = 0.5), wisdom = c(mean = 1, sd = 0.5)),
  seed = 3108,
  verbose = TRUE
)
```

## Arguments

- data_structure:

  Named list specifying site/cohort structure. Each element must be a
  list with:

  `sessions`

  :   Character vector of unique session labels.

  `n_subjects`

  :   Positive integer, number of subjects.

- baseline:

  Named list defining baseline distributions. Supports continuous
  covariates (`c(mean, sd)`) and categorical covariates
  (`c(levels = ...)`). See section **Supported covariates**.

- change:

  Named list of `c(mean, sd)` specifying the per-wave mean shift and
  noise SD for longitudinal covariates. Only continuous variables should
  appear here.

- seed:

  Integer random seed for reproducibility.

- verbose:

  Logical

## Value

A `data.frame` in long format with one row per subject × session. Always
contains columns `site`, `id`, `time`, `folder_id`, + all covariates
declared in `baseline`.

## Supported covariates

Covariates are declared through the `baseline` argument. The type of
covariate is inferred from the names of each list element:

- Continuous (e.g. `age`, `wisdom`):

  Specify `c(mean = m, sd = s)`. Baseline values are drawn from
  \\N(\var{m}, \var{s}^2)\\.

- Categorical (e.g. `sex`):

  Specify `c(levels = c("Male", "Female"))`. Each subject is assigned a
  level with equal probability and this value is held constant across
  sessions.

## Change model

Follow-up values for continuous variables listed in `change` are
computed as: \$\$y\_{i,s} = y\_{i,1} + (s-1)\bar{\delta} + \varepsilon,
\quad \varepsilon \sim N(0, \sigma\_\delta^2)\$\$ where \\s\\ is the
1-based session index. The deviation from baseline is drawn
independently at each session. The variance \\\sigma\_\delta^2\\ is
constant across sessions; only the mean shift accumulates linearly.

## See also

[`simulate_freesurfer_data()`](https://seredef.github.io/verywise/reference/simulate_freesurfer_data.md),
[`simulate_longit_dataset()`](https://seredef.github.io/verywise/reference/simulate_longit_dataset.md)

## Examples

``` r
pheno <- simulate_long_pheno_data(
  data_structure = list(
    GENR = list(sessions = c("01", "02"), n_subjects = 50)
  ),
  baseline = list(age = c(mean = 10, sd = 1), sex = c(levels = c("Male", "Female"))),
  change   = list(age = c(mean = 4, sd = 0.5))
)
#> ⠙ Generate phenotype data
#> ℹ 100 total observations (from 1 site and 2 waves)
#> ⠙ Generate phenotype data
#>   Variables: id, site, time, age, sex, and folder_id
#> ⠙ Generate phenotype data
#> ✔ Generate phenotype data [44ms]
#> 
head(pheno)
#>   id site time       age    sex         folder_id
#> 1  1 GENR   01 11.565062 Female GENR/sub-1_ses-01
#> 2  2 GENR   01  7.527207   Male GENR/sub-2_ses-01
#> 3  3 GENR   01 10.194473 Female GENR/sub-3_ses-01
#> 4  4 GENR   01 11.811407   Male GENR/sub-4_ses-01
#> 5  5 GENR   01 10.131829   Male GENR/sub-5_ses-01
#> 6  6 GENR   01  9.036134 Female GENR/sub-6_ses-01
```
