# Simulate (longitudinal) phenotype data

Generates synthetic phenotype data in long format for multiple
cohorts/sites and multiple timepoints/sessions per subject. Each record
contains:

- Subject ID

- Session/timepoint

- Sex

- Age

- Wisdom

- `folder_id` field matching the FreeSurfer directory structure

## Usage

``` r
simulate_long_pheno_data(
  data_structure = list(cohort1 = list(sessions = c("01", "02"), n_subjects = 100),
    cohort2 = list(sessions = c("01", "02"), n_subjects = 150)),
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

- seed:

  Integer (default = `3108`). Random seed.

- verbose:

  Logical (default = `TRUE`). If `TRUE`, print progress messages.

## Value

A `data.frame` in long format with columns: `site`, `id`, `time`, `sex`,
`age`, `wisdom`, `folder_id`.

## See also

[`simulate_dataset`](https://seredef.github.io/verywise/reference/simulate_dataset.md),
[`simulate_freesurfer_data`](https://seredef.github.io/verywise/reference/simulate_freesurfer_data.md)

## Author

Serena Defina, 2024.
