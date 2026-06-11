# Build a minimal phenotype data frame from a data structure specification

Creates a long-format data frame with one row per subject × session
combination, containing bookkeeping columns used by the simulation
pipeline. No phenotypic variables are added.

## Usage

``` r
build_minimal_pheno(data_structure)
```

## Arguments

- data_structure:

  A named list where each element corresponds to one site/cohort. Every
  element must contain:

  `sessions`

  :   Character vector of unique session labels (e.g. `c("01", "02")`).

  `n_subjects`

  :   Positive integer giving the number of subjects at this site.

## Value

A `data.frame` in long format with columns:

- `id`:

  Integer subject identifier (within-site).

- `site`:

  Site/cohort label (character).

- `time`:

  Session label (character).

- `folder_id`:

  FreeSurfer-style path fragment, e.g. `"cohort1/sub-1_ses-01"`.
