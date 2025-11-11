# Simulate longitudinal FreeSurfer vertex-wise data

Simulates FreeSurfer-formatted brain surface data for multiple
cohorts/sites and multiple timepoints/sessions. The output folder
structure emulates that produced by the FreeSurfer `recon-all` command,
with directories named: `sub-<ID>_ses-<SESSION>`.

The vertex-wise data are stored as `.mgh` files, with optional
simulation of an association with a phenotype variable.

## Usage

``` r
simulate_freesurfer_data(
  path,
  data_structure = list(cohort1 = list(sessions = c("01", "02"), n_subjects = 100),
    cohort2 = list(sessions = c("01", "02"), n_subjects = 150)),
  fs_template = "fsaverage",
  measure = "thickness",
  hemi = "lh",
  fwhmc = "fwhm10",
  vw_mean = 6.5,
  vw_sd = 1.5,
  subj_sd = 0.2,
  site_sd = 0.1,
  roi_subset = c("temporalpole", "frontalpole", "entorhinal"),
  simulate_association = NULL,
  location_association = NULL,
  seed = 3108,
  verbose = TRUE
)
```

## Arguments

- path:

  Character string. Directory where the dataset should be created. Will
  be created if it does not exist.

- data_structure:

  Named list defining cohorts/sites. Each element is a list with:

  `"sessions"`

  :   Character vector of session labels.

  `"n_subjects"`

  :   Integer number of subjects.

- fs_template:

  Character (default = `"fsaverage"`). FreeSurfer template for vertex
  registration. This is used to determine the size of the synthetic
  brain surface data. Options:

  - `"fsaverage"` = 163842 vertices (highest resolution)

  - `"fsaverage6"` = 40962 vertices

  - `"fsaverage5"` = 10242 vertices

  - `"fsaverage4"` = 2562 vertices

  - `"fsaverage3"` = 642 vertices

- measure:

  Character (default = `"thickness"`). Surface measure to simulate. This
  is used for file names.

- hemi:

  Character (default = `"lh"`). Hemisphere code: `"lh"` (left) or `"rh"`
  (right). This is used for file names.

- fwhmc:

  Character (default = `"fwhm10"`). Full-width half maximum smoothing
  parameter. This is used for file names.

- vw_mean:

  Numeric (default = `6.5`). Mean of the simulated vertex-wise values.

- vw_sd:

  Numeric (default = `0.5`). Standard deviation of the simulated
  vertex-wise values.

- subj_sd:

  Numeric (default = `0.2`). Standard deviation of the random intercept
  for subject (relevant in multi-session datasets).

- site_sd:

  Numeric (default = `0.1`). Standard deviation of the random intercept
  for site (relevant in multi-site datasets).

- roi_subset:

  Character vector (default = c('temporalpole', 'frontalpole',
  'entorhinal')). Vertex-wise data is simulated by default only within a
  smaller subset (~1.5%) of the total surface. The rest of the vertex
  values are set to 0, so they won't be analysed, saving time during
  estimation. The region locations are extracted from the annotation
  files in that are distributed with FreeSurfer and saved internally in
  R/sysdata.rda.

- simulate_association:

  Optional. If numeric, must be of length equal to the number of
  generated files; if character, must have the format
  `"<beta> * <variable_name>"`. Associations are injected into one small
  region (the entorhinal cortex).

- location_association:

  Optional string or character vector. If specified, the association is
  only present within these ROIs. The rest of the vertex values will be
  set to have no relationship with any of the predictors. The region
  locations are extracted from the annotation files in that are
  distributed with FreeSurfer and saved internally in R/sysdata.rda.

- seed:

  Integer (default = `3108`). Random seed.

- verbose:

  Logical (default = `TRUE`). If `TRUE`, print progress messages.

## Value

Invisibly returns `NULL`. Vertex-wise data are written to `path` in
FreeSurfer-compatible format.

## See also

[`simulate_dataset`](https://seredef.github.io/verywise/reference/simulate_dataset.md),
[`simulate_long_pheno_data`](https://seredef.github.io/verywise/reference/simulate_long_pheno_data.md)

## Author

Serena Defina, 2024.
