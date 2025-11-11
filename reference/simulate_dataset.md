# Simulate a longitudinal brain surface dataset with associated phenotype data

Generates a synthetic longitudinal dataset for multiple sites/cohorts,
each with multiple timepoints/sessions per subject. The function
produces:

- Brain surface data in FreeSurfer `.mgh` format, organised in a
  verywise folder structure (see vignettes for details).

- A matching `pheno` data frame with mock participant sex and age, saved
  as `"phenotype.csv"` in the `path` directory.

This is useful for testing pipelines or demonstrations where realistic
FreeSurfer-style data and phenotypic information are required.

## Usage

``` r
simulate_dataset(
  path,
  data_structure = list(cohort1 = list(sessions = c("01", "02", "03"), n_subjects = 10),
    cohort2 = list(sessions = c("01", "02"), n_subjects = 20)),
  fs_template = "fsaverage",
  roi_subset = c("temporalpole", "frontalpole", "entorhinal"),
  simulate_association = NULL,
  location_association = NULL,
  overwrite = TRUE,
  seed = 3108,
  verbose = TRUE,
  ...
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

- overwrite:

  Logical (default = `TRUE`). Whether to overwrite an existing phenotype
  file.

- seed:

  Integer (default = `3108`). Random seed.

- verbose:

  Logical (default = `TRUE`). If `TRUE`, print progress messages.

- ...:

  Additional arguments passed to
  [`simulate_freesurfer_data`](https://seredef.github.io/verywise/reference/simulate_freesurfer_data.md).

## Value

Invisibly returns `NULL`. Data and phenotype files are written to
`path`.

## See also

[`simulate_freesurfer_data`](https://seredef.github.io/verywise/reference/simulate_freesurfer_data.md),
[`simulate_long_pheno_data`](https://seredef.github.io/verywise/reference/simulate_long_pheno_data.md)

## Author

Serena Defina, 2024.
