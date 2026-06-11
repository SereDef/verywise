# Simulate a complete longitudinal brain surface dataset

High-level wrapper that:

1.  Generates a longitudinal phenotype `data.frame` via
    [`simulate_long_pheno_data()`](https://seredef.github.io/verywise/reference/simulate_long_pheno_data.md);

2.  Writes it to `path/phenotype.csv`;

3.  Calls
    [`simulate_freesurfer_data()`](https://seredef.github.io/verywise/reference/simulate_freesurfer_data.md)
    for each requested hemisphere.

## Usage

``` r
simulate_longit_dataset(
  path,
  data_structure = list(cohort1 = list(sessions = c("01", "02", "03"), n_subjects = 100),
    cohort2 = list(sessions = c("01", "02"), n_subjects = 150)),
  baseline = list(age = c(mean = 10, sd = 0.5), sex = c(levels = c("Male", "Female")),
    wisdom = c(mean = 0, sd = 1)),
  change = list(age = c(mean = 4, sd = 0.5), wisdom = c(mean = 1, sd = 0.5)),
  roi_associations = list(temporalpole = c(age = 1.3, sex = 0.5), entorhinal = c(age =
    0.9), frontalpole = c(wisdom = 0.7)),
  simulate_other_rois = FALSE,
  hemi = "both",
  measure = "thickness",
  vw_mean = 2.5,
  vw_sd = 0.5,
  subj_sd = 0.2,
  site_sd = 0.1,
  fs_template = "fsaverage",
  fwhmc = "fwhm10",
  seed = 3108,
  verbose = TRUE
)
```

## Arguments

- path:

  Character string; root directory in which to write `.mgh` files.
  Created if it does not exist.

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

- roi_associations:

  Named list of named numeric vectors. Names of the list are
  Desikan-Killiany ROI labels; names of each vector are column names in
  `pheno`; values are beta coefficients. Example:
  `list(temporalpole = c(age = 0.3, wisdom = 0.5))`.

- simulate_other_rois:

  Logical; if `TRUE` all non-association ROIs are simulated under the
  null.

- hemi:

  One of `"lh"`, `"rh"`, or `"both"` (default: `"both"`).

- measure:

  Character; FreeSurfer surface measure, e.g. `"thickness"` or `"area"`.

- vw_mean:

  Numeric; mean of the vertex-level residual noise (i.e. grand-mean
  cortical thickness).

- vw_sd:

  Numeric \\\geq 0\\; SD of vertex-level residual noise.

- subj_sd:

  Numeric \\\geq 0\\; SD of the subject random intercept.

- site_sd:

  Numeric \\\geq 0\\; SD of the site random intercept. Ignored when
  there is only one site (intercept drawn as a single value).

- fs_template:

  Character; FreeSurfer template space, e.g. `"fsaverage"`.

- fwhmc:

  Character; smoothing kernel label, e.g. `"fwhm10"`.

- seed:

  Integer random seed for reproducibility.

- verbose:

  Logical

## Value

`NULL` invisibly. Side effects: `phenotype.csv` and vertex-wise `.mgh`
files written under `path`.

## See also

[`simulate_long_pheno_data()`](https://seredef.github.io/verywise/reference/simulate_long_pheno_data.md),
[`simulate_freesurfer_data()`](https://seredef.github.io/verywise/reference/simulate_freesurfer_data.md)

## Examples

``` r
if (FALSE) { # dir.exists("path/to/simulated/dataset")
simulate_longit_dataset(
  path = 'path/to/simulated/dataset/',
  data_structure = list(
    GENR = list(sessions = c("01", "02", "03"), n_subjects = 50)
  ),
  roi_associations = list(temporalpole = c(age = 1.3)),
  hemi = "lh"
)
}
```
