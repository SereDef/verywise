# Simulate longitudinal FreeSurfer vertex-wise surface data

Writes one `.mgh` file per observation (subject × session) to `path`,
mimicking the FreeSurfer recon-all output layout. Vertex values are
generated from a mixed-effects model with:

- ROI-specific fixed effects specified via `roi_associations`;

- a scalar subject random intercept (shared across all vertices);

- an optional scalar site random intercept.

## Usage

``` r
simulate_freesurfer_data(
  path,
  pheno = NULL,
  data_structure = list(cohort1 = list(sessions = c("01", "02"), n_subjects = 100),
    cohort2 = list(sessions = c("01", "02"), n_subjects = 150)),
  roi_associations = list(),
  simulate_other_rois = FALSE,
  hemi = "lh",
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

- pheno:

  `data.frame` in long format from
  [`simulate_long_pheno_data()`](https://seredef.github.io/verywise/reference/simulate_long_pheno_data.md),
  or `NULL` to auto-generate a minimal bookkeeping frame from
  `data_structure`.

- data_structure:

  Named list; see
  [`simulate_long_pheno_data()`](https://seredef.github.io/verywise/reference/simulate_long_pheno_data.md).
  Ignored when `pheno` is supplied.

- roi_associations:

  Named list of named numeric vectors. Names of the list are
  Desikan-Killiany ROI labels; names of each vector are column names in
  `pheno`; values are beta coefficients. Example:
  `list(temporalpole = c(age = 0.3, wisdom = 0.5))`.

- simulate_other_rois:

  Logical; if `TRUE` all non-association ROIs are simulated under the
  null.

- hemi:

  One of `"lh"` or `"rh"`.

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

  Integer random seed.

- verbose:

  Logical; progress messages.

## Value

`NULL` invisibly. Side effect: `.mgh` files written to
`path/<folder_id>/surf/<hemi>.<measure>.<fwhmc>.<fs_template>.mgh`.

## ROI targeting

- `simulate_other_rois = FALSE`:

  Only ROIs in `roi_associations` are filled; all other vertices remain
  0.

- `simulate_other_rois = TRUE`:

  All Desikan-Killiany cortical ROIs are filled. ROIs absent from
  `roi_associations` are simulated under the null (no fixed effects,
  noise only).

## Data-generating process

For observation \\i\\ and vertex \\v\\ in ROI \\r\\: \$\$y\_{iv} =
\mathbf{x}\_i^\top \boldsymbol{\beta}\_r + b_i + b\_{\text{site}(i)} +
\varepsilon\_{iv}\$\$ where \\b_i \sim N(0, \sigma\_\text{subj}^2)\\,
\\b_s \sim N(0, \sigma\_\text{site}^2)\\, and \\\varepsilon\_{iv} \sim
N(\mu\_\text{vw}, \sigma\_\text{vw}^2)\\. Random effects are scalar (not
vertex-specific).

## See also

[`simulate_long_pheno_data()`](https://seredef.github.io/verywise/reference/simulate_long_pheno_data.md),
[`simulate_longit_dataset()`](https://seredef.github.io/verywise/reference/simulate_longit_dataset.md)

## Examples

``` r
if (FALSE) { # dir.exists("path/to/simulated_brains")
simulate_freesurfer_data(
  path = 'path/to/simulated_brains',
  pheno = pheno,
  roi_associations = list(temporalpole = c(age = 0.3)),
  hemi = "both"
)
}
```
