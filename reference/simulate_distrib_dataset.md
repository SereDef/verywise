# Simulate a distributed (multi-site) cross-sectional brain surface dataset

Generates a synthetic distributed dataset where each site's data lives
in its own subdirectory, mirroring a real multi-site setup where no site
can access another's raw data. For each site \\k\\ and each vertex
\\v\\, the generative model is:

\$\$ y\_{ikv} = \mu_v + \sum\_{j} \beta_j \cdot x\_{ijk} + u\_{kv} +
\varepsilon\_{ikv} \$\$

where \\\mu_v\\ is a vertex-specific intercept drawn once from
\\\mathcal{N}(\code{vw\\mean},\\ \code{vw\\sd}^2)\\, \\\beta_j\\ are
fixed effects shared across all vertices and sites (supplied via
`betas`), \\u\_{kv} \sim \mathcal{N}(0, \tau^2_v)\\ is a site-level
random intercept, and \\\varepsilon\_{ikv} \sim \mathcal{N}(0,
\sigma^2_v)\\ is residual noise.

The output folder layout is:

    <path>/
      <site1>/
        phenotype.csv
        sub-001/surf/<hemi>.<measure>.<fwhmc>.<fs_template>.mgh
        sub-002/surf/...
      <site2>/
        phenotype.csv
        ...

Each site folder can therefore be passed directly to
[`run_vw_fed_local`](https://seredef.github.io/verywise/reference/run_vw_fed_local.md)
as the `subj_dir` and the matching `phenotype.csv` as `pheno`.

## Usage

``` r
simulate_distrib_dataset(
  path,
  site_sizes = c(site1 = 80L, site2 = 120L, site3 = 60L),
  betas = c(sex = -0.2, age = 0.3),
  tau2 = 0.01,
  sigma2 = 0.01,
  fs_template = "fsaverage",
  roi_subset = c("temporalpole", "frontalpole", "entorhinal"),
  location_association = NULL,
  measure = "thickness",
  hemi = "lh",
  fwhmc = "fwhm10",
  vw_mean = 5,
  vw_sd = 0.5,
  overwrite = TRUE,
  seed = 3108,
  verbose = TRUE
)
```

## Arguments

- path:

  Character string. Root directory for all site sub-folders. Created
  recursively if it does not exist.

- site_sizes:

  Named integer vector. Names become site/folder names; values give the
  number of subjects at each site. Default:
  `c(site1 = 80, site2 = 120, site3 = 60)`.

- betas:

  Named numeric vector of fixed-effect coefficients. Each name must
  correspond to a covariate that will be simulated in the phenotype
  data. Currently supported covariate names:

  `age`

  :   Continuous, drawn from \\\mathcal{N}(30, 5^2)\\ and then z-scored
      within each site.

  `sex`

  :   Binary (0/1), drawn from \\\text{Bernoulli}(0.5)\\.

  Example: `c(age = 0.3, sex = -0.2)`. Set any coefficient to `0` to
  include the covariate in the design matrix without injecting a signal
  (useful for null-effect benchmarks).

- tau2:

  Numeric scalar or length-\\V\_{\text{roi}}\\ vector. **Between-site
  variance** of the random intercept \\u\_{kv}\\. A value of `0`
  collapses the model to a fixed-effects OLS. Default: `0.5`.

- sigma2:

  Numeric scalar or length-\\V\_{\text{roi}}\\ vector. **Within-site**
  residual variance \\\varepsilon\_{ikv}\\. Default: `1`. `tau2` and
  `sigma2` together define the intra-class correlation \\\text{ICC} =
  \tau^2 / (\tau^2 + \sigma^2)\\.

- fs_template:

  Character string. FreeSurfer template; determines the total number of
  vertices in each `.mgh` file. Options:

  - `"fsaverage"` = 163842 vertices

  - `"fsaverage6"` = 40962 vertices

  - `"fsaverage5"` = 10242 vertices

  - `"fsaverage4"` = 2562 vertices

  - `"fsaverage3"` = 642 vertices

  Default: `"fsaverage"`.

- roi_subset:

  Character vector of ROI names used to restrict the active vertices
  (the rest are set to `0` and excluded from analysis). Vertex locations
  are read from the FreeSurfer annotation files stored in
  `R/sysdata.rda`. Default:
  `c("temporalpole", "frontalpole", "entorhinal")`.

- location_association:

  Character vector (optional). If supplied, the signal encoded in
  `betas` is injected *only* within these ROIs; the remaining active
  vertices (`roi_subset`) are generated under a null model (\\\beta =
  0\\). Useful for testing spatial localization of discovered effects.

- measure:

  Character string. Surface measure; used for file naming only. Default:
  `"thickness"`.

- hemi:

  Character string. Hemisphere: `"lh"` or `"rh"`. Default: `"lh"`.

- fwhmc:

  Character string. Smoothing label; used for file naming only. Default:
  `"fwhm10"`.

- vw_mean:

  Numeric scalar. Mean of the vertex-specific intercept \\\mu_v \sim
  \mathcal{N}(\code{vw\\mean},\\ \code{vw\\sd}^2)\\. For cortical
  thickness a realistic value is around `2.5` mm. Default: `2.5`.

- vw_sd:

  Numeric scalar. Standard deviation of the vertex-specific intercept
  distribution. Default: `0.5`.

- overwrite:

  Logical. If `FALSE` and `phenotype.csv` already exists in a site
  folder, the existing file is re-used and no new brain surface files
  are written for that site. Default: `TRUE`.

- seed:

  Integer. Random seed for reproducibility. Default: `3108`.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

## Value

Invisibly returns a list with the ground-truth parameters used during
data generation:

- `beta0`:

  Length-\\V\_{\text{roi}}\\ vector of simulated vertex-specific
  intercepts \\\mu_v\\.

- `betas`:

  The `betas` argument as supplied.

- `tau2`:

  Length-\\V\_{\text{roi}}\\ vector of between-site variances (after
  recycling).

- `sigma2`:

  Length-\\V\_{\text{roi}}\\ vector of residual variances (after
  recycling).

- `u`:

  Numeric matrix \\K \times V\_{\text{roi}}\\ of realised site random
  intercepts.

- `icc`:

  Length-\\V\_{\text{roi}}\\ vector of theoretical ICCs: \\\tau^2_v /
  (\tau^2_v + \sigma^2_v)\\.

- `signal_vertices`:

  Logical vector of length \\n\_{\text{verts}}\\ indicating which global
  vertex indices received a non-zero fixed effect (i.e. the intersection
  of `roi_subset` and `location_association`, if supplied).

Data and phenotype files are written to `path` as a side-effect.

## See also

[`run_vw_fed_local`](https://seredef.github.io/verywise/reference/run_vw_fed_local.md)
for the analysis function this feeds into,
[`simulate_longit_dataset`](https://seredef.github.io/verywise/reference/simulate_longit_dataset.md)
for the longitudinal (multi-session) variant.

## Author

Serena Defina, 2026.

## Examples

``` r
truth <- simulate_distrib_dataset(
  path = tempfile("fed_simulation"),
  site_sizes = c(site1 = 80, site2 = 120, site3 = 60),
  betas = c(age = 0.3, sex = -0.2),
  tau2 = 0.05,
  sigma2 = 0.1,
  fs_template = "fsaverage3"  # 642 vertices
)
#>  * the `path` specified ('/tmp/Rtmp8AvUa5/fed_simulation1e5a430645d') does not exist.
#>    I'll try to create it.
#>  * selected 11 vertices (1.7%)
#>  * the `site_dir` specified ('/tmp/Rtmp8AvUa5/fed_simulation1e5a430645d/site1') does not exist.
#>    I'll try to create it.
#>  * [site1] writing 80 surface files...
#>  * the `site_dir` specified ('/tmp/Rtmp8AvUa5/fed_simulation1e5a430645d/site2') does not exist.
#>    I'll try to create it.
#>  * [site2] writing 120 surface files...
#>  * the `site_dir` specified ('/tmp/Rtmp8AvUa5/fed_simulation1e5a430645d/site3') does not exist.
#>    I'll try to create it.
#>  * [site3] writing 60 surface files...
#> Done.

# Ground-truth ICC summary across active vertices
summary(truth$icc)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.3333  0.3333  0.3333  0.3333  0.3333  0.3333 

# Recovered site random intercepts (n_sites x n_vertices matrix)
dim(truth$u)
#> [1]  3 11
```
