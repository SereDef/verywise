# **`verywise`**

### vertex-wise, whole-brain linear mixed models

The goal of `verywise` is to offer a flexible, user-friendly interface
to whole-brain analysis of neuro-imaging data that has been
pre-processed using [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/).

The package was specifically designed for the analysis of *longitudinal*
(e.g. multi-session) and/or *multi-site* neuroimaging data.

Currently, `verywise` allows the estimation of vertex-wise **Linear
Mixed Models**, and **meta-analysis**, but will be extended to other
statistical models in the future.

It can handle imputed (phenotype) data from several packages (`mice`,
`mi`, `amelia`, etc.).

Multiple testing correction is currently achieved using MCZ simulations
from FreeSurfer. This means that you will need FreeSurfer installed and
correctly set up.

You can find more info and extended **tutorials**
[here](https://seredef.github.io/verywise/index.html).

## Installation

You can install the development version of `verywise` from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("SereDef/verywise")

# or 
# install.packages("devtools")
devtools::install_github("SereDef/verywise")
```

## Basic use

``` r
library(verywise)
```

``` r
# Run a linear mixed model
run_vw_lmm(
  formula = vw_thickness ~ sex * age + site + (1 | id), # model formula
  pheno = long_format_data, # An R object already in memory, or "path/to/phentype/data"
  subj_dir = "/path/to/freesurfer/subjects", # Neuro-imaging data location
  outp_dir = "/path/to/output", # Where you want to store results
  hemi = "lh", # (default) or "rh": which hemisphere to run
  n_cores = 4  # parallel processing
)
```

``` r
# Run a meta-analysis
run_vw_meta(
  term = "age", # Which "term" / predictor / effect to pool
  hemi = "lh", # (default) or "rh": which hemisphere to run
  measure = "area", # (default) or any available FreeSurfer metric.
  res_dirs = c("/path/to/study1/results", "/path/to/study2/results"),
  study_names = c("Study1", "Study2"),
  n_cores = 4  # parallel processing
)
```

## Visualization

To inspect and plot your results, you can use our interactive web
application,
[verywiseWIZard](https://github.com/SereDef/verywise-wizard). You can
run this locally or try it out
[here](https://seredef-verywise-wizard.share.connect.posit.cloud/).

Plots can also be generated using `verywise` like so:

``` r
# Plot result brain map (requires FreeSurfer for templates and reticulate for interface with Python-based plotting libraries)
plot_vw_map(
  term = "age",
  hemi = "lh",
  measure = "area",
  res_dir = "/path/to/output",
  outline_rois = c("entorhinal", "precuneus"),
  fs_home = "/path/to/FREESURFER_HOME"
)
```

## Note

This is a spin-off of the [`QDECR`](https://www.qdecr.com/) package,
which handles linear regression models.

## Funders

![Funders](reference/figures/funders.png)

This work was supported by the *FLAG-ERA* grant
[**Infant2Adult**](https://www.infant2adult.com/home) and by The
Netherlands Organization for Health Research and Development (ZonMw,
grant number 16080606).

All feedback and contributions are **very welcome!**.
