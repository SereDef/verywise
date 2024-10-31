
# **`verywise`** <a href="https://seredef.github.io/verywise/"><img src="man/figures/logo.png" align="right" height="139" alt="verywise website" /></a>

### vertex-wise, whole-brain linear mixed models

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of `verywise` is to offer a flexible, user-friendly interface
to whole-brain analysis of neuro-imaging data that has been
pre-processed using [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/).

The package was specifically designed for the analysis of *longitudinal*
(e.g. multi-session) and/or *multi-site* neuroimaging data.

Currently, `verywise` allows the estimation of vertex-wise **Linear
Mixed Models**, but will be extended to other statistical models in the
future.

It can handle imputed (phenotype) data from several packages (`mice`,
`mi`, `amelia`, etc.).

Multiple testing correction is currently achieved using MCZ simulations
from FreeSurfer. This means that you will need FreeSurfer installed and
correctly set up.

## Installation

You can install the development version of `verywise` from
[GitHub](https://github.com/) with:

``` r
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
  subj_dir = "./VeryWiseUserUser/SimulatedExample", # Neuro-imaging data location
  hemi = "both", # (default) which hemispheres to run
  pheno = phenotype, # An R object already in memory or path to file 
)
```

## Visualization

To quickly inspect and plot your results you can use our interactive
application, [BrainMapp](https://github.com/SereDef/brainmapp)
