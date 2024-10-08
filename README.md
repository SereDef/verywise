
<!-- README.md is generated from README.Rmd. Please edit that file -->

# verywise: vertex-wise analysis of neuroimaging data

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of **`verywise`** is to offer a flexible, user-friendly
interface to whole-brain (hemisphere) analysis of FreeSurfer MGH format
data.

The package allows application of several statistical models per vertex,
including linear mixed models for the analysis of longitudinal and
multi-site neuroimaging data.

It can handle imputed data from several packages (`mice`, `mi`,
`amelia`, etc.).

Multiple testing correction is currently achieved using MCZ simulations
from FreeSurfer.

## Installation

You can install the development version of `verywise` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("SereDef/verywise")
```

``` r
# Load it up
library(verywise)
```

## Preparing your dataset

You need FreeSurfer installed and a phenotype file that follows verywise
structure.

## Data poor? Simulate yourself some data, son

``` r
# Simulate FreeSurfer and phenotype dataset
simulate_data(subj_dir = "./VeryWiseUserUser/SimulatedExample")
```

## Basic use

``` r
# Run linear mixed model
run_vw_lmm(vw_thickness ~ sex * age + site + (1 | id), # model formula
  subj_dir = "./VeryWiseUserUser/SimulatedExample", # Data location
  hemi = "lh", # left hemisphere
  n_cores = 1 # Number of cores for parallel processing
)
```
