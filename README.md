
<!-- README.md is generated from README.Rmd. Please edit that file -->

# verywise

<!-- badges: start -->
<!-- badges: end -->

The goal of verywise is to offer a flexible, user-friedly interface to
whole-brain (hemispehere) analysis of FreeSurfer MGH format data. The
package allows application of several statistical models per vertex. It
can handle imputed data from specific packages (mice, mi, amelia, etc.).
Multiple testing correction is achieved using MCZ simulations from
FreeSurfer.

## Installation

You can install the development version of verywise from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("SereDef/verywise")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(verywise)
## basic example code
```
