# Getting ready, step 1: Installation and system requirements

## Installation

You can install the development version of `verywise` from
[GitHub](https://github.com/).

*Note*

: you may need to install the [`devtools`
package](https://devtools.r-lib.org/) first (if it is not already
installed).

``` r
# Install devtools (if not there yet)
if (!requireNamespace("devtools"))
  install.packages('devtools')

# Install verywise
devtools::install_github("SereDef/verywise")
```

This will also download some other R packages that are needed, so give
it a little minute.

Finally, you can attach the package to make it easier to use, like so:

``` r
library(verywise)
#> Welcome, verywise user!
#> This is version: 1.2.2
#> For questions, issues, and bug reports, please see https://github.com/SereDef/verywise
```

## System and software requirements

In order for `verywise` to run smoothly, you will need:

- **`R`** (version \>= **4.1**) installed (*Dah!*).
- A **Unix system**: `verywise` works on Mac and on most Linux
  distributions, but **not on Windows**. This is because we currently
  rely on functions from the FreeSurfer suite which require Unix/Mac.
- [**FreeSurfer**](https://surfer.nmr.mgh.harvard.edu/) installed.

### Installing and setting up FreeSurfer

TODO
