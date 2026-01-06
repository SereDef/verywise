# Running vertex-wise linear mixed models

> **Warning**  
> `verywise` lets you run vertex-wise analyses however you want to. For
> example, you can add interactions and splines / polynomials to your
> models, as many random intercept and slopes as your heart desires…,
> you can use liberal thresholds for multiple test correction and test
> as many hypothesis as you have time for. Remember: this software is
> just a tool, and it is up to you to use it sensibly.

With your set-up and data ready, it is finally time to use some
`verywise` functionality.

First, load the package in our R sesssion, so we call its functions
directly.

``` r
library(verywise)
```

## Run a vertex-wise linear mixed model

The main function in the `verywise` package is `run_vw_lmm` which fits a
vertex-wise linear mixed model to the data. Here is an example:

``` r
# Run a linear mixed model
run_vw_lmm(
  formula = vw_thickness ~ sex * age + site + (1 | id), # model formula
  pheno = phenotype, # An R object already in memory or a path to file
  subj_dir = "./VeryWiseUser/SimulatedExample", # Neuro-imaging data location
  hemi = "both", # (default) which hemispheres to run
  n_cores = 4 # number of cores for parallel processing
  ...
)
```

The function will first check the inputs and then run
[`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) on every vertex
in both hemispheres.

The most important parameters you need to know about:

- **`formula`**: the model formula specifying the linear mixed model.
  This uses `lme4` syntax.
- **`pheno`** : either the phenotype data object (already loaded in the
  global environment) or a string containing a file path. Supported file
  extensions are: rds, csv, txt and sav.
- **`subj_dir`** : a string containing the path to the FreeSurfer data,
  this expects a verywise structure.

You can optionally also specify:

- **`outp_dir`**: a string containing the path, where do you want
  results to be stored. If none is provided, a “results” sub-directory
  will be created inside `subj_dir`.
- **`hemi`**: which hemispheres to run (default = “both”)
- **`seed`**: (default = 3108) random seed.
- **`n_cores`**: (default = 1) number of cores for parallel processing.
- **`FS_HOME`**: FreeSurfer directory, i.e. `$FREESURFER_HOME`.

Check out the output directory.
