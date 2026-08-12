# Optional: simulating data

## Data poor? Simulate yourself some data, son

If you are eager to get started with `verywise` but do not have access
to a dataset (yet), you can generate both a set of **brain surface
files** that mimic real FreeSurfer output and a **phenotype file** to go
with them.

The simulated dataset will automatically be stored into a `verywise`
directory structure and will be ready for you to analyze.

Let’s see how to get there. First, load the package:

``` r

library(verywise)
#> Loaded verywise 1.3.8
```

Now, decide what data you want to generate. In this (small) example, I
will simulate a dataset with 250 subjects, who belong to two cohorts
(including 100 and 150 subjects respectively) and all underwent two MRI
sessions (`"01"` and `"02"`).

These subjects are on average 10 years old at baseline (i.e. at session
`"01"`) and then attendend session `"02"` after 4 years (on average). I
also simulate their sex (50% males) and wisdom (mean = 0 and SD = 1),
which also increases by 1 SD on average over time (i.e., between
sessions).

On the brain side, i simulate surface `"area"` maps across both
hemispheres a. The resolution of these fake surface maps is 163842
vertices (corresponding to the most detailed FreeSurfer template
`"fsaverage"`). The mean area is 2.5 (SD = 0.5)

I also inject a few associations between my covariates and surface area
in my favorite regions from the Desikan-Killiany atlas.

``` r


data_structure <- list(
      cohort1 = list(sessions = c("01", "02"), n_subjects = 100),
      cohort2 = list(sessions = c("01", "02"), n_subjects = 150))

baseline <- list(
  age = c(mean=10, sd=0.5), 
  sex = c(levels=c('Male','Female')),
  wisdom = c(mean = 0, sd = 1))

change <- list(
  age = c(mean=4, sd=0.5),
  wisdom = c(mean = 1, sd = 0.5))

roi_associations < list(
  temporalpole = c(age = 1.3, sex = 0.5), 
  entorhinal = c(age = 0.9), 
  frontalpole = c(wisdom = 0.7))


# Simulate FreeSurfer and phenotype dataset for both hemispheres
simulate_longit_dataset(
    path = "path/to/verywise/simulated_example",
    data_structure = data_structure,
    baseline = baseline,
    change = change,
    roi_associations = roi_associations,
    measure = "area", vw_mean = 2.5, vw_sd = 0.5,
    fs_template = "fsaverage")
```

This should give you everything you need to play around with `verywise`
model fitting. If you do run a model using this data now, you should be
able to recover an association with age located in two clusters
(perfectly overlapping with the temporalpole and entorhinal regions),
same goes for sex and age. How *neat*.

Under the hood,
[`simulate_longit_dataset()`](https://seredef.github.io/verywise/reference/simulate_longit_dataset.md)
will call two functions:

- [`simulate_long_pheno_data()`](https://seredef.github.io/verywise/reference/simulate_long_pheno_data.md)
  to generate the phenotype data
- [`simulate_freesurfer_data()`](https://seredef.github.io/verywise/reference/simulate_freesurfer_data.md)
  to generate the brain surface files

You can also call these two functions separately, in case you only need
one of the two data sources.

### Only generate the brain data

If you have a phenotype dataset already, but you are waiting for the
FreeSurfer data to cook.

In this case the phenotype data is ready and stored in a CSV file
located at `"path/to/my_phenotype_data.rds"`. This dataset should be in
long format and it should contain, a “folder_id” variable as well as the
the two variables “good_variable” and “bad_variable” that I am going to
use to simulate a positive and a negative association respectively.

The following code will simulate left hemisphere thickness maps,
smoothed at `10mm` FWHM, with an overall mean of `6.5` (mm thickness)
and a standard deviation of `1.5`. This time, to save us some time, I
will use a smaller surface template (`"fsaverage3"`), which has only
`40962` vertices in total.

``` r


phenotype_data <- readRDS("path/to/my_phenotype_data.rds")

# Simulate FreeSurfer dataset
simulate_freesurfer_data(path = "path/to/simulated_FreeSurfer_output",
                         pheno = phenotype_data,
                         roi_associations =  list(
                         frontalpole = c(good_variable = 0.7, bad_variable = -0.3)),
                         measure = "thickness",
                         hemi = "lh",
                         fs_template = "fsaverage3", # smaller template
                         vw_mean = 6.5,
                         vw_sd = 1.5,
                         subj_sd = 0.2, 
                         site_sd = 0.1)
```

### Only generate the phenotype data

On the other hand, if you already have some surfaces to use but you
would like a phenotype to go with it, you can generate a minimal “long
format” data.frame with the variables you like (though you may need to
adapt the “folder_id” variable to your brain data structure).

``` r

# Simulate phenotype dataset, using the same data_structure as above
phenotype_data <- simulate_long_pheno_data(data_structure = data_structure, 
                                           baseline = baseline, 
                                           change = change,
                                           seed = my_random_seed_I_did_not_forget_to_set)
```

### Other parameters

[`simulate_longit_dataset()`](https://seredef.github.io/verywise/reference/simulate_longit_dataset.md)
also takes a handful of other (optional) parameters that you may find
useful. Please refer to the function documentation for more details.

## Simulate “distributed” data

`verywise` also handles datasets that are not accessible by a single
analyst (distributed framework). To simulate a similar dataset for
testing:

``` r

site_sizes = c(
   site1 = 50, 
   site2 = 100, 
   site3 = 5
)

for (hemi in c('rh','lh')) {

  true_estimates <- simulate_distrib_dataset(
    path = "path/to/testing/folder",
    site_sizes = site_sizes,
    fs_template = 'fsaverage',
    measure = 'area',
    hemi = hemi,
    fwhmc = 'fwhm10',
    vw_mean = 5,
    vw_sd = 0.1,
    tau2 = 0.1, # Between-site variance of the random intercept
    sigma2 = 0.05, # Within-site residual variance
    betas = c(sex = -0.5, age = 0.8),
    roi_subset = c('temporalpole', 'frontalpole', 'entorhinal'),
    location_association = 'frontalpole',
    overwrite = TRUE,
    seed = 42,
    verbose = TRUE)
  
}
```

## 

Next article: [Run a vertex-wise linear mixed
model](https://seredef.github.io/verywise/articles/03-run-vw-lmm.md)
