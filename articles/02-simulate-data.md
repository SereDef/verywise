# Optional: simulating data

## Data poor? Simulate yourself some data, son

If you are eager to get started with `verywise` but do not have access
to a dataset (yet), you can generate both a set of ready-made FreeSurfer
**brain surface files** and a **phenotype file** to go with them.

The simulated dataset will automatically be stored into a `verywise`
directory structure and will be ready for you to analyse.

Let’s see how to get there. First, load the package:

``` r
library(verywise)
#> Welcome, verywise user!
#> This is version: 1.2.4
#> For questions, issues, and bug reports, please see https://github.com/SereDef/verywise
```

Now, decide what data you want to generate. In this (small) example, I
will simulate a dataset with 250 subjects, who belong to two cohorts
(including 100 and 150 subjects respectively) and all underwent two MRI
sessions (`"01"` and `"02"`). I am only interested in surface area, so I
will only be genetating `"area"` maps. The resolution of these fake
surface maps is 163842 vertices (corresponding to the most detailed
FreeSurfer template `"fsaverage"`).

To make this more realistic, I want to inject a positive association
between `"age"` and surface area in my data (beta = 3.5). This
association will only be present in the `"frontalpole"` region from the
Desikan-Killiany atlas.

``` r

data_structure <- list("cohort1" = list("n_subjects" = 100, 
                                        "sessions" = c("01","02")),
                       "cohort2" = list("n_subjects" = 150, 
                                        "sessions" = c("01","02")))

# Simulate FreeSurfer and phenotype dataset for both hemispheres
for (hemi in c('lh','rh')) {
  simulate_dataset(path = "path/to/verywise/simulated_example",
                   data_structure = data_structure,
                   measure = "area",
                   hemi = hemi,
                   fs_template = "fsaverage", # FreeSurfer template (largest available: 163842 vertices)
                   vw_mean = 5, # mean surface area value
                   vw_sd = 0.1, # standard deviation of surface area values
                   subj_sd = 0.02, # variability between subjects
                   site_sd = 0.01, # variability between sites
                   # (optionally) only simulate data within certain ROIs:
                   # roi_subset = c('temporalpole', 'frontalpole', 'entorhinal'), 
                   simulate_association = "3.5 * age",
                   location_association = "frontalpole")
}
```

This should give you everything you need to play around with `verywise`
model fitting. If you do run a model using this data now, you should be
able to recover an association with age located in a single cluster,
perfectly overlapping with the frontal pole region. How *neat*.

Under the hood,
[`simulate_dataset()`](https://seredef.github.io/verywise/reference/simulate_dataset.md)
will call two functions:

- [`simulate_freesurfer_data()`](https://seredef.github.io/verywise/reference/simulate_freesurfer_data.md)
  to generate the brain surface files, and
- [`simulate_long_pheno_data()`](https://seredef.github.io/verywise/reference/simulate_long_pheno_data.md)
  to generate the phenotype data

You can also call these two functions separately, in case you only need
one of the two data sources.

## Only generate the brain data

If you have a phenotype dataset already, but you are waiting for the
FreeSurfer data to cook.

In this case I have a sigle dataset / cohort called `"my_site"` with 100
subjects who underwent three MRI sessions (`"baseline"`, `"F1"` and
`"F2"`). The phenotype data is ready and stored in a CSV file located at
`"path/to/my_phenotype_data.rds"`. This dataset should contain 300 rows
(100 subjects x 3 sessions), and a few variables of interest, including
`"wisdom"`.

Let’s simulate left hemisphere thickness maps, smoothed at `10mm` FWHM,
with an overall mean of `6.5` (mm thickness) and a standard deviation of
`1.5`. This time, to save us some time, I will use a smaller surface
template (`"fsaverage3"`), which has only `40962` vertices in total.

I also want to inject a positive association between `"wisdom"` and
thickness (beta = 10) in the `"entorhinal"` and `"temporalpole"`
regions.

``` r

phenotype_data <- readRDS("path/to/my_phenotype_data.rds")

data_structure <- list("my_site" = list("sessions" = c("baseline", "F1", "F2"),
                                        "n_subjects" = 100))

# Simulate FreeSurfer dataset
simulate_freesurfer_data(path = "path/to/simulated_FreeSurfer_output",
                         data_structure = data_structure,
                         fs_template = "fsaverage3",
                         measure = "thickness",
                         hemi = "lh",
                         fwhmc = "fwhm10", # smoothness level
                         vw_mean = 6.5,
                         vw_sd = 1.5,
                         subj_sd = 0.2, 
                         site_sd = 0.1,
                         simulate_association = 10 * phenotype$wisdom,
                         location_association = c('temporalpole', 'entorhinal'),
                         seed = my_random_seed_I_did_not_forget_to_set)
```

## Only generate the phenotype data

On the other hand, if you already have some surfaces to use but you
would like a phenotype to go with it, you can generate a minimal “long
format” dataframe (with variables: `folder_id`, `id`, `site`, `age`,
`sex`, `wisdom`).

``` r
# Simulate phenotype dataset, using the same data_structure as above
phenotype_data <- simulate_long_pheno_data(data_structure = data_structure,
                                           seed = my_random_seed_I_did_not_forget_to_set)
```

## 

Next article: [Run a vertex-wise linear mixed
model](https://seredef.github.io/verywise/articles/03-run-vw-lmm.md)
