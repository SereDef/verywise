# Run many \`verywise\` analyses in parallel (on SLURM)

The `run_vw_*()"` functions from `verywise` are designed to run a single
vertex-wise analysis: i.e. a single model and in only one hemisphere.
This is done on purpose, to keep things simple and modular, and maximise
computational efficiency, while avoiding issues with nested
parallelization.

However, in practice, you will probably want to run more that a single
hemisphere vertex-wise analysis, for example, you likely want to analyse
both left and right hemispheres, and possibly assess multiple models
(i.e. looking a more than one brain outcome, or at more than one set of
predictors and so on).

You can do this *sequentially* of course, just by calling the
`run_vw_*()"` multiple times with varying specifications. But if you
have a large number of analyses to run, this can quickly become
time-consuming.

So here are some tips on running multiple `verywise` analyses *in
parallel*, using **SLURM job arrays** (since this is a common setup on
HPC clusters).

### 1. Create a list of analysis specifications

First, let’s design our `analysis.R` script, which will be called by
each job in the SLURM array. This will change quite a bit depending on
your specific analyses, but here is an example of what it could look
like.

We start by defining some constants: things that will be the same across
all the analysises (note: you don’t really *need* to set these up as
variables, but I think it makes the code cleaner and easier to maintain,
and I am in charge here so we do this how I like it).

``` r
# Load verywise
library(verywise)

# Define project paths and constants ----------------------------------------------

proj_dir <- "path/to/main/project/directory"

fs_home <- "/path/to/FREESURFER_HOME"
subj_dir <- "path/to/freesurfer/data"

outp_dir <- file.path(proj_dir, "results")

pheno_filepath <- file.path(proj_dir, 'phenotype_clean.rds')

# this part of the formula will be constant across all analyses
covariates <- 'birth_weight + SES + ethnicity' 
random_effects <- '(1 | id)'

# Get SLURM parameters ------------------------------------------------------------

analysis_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID', unset = 1))
n_cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK', unset = 1))
```

Note the `analysis_id` and `n_cores` variables, which catch the array
task ID and number of CPUs from the SLURM environment. These are used to
select *which of our many analyses* to run, and how many resources can
be used to run it. We will define these values in the SLURM job script
below.

Now, let’s organise all our analyses into a ***grid of analysis
specifications*** (i.e. combinations of hemispheres, outcomes, models,
datasets etc). Each row of this grid will correspond to a single
analysis.

In this example I want to analyse both hemispheres (`lh` and `rh`),
three different outcomes (`thickness`, `area` and `w_g.pct`), and a two
model specifications: one including an interaction between sex and age
(`*`), and one without this interaction (`+`).

This gives us a total of *12 analyses* (2 hemispheres x 3 outcomes x 2
models).

The base R function
[`expand.grid()`](https://rdrr.io/r/base/expand.grid.html) is a useful
way to create such a grid of all possible parameter combinations, but
you can also create this manually if you prefer (for example, if only
certain combinations make sense for you).

We then use the `analysis_id` variable to select one specific parameter
combination (i.e. a row in the grid) and use that to define the current
analysis. For example, we: a. take the value of hemisphere (just making
sure it is a character string and not a factor); b. format the outcome
variable name correctly (adding the `vw_` prefix); c. define the model
formula for the current analysis (e.g. with or without interaction
term).

``` r
# Define parameter space ==========================================================

hemis <- c('lh','rh')
outcs <- c('thickness', 'area', 'w_g.pct')
model <- c('*', '+')

param_grid <- expand.grid('hemi'=hemis, 'outc'=outcs, 'mode'=model)

# Select and process current parameter combination --------------------------------

params <- param_grid[analysis_id, ]

hemi <- as.character(params$hemi)

outc <- paste0('vw_', params$outc)

model_spec <- as.formula(paste(outc, '~ age', params$mode, 'sex +', covariates, '+', random_effects))
```

Finally, let’s run the analysis with all our “current” parameters:

``` r
# Run analysis --------------------------------------------------------------------

output <- run_vw_lmm(formula = model_spec,
                     subj_dir = subj_dir,
                     outp_dir = outp_dir,
                     hemi = hemi,
                     pheno = pheno_filepath, 
                     n_cores = n_cores,
                     FS_HOME = fs_home,
                     save_ss = TRUE)
```

### 2. Create a job script

All of this was fun, but we still need to create a SLURM job script
(e.g. `run_analyses.sh`), which will submit the job array to the
scheduler and stuff done for us. Here we will specify how many analyses
we want to run (i.e. the size of the array), and how many resources to
allocate to each job.

``` sh
#!/bin/bash

#SBATCH --job-name=verywise_job_array
#SBATCH --array=1-12 # 12 models in total (see above)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64   # NOTE R limit parallel processes = 124
#SBATCH --time=1-00:00:00    # 1 day time limit, just to be on the safe side
#SBATCH --error=analysis_log%a
#SBATCH --output=analysis_log%a
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@email.com

# Load modules (not always but sometimes necessary on some HPCs) -------------
module load R

# Run analyses (in parallel) -------------------------------------------------
Rscript analysis.R
```

Then, inside your project directory, submit the job array to SLURM with:
`sbatch run_analyses.sh`.

### Performance tip 1: avoid implicit parallelization

On some systems, *implicit parallelism* in low-level matrix algebra
libraries (BLAS/LAPACK) can interfere with explicit parallelization used
internally by `verywise`, secretly slowing your analyses down. If you
feel like processing is taking too long for your liking, I recommend
disabling these implicit threading libraries (*before starting R*).

You can do so, for example by setting the following environment
variables in your shell or in your SLURM job script (before calling
`Rscript [...]`):

``` sh
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
```

Also note that using a very large number of cores (e.g. **\> 120**) may
sometimes cause worker initialization or other issues (e.g. R parallel
processes limits).

### Performance tip 2: precompute and re-use the super-subject matrix

If you find yourself running multiple analyses on the same dataset
(i.e. *same dataset* and *same brain outcome*), you can save some
computation time by pre-computing the “super-subject matrix” only once,
and re-using it across all your analyses.

``` r
library(verywise)

supsubj_dir = "/path/to/save/ss"

ss <- build_supersubject(
                  subj_dir = "/path/to/freesurfer/subjects",
                  folder_ids = pheno[, 'folder_id'], # assuming 'pheno' is your phenotype data frame
                  supsubj_dir = supsubj_dir,
                  measure = "thickness",
                  hemi = "lh",
                  fs_template = "fsaverage",
                  n_cores = 4)
```

Then simply use `supsubj_dir` as your `subj_dir` argument in
[`run_vw_lmm()`](https://seredef.github.io/verywise/reference/run_vw_lmm.md),
and we will do the rest :)

## 

Next article: [Inspect and *visualize* verywise
results](https://seredef.github.io/verywise/articles/05-visualize-results.md)
