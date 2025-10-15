# library(lme4)
library(profvis)

# Example data
set.seed(3108)

temp_test_dir <- system.file("tests/testthat/fixtures", package = "verywise")

subj_dir <- file.path(temp_test_dir, "fs7")
outp_dir <- file.path(subj_dir, 'results')

test_formula <- vw_area ~ sex + age + wisdom + (1 | id)

fs_home = "/Applications/freesurfer/7.4.1" # mac only

run_vw_lmm(
    formula = test_formula,
    pheno = read.csv(file.path(subj_dir, "phenotype.csv")),
    subj_dir = subj_dir,
    outp_dir = outp_dir,
    hemi = "lh",
    fs_template = "fsaverage",
    apply_cortical_mask = TRUE,
    folder_id = "folder_id",
    tolerate_surf_not_found = 20,
    use_model_template = TRUE,
    weights = NULL,
    seed = 42,
    n_cores = 1,
    chunk_size = 1000,
    FS_HOME = fs_home,
    fwhm = 10,
    mcz_thr = 30,
    cwp_thr = 0.025,
    save_optional_cluster_info = FALSE,
    save_ss = FALSE,
    save_residuals = FALSE,
    verbose = TRUE
)

# NOTE: find profiler that works well with parallel settings 
# profvis::profvis({
#   run_vw_lmm(
#       formula = test_formula,
#       pheno = pheno,
#       subj_dir = subj_dir,
#       outp_dir = outp_dir,
#       hemi = "lh",
#       fs_template = "fsaverage",
#       apply_cortical_mask = TRUE,
#       folder_id = "folder_id",
#       tolerate_surf_not_found = 20,
#       use_model_template = TRUE,
#       weights = NULL,
#       seed = 42,
#       n_cores = 1,
#       chunk_size = 1000,
#       FS_HOME = fs_home,
#       fwhm = 10,
#       mcz_thr = 30,
#       cwp_thr = 0.025,
#       save_optional_cluster_info = FALSE,
#       save_ss = FALSE,
#       save_residuals = FALSE,
#       verbose = TRUE
#   )
# })
