# Create a temp directory for shared (simualted) test data
temp_test_dir <- system.file("tests/testthat/fixtures", package = "verywise")
message('Setting up test data at', temp_test_dir)
# temp_test_dir <- withr::local_tempdir()
# temp_test_dir <- file.path(tempdir(), "verywise-testdata")
# dir.create(temp_test_dir, showWarnings = FALSE, recursive = TRUE)

# Setup
measure = "area"
fwhmc = "fwhm10"
n_subs = 3

verbosity = TRUE

# ------------------------------------------------------------------------------
# 1. test all functionality except the FreeSurfer-related ones using fsaverage3
#    (642 vertices)
fs3_data_path <- file.path(temp_test_dir, "fs3")

data_structure = list("site1" = list("sessions" = c("01",'02'),
                                     "n_subjects" = n_subs),
                      "site2" = list("sessions" = c("01",'02'),
                                     "n_subjects" = n_subs))

for (hemi in c('lh','rh')) {
  verywise::simulate_dataset(path = fs3_data_path,
                   data_structure = data_structure,
                   fs_template = "fsaverage3",
                   hemi = hemi,
                   measure = measure,
                   fwhmc = fwhmc,
                   vw_mean = 6.5,
                   vw_sd = 0.5,
                   simulate_association = NULL,
                   verbose = verbosity)
}

pheno <- readr::read_csv(file.path(fs3_data_path, "phenotype.csv"))

# Save as rds as well
saveRDS(pheno, file.path(fs3_data_path, "phenotype.rds"))

# Inject 25% missingness at random
pheno_miss <- pheno

set.seed(3108)
pheno_miss[sample(nrow(pheno), round(nrow(pheno)*0.25)), 'sex'] <- NA
pheno_miss[sample(nrow(pheno), round(nrow(pheno)*0.25)), 'age'] <- NA

pheno_imp <- mice::mice(pheno_miss, m = 3, maxit = 3)

saveRDS(pheno_imp, file.path(fs3_data_path, "phenotype_mids.rds"))

# Super-subject
for (hemi in c('lh','rh')) {
  ss <- verywise::build_supersubject(
    subj_dir = fs3_data_path,
    folder_ids = pheno$folder_id,
    supsubj_dir = file.path(fs3_data_path, "ss"),
    measure = measure,
    hemi = hemi,
    n_cores = 1,
    fwhmc = fwhmc,
    fs_template = "fsaverage3",
    save_rds = TRUE,
    verbose = verbosity)
}


# ------------------------------------------------------------------------------
# 2. test FreeSurfer-related functionality using the "full" fsaverage template
#    (163842 vertices)
fs7_data_path <- file.path(temp_test_dir, "fs7")

n_subs <- 10

data_structure = list("site1" = list("sessions" = c("01",'02'),
                                     "n_subjects" = n_subs))

for (hemi in c('lh','rh')) {
  verywise::simulate_dataset(path = fs7_data_path,
                   data_structure = data_structure,
                   fs_template = "fsaverage",
                   measure = measure,
                   hemi = hemi,
                   fwhmc = fwhmc,
                   vw_mean = 0.6,
                   vw_sd = 0.1,
                   simulate_association = "0.5 * age",
                   verbose = verbosity)
}

pheno <- readr::read_csv(file.path(fs7_data_path, "phenotype.csv"))

# Super-subject
for (hemi in c('lh','rh')) {
  ss <- verywise::build_supersubject(
    subj_dir = fs7_data_path,
    folder_ids = pheno$folder_id,
    supsubj_dir = file.path(fs7_data_path,"ss"),
    measure = measure,
    hemi = hemi,
    n_cores = 1,
    fwhmc = fwhmc,
    fs_template = "fsaverage",
    save_rds = TRUE,
    verbose = verbosity)
}
