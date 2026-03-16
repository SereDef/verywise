# library(lme4)
library(verywise)
library(profvis)

# Example data
set.seed(3108)

temp_test_dir <- system.file("tests/testthat/fixtures", package = "verywise")

subj_dir <- file.path(temp_test_dir, "fs7")
outp_dir <- file.path('~/Desktop/vw_tmp_results')

test_formula <- vw_area ~ sex + age + wisdom + (1 | id)

fs_home = "/Applications/freesurfer/7.4.1" # mac only

pheno = read.csv(file.path(subj_dir, "phenotype.csv"))
measure = 'area'
hemi = 'lh'
fs_template = "fsaverage"
fwhm = 10

ss1 <- build_supersubject(
      subj_dir = subj_dir,
      folder_ids = pheno[,'folder_id'],
      supsubj_dir = file.path(outp_dir, 'ss1'),
      measure = measure,
      hemi = hemi,
      fs_template = fs_template,
      n_cores = 1,
      fwhmc = paste0("fwhm", fwhm),
      save_rds = FALSE,
      error_cutoff = 20,
      verbose = TRUE
    )
# ss2 <- build_supersubject(
#       subj_dir = subj_dir,
#       folder_ids = pheno[,'folder_id'],
#       supsubj_dir = file.path(outp_dir, 'ss2'),
#       measure = measure,
#       hemi = hemi,
#       fs_template = fs_template,
#       n_cores = 4,
#       fwhmc = paste0("fwhm", fwhm),
#       save_rds = FALSE,
#       error_cutoff = 20,
#       verbose = TRUE
#     )

# dim(ss1) 
# dim(ss2)

# expect_equal(ss1[,2], ss2[, 2], tolerance = 1e-10)

# Test the main function runs ================================================
res = run_vw_lmm(
    formula = test_formula,
    pheno = read.csv(file.path(subj_dir, "phenotype.csv")),
    subj_dir = subj_dir,
    outp_dir = outp_dir,
    hemi = "lh",
    fs_template = "fsaverage",
    apply_cortical_mask = TRUE,
    folder_id = "folder_id",
    tolerate_surf_not_found = 20,
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
    verbose = TRUE)

# Profile it (sequential required) ===========================================
p <- profvis::profvis(
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
    verbose = TRUE)
)
htmlwidgets::saveWidget(p, file.path(outp_dir, "profvis.html"))
browseURL(file.path(outp_dir, "profvis.html"))
# Check tresults are as expected ---------------------------------------------------
for (term in c("wisdom")) { #, "age", "sexMale", "(Intercept)")) {

  p = plot_vw_map(res_dir = outp_dir,
                  term = term,
                  hemi = "lh",
                  measure = "area",
                  outline_rois = c('temporalpole', 'frontalpole','entorhinal'), 
                  fs_home = fs_home)
}

# ==================================================================================
library(lme4)
library(microbenchmark)

model_template <- precompile_model(formula = test_formula, 
  tmp_data = pheno, 
  tmp_y = ss1[, 1], measure = 'area',
  lmm_control = lme4::lmerControl(), 
  REML = TRUE, verbose = TRUE)

vertex <- ss1[, which(locate_roi('temporalpole'))[1] ]

n_obs <- nrow(pheno)
test_y <- rnorm(n_obs, mean = 6.5, sd = 0.2)

pheno['vw_area'] <- test_y

test_formula <- vw_area ~ sex + age + (1 | id)

test_template <- lmer(test_formula, data = pheno)

  # Use the template for single_lmm
result <- single_lmm(pheno, test_y, 'vw_area', model_template = test_template)


ps <- profvis::profvis({
single_lmm(pheno, y = vertex,
  y_name = "vw_area", model_template = model_template, weights = NULL)
})
htmlwidgets::saveWidget(ps, file.path(outp_dir, "profvis_single.html"))
browseURL(file.path(outp_dir, "profvis_single.html"))

pooled_stats <- vw_pool(list(outp, outp), m = 2, pvalue_method="t-as-z")

# mb <- microbenchmark(
#   se_from_summary = {
#     se1 <- as.matrix(summary(fit)$coefficients[, "Std. Error"])
#   },
#   se_from_vcov = {
#     se2 <- as.matrix(sqrt(diag(as.matrix(stats::vcov(fit)))))
#   },
#   times = 1000L,
#   check = 'identical'
# )

# print(mb)
