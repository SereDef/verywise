library(lme4)

# Example data
set.seed(3108)

subj_dir <- test_path("fixtures", "fs7")
pheno <- read.csv(file.path(subj_dir, "phenotype.csv"))

test_formula <- vw_area ~ sex + age + (1 | id)
fs_home = "/Applications/freesurfer/7.4.1" # mac only

# ==============================================================================
test_that("run_vw_lmm runs end-to-end with simulated data", {
  skip_on_os("windows")

  if (!dir.exists(fs_home)) {
    skip("FreeSurfer not found in FREESURFER_HOME")
  }

  # outp_dir <- withr::local_tempdir()
  # If you need to inspect results (note: running test() instead of check())
  outp_dir <- file.path(subj_dir, 'results')

  # Run function
  result <- run_vw_lmm(
    formula = test_formula,
    pheno = pheno,
    subj_dir = subj_dir,
    outp_dir = outp_dir,
    hemi = "lh",
    fs_template = "fsaverage",
    apply_cortical_mask = TRUE,
    folder_id = "folder_id",
    tolerate_surf_not_found = 20,
    use_model_template = TRUE,
    weights = NULL,
    # prioritize speed over accuracy
    lmm_control = lme4::lmerControl(calc.derivs = FALSE,
                                    use.last.params = TRUE,
                                    check.rankX = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nlev = "ignore",
                                    check.nlev.gtreq.5 = "ignore",
                                    check.nlev.gtr.1 = "ignore",
                                    check.nobs.vs.nRE = "ignore",
                                    check.formula.LHS = "ignore",
                                    check.scaleX = "ignore",
                                    check.conv.grad = "ignore",
                                    check.conv.singular = "ignore",
                                    check.conv.hess = "ignore"),
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

  # Structure tests
  expect_type(result, "list")
  expect_named(result, c("coef", "se", "p", "resid"))
  expect_s4_class(result$coef, "FBM")
  expect_s4_class(result$se, "FBM")
  expect_s4_class(result$p, "FBM")
  expect_s4_class(result$resid, "FBM")

  expect_true(file.exists(
    file.path(outp_dir, 'lh.area.stack1.coef.mgh')))
  expect_true(file.exists(
    file.path(outp_dir, 'lh.area.stack2.p.mgh')))
  expect_true(file.exists(
    file.path(outp_dir, 'lh.area.stack3.cache.th30.abs.sig.ocn.mgh')))

})
