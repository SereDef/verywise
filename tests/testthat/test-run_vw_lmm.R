library(lme4)

# Example data
set.seed(3108)

subj_dir <- test_path("fixtures", "fs7")
pheno <- read.csv(file.path(subj_dir, "phenotype.csv"))

test_formula <- vw_area ~ sex + age + wisdom + (1 | id)

quick_run <- lme4::lmerControl(calc.derivs = FALSE, use.last.params = TRUE,
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
  check.conv.hess = "ignore")

# ==============================================================================
test_that("run_vw_lmm runs end-to-end with simulated data", {

  # skip_on_os("windows") # FreeSurfer not supported
  # skip_on_os("linux") # System settings implicit parallelism
  fs_home <- skip_if_no_freesurfer()

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
    weights = NULL,
    # prioritize speed over accuracy
    lmm_control = quick_run,
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
  expect_named(result, c("coef", "se", "p", "mfit", "resid", "clust"))
  expect_s4_class(result$coef, "FBM")
  expect_s4_class(result$se, "FBM")
  expect_s4_class(result$p, "FBM")
  expect_s4_class(result$mfit, "FBM")
  expect_s4_class(result$resid, "FBM")
  expect_s4_class(result$clust, "FBM")

  expect_true(file.exists(
    file.path(outp_dir, 'lh.area.stack1.coef.mgh')))
  expect_true(file.exists(
    file.path(outp_dir, 'lh.area.stack2.p.mgh')))
  expect_true(file.exists(
    file.path(outp_dir, 'lh.area.stack3.cache.th30.abs.sig.ocn.mgh')))

  # Residuals should NOT be persisted to disk when save_residuals = FALSE
  expect_false(file.exists(result$resid$rds))

  mgh_resid_files <- list.files(outp_dir, pattern = "residuals\\.mgh$",
                                 recursive = TRUE, full.names = TRUE)
  expect_length(mgh_resid_files, 0)

})

# ── save_cov: term covariance extraction through run_vw_lmm() ────────────────

test_that("run_vw_lmm rejects save_cov with non-existent or wrong number of terms", {

  fs_home <- skip_if_no_freesurfer()

  outp_dir <- withr::local_tempdir()

  expect_error(
    run_vw_lmm(
      formula = test_formula,
      pheno = pheno,
      subj_dir = subj_dir,
      outp_dir = outp_dir,
      hemi = "lh",
      save_cov = "age", # only one term, should error before any modeling starts
      verbose = FALSE
    ),
    regexp = "Only one covariance can be extracted"
  )

  expect_error(
    run_vw_lmm(
      formula = test_formula,
      pheno = pheno,
      subj_dir = subj_dir,
      outp_dir = outp_dir,
      hemi = "lh",
      save_cov = c("age", "not_a_real_term"),
      verbose = FALSE
    ),
    regexp = "not present in the model"
  )
})

test_that("run_vw_lmm returns and saves extracted term covariance when save_cov is set", {

  fs_home <- skip_if_no_freesurfer()

  outp_dir <- withr::local_tempdir()

  result <- run_vw_lmm(
    formula = test_formula,
    pheno = pheno,
    subj_dir = subj_dir,
    outp_dir = outp_dir,
    hemi = "lh",
    fs_template = "fsaverage",
    save_cov = c("age", "wisdom"),
    lmm_control = quick_run,
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

  expect_true("cov" %in% names(result))
  if ("cov" %in% names(result)) {
    expect_s4_class(result$cov, "FBM")
  }

  cov_mgh_files <- list.files(outp_dir, pattern = "\\.cov\\.mgh$",
                               recursive = TRUE, full.names = TRUE)
  expect_true(length(cov_mgh_files) >= 1)
})

# ── save_residuals: residual matrix persistence through run_vw_lmm() ────────

test_that("run_vw_lmm persists residuals matrix to disk when save_residuals = TRUE", {

  fs_home <- skip_if_no_freesurfer()

  outp_dir <- withr::local_tempdir()

  result <- run_vw_lmm(
    formula = test_formula,
    pheno = pheno,
    subj_dir = subj_dir,
    outp_dir = outp_dir,
    hemi = "lh",
    fs_template = "fsaverage",
    lmm_control = quick_run,
    seed = 42,
    n_cores = 1,
    chunk_size = 1000,
    FS_HOME = fs_home,
    fwhm = 10,
    mcz_thr = 30,
    cwp_thr = 0.025,
    save_optional_cluster_info = FALSE,
    save_ss = FALSE,
    save_residuals = TRUE,
    verbose = FALSE
  )

  # The FBM backing file (.rds) should now exist on disk
  expect_true(file.exists(result$resid$rds))
  expect_true(file.exists(result$resid$backingfile))

  # residuals.mgh should still not be written
  mgh_resid_files <- list.files(outp_dir, pattern = "residuals\\.mgh$",
                                 recursive = TRUE, full.names = TRUE)
  expect_length(mgh_resid_files, 0)

  # The saved FBM should be re-loadable and match the in-memory matrix
  reloaded <- bigstatsr::big_attach(result$resid$rds)
  expect_equal(reloaded[], result$resid[])
})