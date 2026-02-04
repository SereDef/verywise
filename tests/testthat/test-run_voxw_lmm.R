set.seed(1)

# Minimal phenotype: 4 observations, 2 subjects, 2 timepoints
pheno <- data.frame(
  folder_id      = rep(paste0("sub", 1:2), each = 2),
  participant_id = rep(1:2, each = 2),
  time           = rep(1:2, times = 2),
  age            = rnorm(4, 10, 1),
  sex            = factor(rep(c("M", "F"), each = 2)),
  obs_id         = 1:4)

# Minimal super‑subject file: 4 observations x 3 voxels
ss_mat <- matrix(rnorm(4 * 3), nrow = 4, ncol = 3)

test_that("run_voxw_lmm runs on minimal data", {

  ss_file <- tempfile(fileext = ".csv")
  write.csv(ss_mat, ss_file, row.names = FALSE)

  out_dir <- tempfile("voxel_results")
  dir.create(out_dir)

  # Disable verbose to keep test output clean and avoid parallel log noise
  res <- run_voxw_lmm(
    formula    = vw_value ~ age + sex + (1 | participant_id),
    pheno      = pheno,
    ss_file    = ss_file,
    brain_template = 3,
    outp_dir   = out_dir,
    n_cores    = 1, 
    chunk_size = 2,
    verbose    = FALSE
  )

  # Structure of result
  expect_type(res, "list")
  expect_named(res, c("coef", "se", "p", "resid"))

  # Check that each element is an FBM and has the expected dimensions
  coef_fbm <- res$coef
  p_fbm    <- res$p
  r_fbm    <- res$resid

  expect_s4_class(coef_fbm, "FBM")
  expect_s4_class(p_fbm,    "FBM")
  expect_s4_class(r_fbm,    "FBM")
  
  expect_equal(dim(p_fbm),   dim(coef_fbm))
  expect_equal(dim(coef_fbm)[2], ncol(ss_mat))  # columns = voxels

  expect_equal(dim(r_fbm)[2], ncol(ss_mat))  # residuals: rows = obs, cols = voxels
  expect_equal(dim(r_fbm)[1], nrow(pheno))

  # Backing files are still present for results [TMP]
  expect_true(file.exists(coef_fbm$backingfile))

  # But super‑subject backing file "ss.bk" should be removed
  expect_false(file.exists(file.path(out_dir, "ss.bk")))
})

test_that("run_voxw_lmm errors on invalid inputs", {

  ss_file <- tempfile(fileext = ".csv")
  write.csv(matrix(rnorm(4 * 3), 4, 3), ss_file, row.names = FALSE)

  out_dir <- tempfile("bad_voxel_results")
  dir.create(out_dir)

  bad_pheno <- pheno
  bad_pheno$age <- NULL

  # Missing phenotype variable (age)
  expect_error(
    run_voxw_lmm(
      formula  = vw_value ~ age + sex + (1 | participant_id),
      pheno    = bad_pheno,
      ss_file  = ss_file,
      brain_template = 3,
      outp_dir = out_dir,
      n_cores  = 1,
      verbose  = FALSE
    ),
    regexp = "age", fixed = FALSE
  )

  # Invalid chunk_size
  expect_error(
    run_voxw_lmm(
      formula    = voxel_value ~ age + sex + (1 | participant_id),,
      pheno      = pheno,
      ss_file    = ss_file,
      outp_dir   = out_dir,
      chunk_size = 0,
      n_cores    = 1,
      verbose    = FALSE
    )
  )
})