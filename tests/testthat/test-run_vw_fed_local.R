# Example data
set.seed(3108)

# ── Shared test fixture ────────────────────────────────────────────────────────
subj_dir <- test_path('fixtures', 'fed')
# pheno <- read.csv(file.path(subj_dir, 'phenotype.csv'))

test_formula <- vw_area ~ sex + age
# fs_home = '/Applications/freesurfer/7.4.1' # mac only
outp_dir <- file.path(subj_dir, 'results')

# ── 1. Input validation ────────────────────────────────────────────────────────

test_that('missing formula or unknown outcome throws an error', {
  expect_error(
    run_vw_fed_local(
      site_name = 'site1',
      formula = 'some ~ silly + formula',
      pheno = file.path(subj_dir, 'site1', 'phenotype.csv'),
      subj_dir = file.path(subj_dir, 'site1'),
      outp_dir = outp_dir),
    regexp = 'must be a valid R formula object'
  )

  expect_error(
    run_vw_fed_local(
      site_name = 'site1',
      formula = vw_pokemon ~ sex + age,
      pheno = file.path(subj_dir, 'site1', 'phenotype.csv'),
      subj_dir = file.path(subj_dir, 'site1'), 
      outp_dir = outp_dir),
    regexp = 'should be a brain surface metric'
  )
})

test_that('missing subj_dir throws an error', {
  expect_error(
    run_vw_fed_local(
      site_name = 'site1',
      formula = test_formula,
      pheno = file.path(subj_dir, 'site1', 'phenotype.csv'),
      subj_dir = file.path(subj_dir, 'nonexistent_dir'),
      outp_dir = outp_dir),
    regexp = 'does not exist'  
  )
})

test_that('pheno missing folder_id column throws an error', {
  pheno <- utils::read.csv(file.path(subj_dir, 'site1', 'phenotype.csv'))
  pheno$folder_id <- NULL

  expect_error(
    run_vw_fed_local(
        site_name = 'site1',
        formula = test_formula,
        pheno = pheno,
        subj_dir = file.path(subj_dir, 'site1'),
        outp_dir = outp_dir),
    regexp = 'not found in data'  
  )
})

# ── 2. Successful execution ───────────────────────────────────────────────────

test_that('run_vw_fed_local completes without error', {

  site <- 'site1'

  local_test <- run_vw_fed_local(
      site_name = site,
      formula = test_formula,
      pheno = file.path(subj_dir, site, 'phenotype.csv'),
      subj_dir = file.path(subj_dir, site),
      outp_dir = outp_dir
    )
  
  expect_type(local_test, 'list')
  expect_named(local_test, 
    c('n_obs', 'terms', 'XtX', 'X1', 'fs_template', 'n_good_vx', 
      'date_created', 'verywise_version', 'XtY', 'psumsY'),
    ignore.order = TRUE)
  
  expect_type(local_test$n_obs, 'integer')
  expect_equal(local_test$n_obs, 30L)
  expect_equal(local_test$fs_template, 'fsaverage')
  expect_equal(local_test$terms, c('(Intercept)', 'sex', 'age'))
  expect_true(local_test$n_good_vx > 0L)
  # expect_lte(local_test$n_good_vx, 10242L)  # fsaverage5 total vertices

  p <- length(local_test$terms)
  n_verts <- count_vertices('fsaverage')

  expect_length(local_test$X1, p)
  expect_true(is.matrix(local_test$XtX))
  expect_equal(dim(local_test$XtX), c(p, p))

  expect_true(inherits(local_test$XtY, 'FBM'))
  expect_equal(local_test$XtY$nrow, p)
  expect_equal(local_test$XtY$ncol, n_verts)

  expect_true(inherits(local_test$psumsY, 'FBM'))
  expect_equal(local_test$psumsY$nrow, 2L)
  expect_equal(local_test$psumsY$ncol, n_verts)

  bk_files <- list.files(outp_dir, 
                         pattern = '\\.bk$', recursive = TRUE)
  expect_gte(length(bk_files), 2L)  # at least XtY and psumsY

  rds_files <- list.files(outp_dir, pattern = '\\.static\\.rds$', recursive = TRUE)
  expect_length(rds_files, 1L)

})


# # ── 4. Numerical sanity ───────────────────────────────────────────────────────

# test_that('XtX is symmetric positive-definite', {
#   outp   <- withr::local_tempdir()
#   result <- run_local(outp)

#   expect_equal(result$XtX, t(result$XtX), tolerance = 1e-6)
#   expect_true(all(eigen(result$XtX, only.values = TRUE)$values > 0))
# })

# test_that('XtY is non-zero over known signal vertices', {
#   outp   <- withr::local_tempdir()
#   result <- run_local(outp)

#   signal_verts <- which(sim_truth$signal_vertices)
#   XtY_signal   <- result$XtY[, signal_verts]

#   expect_true(any(XtY_signal != 0),
#               label = 'XtY should be non-zero at signal vertices')
# })

# test_that('psumsY row 1 (1tY) is positive over signal vertices', {
#   outp   <- withr::local_tempdir()
#   result <- run_local(outp)

#   signal_verts <- which(sim_truth$signal_vertices)
#   col_sums     <- result$psumsY[1, signal_verts]

#   expect_true(all(col_sums > 0),
#               label = 'column sums of Y should be positive (thickness > 0)')
# })

# test_that('psumsY row 2 (YtY) is strictly greater than row 1^2 / n', {
#   # Cauchy-Schwarz: sum(y^2) >= (sum(y))^2 / n
#   outp   <- withr::local_tempdir()
#   result <- run_local(outp)

#   n            <- result$n_obs
#   signal_verts <- which(sim_truth$signal_vertices)
#   sum_y        <- result$psumsY[1, signal_verts]
#   sum_y2       <- result$psumsY[2, signal_verts]

#   expect_true(all(sum_y2 >= sum_y^2 / n - 1e-4),
#               label = 'Cauchy-Schwarz: YtY >= (1tY)^2 / n')
# })


# # ── 5. save_ss behaviour ──────────────────────────────────────────────────────

# test_that('save_ss = TRUE keeps ss directory after function exits', {
#   outp <- withr::local_tempdir()
#   run_local(outp, save_ss = TRUE)

#   ss_dir <- file.path(outp, 'ss')
#   expect_true(dir.exists(ss_dir))
# })

# test_that('save_ss = FALSE cleans up ss directory after function exits', {
#   outp <- withr::local_tempdir()
#   run_local(outp, save_ss = FALSE)

#   ss_dir <- file.path(outp, 'ss')
#   expect_false(dir.exists(ss_dir))
# })

# test_that('save_ss as a custom path saves ss there', {
#   outp    <- withr::local_tempdir()
#   ss_path <- withr::local_tempdir(pattern = 'custom_ss')
#   run_local(outp, save_ss = ss_path)

#   expect_true(dir.exists(ss_path))
#   expect_gt(length(list.files(ss_path, recursive = TRUE)), 0L)
# })


# # ── 7. Reproducibility ────────────────────────────────────────────────────────

# test_that('two runs with the same seed produce identical XtY matrices', {
#   outp1 <- withr::local_tempdir()
#   outp2 <- withr::local_tempdir()

#   r1 <- run_local(outp1, seed = 1L)
#   r2 <- run_local(outp2, seed = 1L)

#   signal_verts <- which(sim_truth$signal_vertices)
#   expect_equal(r1$XtY[, signal_verts], r2$XtY[, signal_verts],
#                tolerance = 1e-5)
# })
