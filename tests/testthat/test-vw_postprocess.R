library(bigstatsr)

set.seed(42)

# Create small example FBMs
n_verts <- 10
n_terms <- 3
n_obsvt <- 10
mat_type <- "float"

coefs <- c(1.5, 2.5, 0.00000001)
ses <- c(0.2, 0.01, 2)
ps <- c(1, 0.00001, 0.00000000000000000000001)
resids <- rnorm(n_obsvt)

logps <- -1 * log10(ps)

vw_results <- list(
  coef = FBM(nrow=n_terms, ncol=n_verts, type=mat_type, init = rep(coefs, n_verts)),
    se = FBM(nrow=n_terms, ncol=n_verts, type=mat_type, init = rep(ses, n_verts)),
     p = FBM(nrow=n_terms, ncol=n_verts, type=mat_type, init = rep(ps, n_verts)),
 resid = FBM(nrow=n_obsvt, ncol=n_verts, type=mat_type, init = rep(resids, n_verts))
)

test_that("convert_to_mgh writes correct output for basic stats", {

  tmpdir <- withr::local_tempdir()
  dir.create(tmpdir, recursive = TRUE)

  result_path <- file.path(tmpdir, "stats/hemi.measure")

  # Mock fbm2mgh to capture calls instead of writing .mgh
  called <- list()
  mock_fbm2mgh <- function(fbm, fnames, mode) {
    called <<- append(called, list(list(fnames = fnames, mode = mode, data = fbm[])))
    invisible(NULL)
  }

  # Run with mock
  with_mock(
    fbm2mgh = mock_fbm2mgh,
    convert_to_mgh(vw_results, result_path, stacks = 1:n_terms)
  )

  # Check that residuals use correct mode
  resid_call <- called[[which(vapply(called, function(x) any(
    grepl("residuals.mgh", x$fnames)), logical(1)))]]
  expect_equal(resid_call$mode, "allrows.1file")

  # Check that coef and se are untransformed
  coef_call <- called[[which(vapply(called, function(x) any(
    grepl("coef", x$fnames)), logical(1)))]]
  expected <- matrix(rep(coefs, n_verts), nrow = n_terms)
  expect_equal(coef_call$data, expected)

  # Check that -log10p transformation works
  log10p_call <- called[[which(vapply(called, function(x) any(
    grepl("-log10p", x$fnames)), logical(1)))]]
  expected <- matrix(rep(logps, n_verts), nrow = n_terms)
  expect_equal(log10p_call$data, expected, tolerance = 1e-10)
})

# ==============================================================================

test_that("move_result_files correctly matches and copies files", {
  # Create temporary from/to directories
  from_dir <- file.path(tempdir(), "from_test")
  dir.create(from_dir, recursive = TRUE, showWarnings = FALSE)

  to_dir   <- file.path(tempdir(), "to_test")
  dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)

  # Create mock file structure in from_dir
  files <- c(
    "lh.area.stack1.cache.th30.abs.sig.ocn.mgh", # matches clusters pattern
    "rh.volume.stack12.coef.mgh",                # matches coef pattern
    "stack_names.txt",                           # matches names file
    "random.txt"                                 # should not match
  )

  # Create them inside some sub-directories
  for (subdir in c('model1', 'model2')) {
    for (file in files) {
      dir.create(file.path(from_dir, subdir), recursive = TRUE, showWarnings = FALSE)
      writeLines("testdata", file.path(from_dir, subdir, file))
    }
  }

  # Run the function
  result <- move_result_files(from_dir, to_dir)

  # All matching files should have been copied
  expect_true(all(result))

  # Check that only matched files exist in to_dir
  copied_files <- list.files(to_dir, recursive = TRUE)
  expect_true(all(copied_files != "random.txt"))

  # Check that sub-directory structure is preserved
  expect_true(file.exists(file.path(to_dir, "model1", "stack_names.txt")))
  expect_true(file.exists(file.path(to_dir, "model2", "stack_names.txt")))

  expect_true(file.exists(file.path(to_dir, "model1", "lh.area.stack1.cache.th30.abs.sig.ocn.mgh")))
  expect_true(file.exists(file.path(to_dir, "model2", "rh.volume.stack12.coef.mgh")))

})

test_that("move_result_files returns invisible logical vector", {

  from_dir <- file.path(tempdir(), "from_empty_test")
  dir.create(from_dir, showWarnings = FALSE)

  to_dir   <- file.path(tempdir(), "to_empty_test")
  dir.create(to_dir, showWarnings = FALSE)

  ret <- move_result_files(from_dir, to_dir)
  expect_type(ret, "logical")
  expect_length(ret, 0)
})
