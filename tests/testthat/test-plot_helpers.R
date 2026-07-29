# ── Fixtures ──────────────────────────────────────────────────────────────────
subj_dir <- test_path("fixtures", "fs7")
res_dir <- file.path(subj_dir, "results")
lh_coef_file <- list.files(res_dir, pattern = "^lh\\..*\\.coef\\.mgh$", full.names = TRUE)[1]

# ── check_hemi & check_mask ───────────────────────────────────────────────────

test_that("check_hemi resolves file paths correctly", {
  skip_if_not(file.exists(lh_coef_file), "Missing lh coef fixture")
  
  res <- check_hemi(lh_coef_file, "fsaverage")
  expect_type(res, "double")
  expect_gt(length(res), 0)
})

test_that("check_hemi errors on nonexistent files", {
  
  expect_snapshot(
    check_hemi("nonexistent_file.mgh", "fsaverage"),
    error = TRUE
  )
})

test_that("check_hemi errors if vector length < declared fsaverage", {
  short_vec <- rep(1, 100)
  expect_error(
    res <- check_hemi(short_vec, "fsaverage")
  )
})

test_that("check_hemi warns on mismatched vector lengths", {
  long_vec <- rep(1, 10000)
  expect_snapshot(
    res <- check_hemi(long_vec, "fsaverage3")
  )
  expect_equal(res, long_vec)
})

test_that("check_mask validates logical inputs", {
  valid_mask <- rep(TRUE, 163842)
  expect_equal(check_mask(valid_mask, "fsaverage"), valid_mask)
  
  expect_snapshot(
    check_mask(c(1, 2, 3), "fsaverage"),
    error = TRUE
  )
})

# ── load_and_mask_coef ────────────────────────────────────────────────────────

test_that("load_and_mask_coef warns and skips missing files", {
  expect_snapshot(
    res <- load_and_mask_coef("lh", "area", "stack999", res_dir, threshold = NULL)
  )
  expect_null(res$coef)
  expect_null(res$mask)
})

test_that("load_and_mask_coef handles numeric thresholds", {
  skip_if_not(file.exists(lh_coef_file), "Missing lh coef fixture")
  
  # We manually pull the stack ID from the file name to mock the call
  stack_id <- gsub(".*\\.(stack[0-9]+)\\.coef\\.mgh", "\\1", basename(lh_coef_file))
  
  res <- load_and_mask_coef("lh", "area", stack_id, res_dir, threshold = 0.5)
  expect_false(is.null(res$coef))
  expect_type(res$mask, "logical")
})

test_that("load_and_mask_coef applies CWS thresholds appropriately", {
  skip_if_not(file.exists(lh_coef_file), "Missing lh coef fixture")
  stack_id <- gsub(".*\\.(stack[0-9]+)\\.coef\\.mgh", "\\1", basename(lh_coef_file))
  
  # When threshold="cws" and OCN file exists, it applies mask. 
  # If no OCN file exists, it will trigger a message and return unmasked.
  expect_snapshot({
    res <- load_and_mask_coef("lh", "area", stack_id, res_dir, threshold = "cws")
  })
  
  expect_false(is.null(res$coef))
})

# ── resolve_mesh ──────────────────────────────────────────────────────────────

test_that("resolve_mesh finds local FreeSurfer dir if it exists", {
  # Mock a local FreeSurfer subjects directory
  mock_fs_home <- withr::local_tempdir()
  surf_dir <- file.path(mock_fs_home, "subjects", "fsaverage", "surf")
  dir.create(surf_dir, recursive = TRUE)

  scrub_tmp <- function(x) {
    gsub(paste0(mock_fs_home), "<FS_HOME>", x, fixed = TRUE)
  }

  expect_snapshot({
    res <- resolve_mesh("fsaverage", mock_fs_home)
  }, transform = scrub_tmp)

  expect_equal(res, mock_fs_home)

})

test_that("resolve_mesh returns NULL (Python download flag) if mesh is missing", {
  mock_fs_home <- withr::local_tempdir()
  
  expect_snapshot({
    res <- resolve_mesh("fsaverage", mock_fs_home)
  })
  expect_null(res)
})