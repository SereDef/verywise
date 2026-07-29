
# ── Fixtures ──────────────────────────────────────────────────────────────────
subj_dir <- test_path("fixtures", "fs7")
res_dir <- file.path(subj_dir, "results")

# Pre-load real fixture data for testing the renderer directly
lh_coef_file <- list.files(res_dir, pattern = "^lh\\..*\\.coef\\.mgh$", full.names = TRUE)[1]
rh_coef_file <- list.files(res_dir, pattern = "^rh\\..*\\.coef\\.mgh$", full.names = TRUE)[1]

# ── Setup ─────────────────────────────────────────────────────────────────────
skip_if_not_installed("reticulate")

# ── Input validation ──────────────────────────────────────────────────────────

test_that("plot_vw_surf validates missing inputs", {
  expect_snapshot(
    plot_vw_surf(),
    error = TRUE
  )
})

test_that("plot_vw_surf validates output extension", {
  expect_snapshot(
    plot_vw_surf(lh = lh_coef_file, to_file = "plot.pdf"),
    error = TRUE
  )
})

test_that("plot_vw_surf validates views and ROIs", {
  expect_snapshot(
    plot_vw_surf(lh = lh_coef_file, views = c("lateral", "invalid_view")),
    error = TRUE
  )
  
  expect_snapshot(
    plot_vw_surf(lh = lh_coef_file, roi_outline = "fake_roi", fs_home = subj_dir),
    error = TRUE
  )
})

# ── File Generation ───────────────────────────────────────────────────────────

test_that("generates an interactive HTML file when to_file is NULL", {
  res <- plot_vw_surf(lh = lh_coef_file, fs_template = "fsaverage")
  
  expect_type(res, "character")
  expect_true(file.exists(res))
  expect_match(basename(res), "\\.html$")
})

test_that("generates a static PNG file when requested", {
  tmp_png <- withr::local_tempfile(fileext = ".png")
  
  res <- plot_vw_surf(
    lh = lh_coef_file,
    rh = if (!is.na(rh_coef_file)) rh_coef_file else NULL,
    to_file = tmp_png,
    views = c("lateral", "medial"),
    cmap = "RdBu_r",
    title = "Test Static Plot",
    dpi = 50L # Low DPI for faster testing
  )
  
  expect_equal(res, tmp_png)
  expect_true(file.exists(res))
  expect_gt(file.size(res), 0)
})