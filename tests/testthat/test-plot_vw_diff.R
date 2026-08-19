# Plotting tests require a stable local Python/nilearn/FreeSurfer setup
skip_on_ci()
skip_if_not_installed("reticulate")
skip_if_no_freesurfer()

# ── Fixtures ──────────────────────────────────────────────────────────────────
subj_dir <- test_path("fixtures", "fs7")
res_dir <- file.path(subj_dir, "results")

# ── Setup ─────────────────────────────────────────────────────────────────────

# Find two MGH files from fixtures to compute diffs if available
lh_files <- list.files(res_dir, pattern = "^lh\\..*\\.coef\\.mgh$", full.names = TRUE)

n_vert <- 10242 # fsaverage length

# Dummy maps for numeric vector testing
lh_a_dummy <- rep(1.0, n_vert) 
lh_b_dummy <- rep(0.5, n_vert)
rh_a_dummy <- rep(2.0, n_vert)
rh_b_dummy <- rep(1.0, n_vert)

# ── Input validation & Errors ─────────────────────────────────────────────────

test_that("vw_diff errors when maps have mismatched lengths", {
  expect_snapshot(
    plot_vw_diff(lh_a = lh_a_dummy, lh_b = c(1, 2, 3)),
    error = TRUE
  )
})

test_that("vw_diff errors when one side is NULL but the other is provided", {
  expect_snapshot(
    plot_vw_diff(lh_a = lh_a_dummy, lh_b = NULL),
    error = TRUE
  )
})

test_that("vw_diff errors with invalid object types", {
  expect_snapshot(
    plot_vw_diff(lh_a = data.frame(a = 1), lh_b = lh_b_dummy),
    error = TRUE
  )
})

# ── Computation & Plot Orchestration ──────────────────────────────────────────

test_that("computes accurate difference maps from numeric vectors", {
  # Capture the console output using expect_snapshot, which will also run the function
  expect_snapshot({
    res <- plot_vw_diff(
      lh_a = lh_a_dummy, 
      lh_b = lh_b_dummy,
      rh_a = rh_a_dummy,
      rh_b = rh_b_dummy,
      label_a = "Condition A",
      label_b = "Condition B",
      fs_template = 'fsaverage5'
    )
  })
})

test_that("computes difference maps from MGH file paths", {
  skip_if(length(lh_files) < 2, "Need at least two left-hemisphere MGH fixtures")
  
  # We test the static output generation using real MGH files
  tmp_png <- withr::local_tempfile(fileext = ".png")
  
  res_plot <- plot_vw_diff(
    lh_a = lh_files[1],
    lh_b = lh_files[2],
    label_a = "Model 1",
    label_b = "Model 2",
    fs_template = 'fsaverage5',
    to_file = tmp_png,
    dpi = 50L
  )
  
  expect_true(file.exists(res_plot))
  expect_match(basename(res_plot), "\\.png$")
  expect_gt(file.size(res_plot), 0)
})