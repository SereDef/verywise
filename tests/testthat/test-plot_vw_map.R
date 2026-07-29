# ── Fixtures ──────────────────────────────────────────────────────────────────
subj_dir <- test_path("fixtures", "fs7")
res_dir <- file.path(subj_dir, "results")

# ── Setup ─────────────────────────────────────────────────────────────────────
skip_if_not_installed("reticulate")

# ── Input validation & File Resolution ────────────────────────────────────────

test_that("errors if res_dir does not exist", {
  expect_snapshot(
    plot_vw_map("/nonexistent/path", term = "age"),
    error = TRUE
  )
})

test_that("errors if term is absent from stack_names.txt", {
  expect_snapshot(
    plot_vw_map(res_dir, term = "__no_such_term__"),
    error = TRUE
  )
})

test_that("errors when coef file is missing for requested hemisphere", {
  dir <- withr::local_tempdir()
  file.copy(file.path(res_dir, "stack_names.txt"), dir)
  
  expect_snapshot(
    plot_vw_map(dir, term = "age", hemi = "lh", threshold = NULL),
    error = TRUE,
    transform = function(x) gsub(dir, "<TEMP_DIR>", x, fixed = TRUE)
  )
})

test_that("warns and continues when stack_names.txt is missing (meta-analysis mode)", {
  dir <- withr::local_tempdir()

  file.copy(
    file.path(res_dir, "lh.area.stack3.coef.mgh"),
    file.path(dir, "lh.area.age.coef.mgh")
  )

  expect_snapshot(
    res <- plot_vw_map(dir, term = "age", hemi = "lh", threshold = NULL),
    transform = function(x) gsub(dir, "<TEMP_DIR>", x, fixed = TRUE)
  )

  expect_true(file.exists(res))
})

# ── End-to-end Orchestration ──────────────────────────────────────────────────

test_that("successfully generates an interactive plot with CWS masking", {
  # This tests the default threshold = 'cws' routing without mocking
  res <- plot_vw_map(
    res_dir = res_dir, 
    term = "age", 
    hemi = "lh",
    fs_template = "fsaverage5"
  )
  
  expect_true(file.exists(res))
  expect_match(basename(res), "\\.html$")
})

test_that("successfully generates a static plot with numeric thresholding", {
  tmp_png <- withr::local_tempfile(fileext = ".png")
  
  res <- plot_vw_map(
    res_dir = res_dir, 
    term = "age", 
    hemi = "both", 
    threshold = 0.05,
    to_file = tmp_png,
    cmap = "viridis",
    dpi = 50L
  )
  
  expect_true(file.exists(res))
  expect_match(basename(res), "\\.png$")
  expect_gt(file.size(res), 0)
})