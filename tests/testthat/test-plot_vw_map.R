# tests/testthat/test-plot_vw_map.R

# ── fixtures ─────────────────────────────────────────────────────────────────
subj_dir <- test_path("fixtures", "fs7")
res_dir <- file.path(subj_dir, "results")

# Spy that captures args passed to plot_vw_surf without invoking Python
fake_plot_vw_surf <- function(lh = NULL, rh = NULL, ...) {
  list(lh = lh, rh = rh, ...)
}

# ── helpers ───────────────────────────────────────────────────────────────────

capture_surf_args <- function(...) {
  result <- NULL
  local_mocked_bindings(
    plot_vw_surf = function(...) { result <<- list(...); invisible("fake.html") }
  )
  plot_vw_map(res_dir = res_dir, ...)
  result
}


# ── Input validation ──────────────────────────────────────────────────────────

test_that("errors if res_dir does not exist", {
  expect_snapshot(
    plot_vw_map("/nonexistent/path", term = "age"),
    error = TRUE
  )
})

test_that("errors if stack_names.txt is missing", {
  dir <- withr::local_tempdir()
  expect_snapshot(
    plot_vw_map(dir, term = "age"),
    error = TRUE
  )
})

test_that("errors if term is absent from stack_names.txt", {
  expect_snapshot(
    plot_vw_map(res_dir, term = "__no_such_term__"),
    error = TRUE
  )
})


# ── File resolution ───────────────────────────────────────────────────────────

test_that("errors when coef file is missing for requested hemisphere", {
  dir <- withr::local_tempdir()
  file.copy(file.path(res_dir, "stack_names.txt"), dir)

  local_mocked_bindings(
    plot_vw_surf = function(...) invisible("fake.html")
  )

  expect_snapshot(
    plot_vw_map(dir, term = "age", hemi = "lh", threshold = NULL),
    error = TRUE
  )
})

test_that("warns and plots unmasked when OCN file is absent", {
  dir <- withr::local_tempdir()
  coef_files <- list.files(res_dir, pattern = "\\.coef\\.mgh$",
                           full.names = TRUE)
  file.copy(c(coef_files, file.path(res_dir, "stack_names.txt")), dir)

  local_mocked_bindings(
    plot_vw_surf = function(...) invisible("fake.html")
  )

  expect_snapshot(
    plot_vw_map(dir, term = "age", hemi = "lh", threshold = "cws")
  )
})


# ── CWS masking ───────────────────────────────────────────────────────────────

test_that("cws masking sets non-significant vertices to NA", {
  args <- capture_surf_args(term = "age", hemi = "lh", threshold = "cws")

  expect_snapshot(which(is.na(args$lh)))
  expect_snapshot(args$lh[!is.na(args$lh)])
})

test_that("cws masking preserves coefficient values for significant vertices", {
  ids <- utils::read.table(
    file.path(res_dir, "stack_names.txt"),
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
  )
  stack <- paste0("stack", ids[ids$stack_name == "age", "stack_number"])

  raw_coef <- load.mgh(
    file.path(res_dir, paste("lh", "area", stack, "coef.mgh", sep = "."))
  )$x

  ocn_file <- list.files(res_dir,
    pattern = paste0("^lh\\.area\\.", stack, "\\.cache\\..*\\.sig\\.ocn\\.mgh$"),
    full.names = TRUE
  )
  raw_ocn <- load.mgh(ocn_file)$x

  args <- capture_surf_args(term = "age", hemi = "lh", threshold = "cws")

  sig_idx <- which(raw_ocn != 0)
  expect_equal(args$lh[sig_idx], raw_coef[sig_idx])
})


# ── Hemisphere routing ────────────────────────────────────────────────────────

test_that("hemi = 'lh' loads only left hemisphere", {
  args <- capture_surf_args(term = "age", hemi = "lh", threshold = NULL)


  expect_null(args$rh)
  expect_snapshot(length(args$lh))
})

test_that("hemi = 'rh' loads only right hemisphere", {

  rh_files <- list.files(res_dir, pattern = "^rh\\..*\\.coef\\.mgh$")
  skip_if(length(rh_files) == 0, "No rh coef fixture available")

  args <- capture_surf_args(term = "age", hemi = "rh", threshold = NULL)
  expect_null(args$lh)
  expect_snapshot(length(args$rh))
})

test_that("hemi = 'both' loads both hemispheres", {
  args <- capture_surf_args(term = "age", hemi = "both", threshold = NULL)

  # lh is always present in fixture
  expect_false(is.null(args$lh))
  expect_snapshot(length(args$lh))

  rh_files <- list.files(res_dir, pattern = "^rh\\..*\\.coef\\.mgh$")
  if (length(rh_files) > 0) {
    expect_false(is.null(args$rh))
    expect_snapshot(length(args$rh))
  } else {
    expect_null(args$rh)
  }
})


# ── Threshold forwarding ──────────────────────────────────────────────────────

test_that("numeric threshold is forwarded unchanged", {
  args <- capture_surf_args(term = "age", hemi = "lh", threshold = 0.05)
  expect_snapshot(args$threshold)
})

test_that("threshold = NULL forwards NULL", {
  args <- capture_surf_args(term = "age", hemi = "lh", threshold = NULL)
  expect_null(args$threshold)
})

test_that("threshold = 'cws' forwards NULL to plot_vw_surf", {
  args <- capture_surf_args(term = "age", hemi = "lh", threshold = "cws")
  expect_null(args$threshold)
})


# ── Title and ... forwarding ──────────────────────────────────────────────────

test_that("default title is built from term and measure", {
  args <- capture_surf_args(term = "age", measure = "area",
                            hemi = "lh", threshold = NULL)
  expect_snapshot(args$title)
})

test_that("user title via ... overrides default and is not duplicated", {
  expect_no_error(
    args <- capture_surf_args(term = "age", hemi = "lh",
                              threshold = NULL, title = "My title")
  )
  expect_snapshot(args$title)
})

test_that("extra ... args reach plot_vw_surf", {
  args <- capture_surf_args(term = "age", hemi = "lh",
                            threshold = NULL, cmap = "RdBu_r", dpi = 300L)
  expect_snapshot(args$cmap)
  expect_snapshot(args$dpi)
})

test_that("surface is forwarded to plot_vw_surf", {
  args <- capture_surf_args(term = "age", hemi = "lh",
                            threshold = NULL, surface = "inflated")
  expect_snapshot(args$surface)
})