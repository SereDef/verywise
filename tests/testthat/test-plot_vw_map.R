test_that("plot_vw_map basic input validation works", {
  skip_if(is.null(Sys.getenv("FREESURFER_HOME")))
  skip_if_not_installed("reticulate")

  # Skip if required Python modules are not available
  skip_if_not(reticulate::py_module_available("nilearn"), "nilearn not available")
  skip_if_not(reticulate::py_module_available("nibabel"), "nibabel not available")

  res_dir <- test_path("fixtures", "fs7", "results")
  # test_formula <- vw_area ~ sex + age + wisdom + (1 | id)
  fs_home = Sys.getenv("FREESURFER_HOME") # "/Applications/freesurfer/7.4.1" # mac only

  # Function should error if invalid term provided
  expect_error(plot_vw_map(res_dir, "nonexistent_term", fs_home = fs_home),
                           regexp = "Term 'nonexistent_term' not found")

  # Function should return a Python object when everything is valid
  result <- plot_vw_map(
    res_dir = res_dir,
    term = "wisdom",
    hemi = "lh",
    measure = "area",
    surface = "inflated",
    fs_home = fs_home,
    show_in_browser = FALSE
  )

  expect_true(inherits(result, "python.builtin.object"))
  expect_true(inherits(result, "nilearn.plotting.displays._figures.SurfaceFigure"))
})
