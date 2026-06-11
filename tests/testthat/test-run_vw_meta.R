
library(metafor)

# Example data
set.seed(3108)

res_dir <- test_path("fixtures", "fs7", "results")
out_dir <- test_path("fixtures", "rma")

test_that("run_vw_meta runs end-to-end with random terms", {

  # I meta-analyse two terms from the same model, which makes no sense but 
  # only for the pipeline to run 

  result <- run_vw_meta(
      term = c("age", "wisdom"), 
      hemi = "lh", 
      measure = "area",
      study_names = c("S1", "S2"),
      study_weights = NULL,
      res_dirs = c(res_dir, res_dir), 
      outp_dir = out_dir,
      fs_template = "fsaverage",
      n_cores = 1, 
      verbose = TRUE
    )
  
  # Structure tests
  expect_type(result, "list")
  expect_named(result, c("coef", "se", "p"))
  expect_s4_class(result$coef, "FBM")
  expect_s4_class(result$se, "FBM")
  expect_s4_class(result$p, "FBM")

  expect_true(file.exists(
    file.path(out_dir, 'lh.area.age.coef.mgh')))
  expect_true(file.exists(
    file.path(out_dir, 'lh.area.age.p.mgh')))
  expect_true(file.exists(
    file.path(out_dir, 'lh.area.age.fdr.mgh')))
  
  
})