
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

# test_that("run_vw_meta handles medial wall NAs and valid vertices correctly", {
  
#   # 1. Setup minimal directory structure
#   dir1 <- withr::local_tempdir()
#   dir2 <- withr::local_tempdir()
  
#   stack_df <- data.frame(stack_number = 1, stack_name = "age", stringsAsFactors = FALSE)
#   write.table(stack_df, file.path(dir1, "stack_names.txt"), sep = "\t", row.names = FALSE)
#   write.table(stack_df, file.path(dir2, "stack_names.txt"), sep = "\t", row.names = FALSE)
  
#   file.create(file.path(dir1, "lh.area.stack1.coef.mgh"))
#   file.create(file.path(dir1, "lh.area.stack1.se.mgh"))
#   file.create(file.path(dir2, "lh.area.stack1.coef.mgh"))
#   file.create(file.path(dir2, "lh.area.stack1.se.mgh"))
  
#   # 2. Mock internal verywise utilities
#   local_mocked_bindings(
#     vw_init_message = function(...) invisible(),
#     vw_message = function(...) invisible(),
#     vw_error = function(msg) stop(msg),
#     check_measure = function(...) TRUE,
#     check_path = function(p, ...) p,
#     outcome_name = function(...) "test_outcome",
#     check_numeric_param = function(...) TRUE,
#     check_cores = function(...) 1,
#     mask_cortex = function(...) rep(c(TRUE, FALSE), each = 5), # 10 vertices: 5 cortex, 5 medial wall
#     build_output_bks = function(path, res_bk_names, ...) {
#       setNames(replicate(length(res_bk_names), tempfile(fileext = ".bk")), res_bk_names)
#     },
#     make_chunk_sequence = function(verts, ...) list(verts),
#     with_parallel = function(..., expr) eval(substitute(expr)),
#     init_progress_tracker = function(...) invisible(),
#     convert_to_mgh = function(...) invisible()
#   )
  
#   # 3. Mock load.mgh to return NAs for some vertices (simulating missing study data)
#   local_mocked_bindings(
#     load.mgh = function(file) {
#       val <- rep(0.5, 10)
#       if (grepl("se.mgh", file)) val <- rep(0.1, 10)
#       val[1] <- NA # First vertex missing in all studies
#       return(val)
#     }
#   )
  
#   # Mock cli functions 
#   withr::local_namespace("cli")
#   local_mocked_bindings(
#     cli_progress_step = function(...) invisible(),
#     cli_progress_done = function(...) invisible()
#   )
  
#   # 4. Run Function
#   out <- run_vw_meta(
#     term = "age", hemi = "lh", measure = "area",
#     res_dirs = c(dir1, dir2), study_names = c("S1", "S2"),
#     outp_dir = withr::local_tempdir(), n_cores = 1, verbose = FALSE
#   )
  
#   # 5. Assertions
#   expect_named(out, c("coef", "se", "p"))
  
#   # The first vertex was explicitly NA in the mock, so metafor should be skipped
#   expect_true(is.na(out$coef[1, 1]))
  
#   # The second vertex is valid (cortex, no NA), should have meta-analytic results
#   expect_false(is.na(out$coef[1, 2]))
#   expect_true(out$p[1, 2] >= 0 && out$p[1, 2] <= 1)
  
#   # The 6th vertex is medial wall (is_cortex = FALSE), should remain NA
#   expect_true(is.na(out$coef[1, 6]))
# })