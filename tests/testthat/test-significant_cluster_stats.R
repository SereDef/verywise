# Unit tests for the significant_cluster_stats function
test_that("significant_cluster_stats handles invalid direct ories", {
  subj_dir <- test_path("fixtures", "fs7")

  expect_error(significant_cluster_stats("mean", "invalid_ss_dir", 
                                         file.path(subj_dir, 'results'), 
                                         "wisdom", "area", "lh"), 
               "Super-subject directory does not exist")
  
  expect_error(significant_cluster_stats("mean", file.path(subj_dir, 'ss'),
                                         "invalid_res_dir", "wisdom", "thickness", "rh"), 
               "Results directory does not exist")
})

test_that("We significant_cluster_stats with missing stack_names.txt", {
  subj_dir <- test_path("fixtures", "fs7")

  expect_error(significant_cluster_stats("mean", file.path(subj_dir, 'ss'),
                                         file.path(subj_dir), 
                                         "wisdom", "area", "lh"), 
               "Cannot find the `stack_names.txt`")
})

test_that("significant_cluster_stats with invalid term", {
  subj_dir <- test_path("fixtures", "fs7")

  expect_error(significant_cluster_stats("mean", file.path(subj_dir, 'ss'),
                                         file.path(subj_dir, 'results'), 
                                         "invalid_term", "thickness", "lh"), 
               "Term 'invalid_term' not found in stack_names.txt")
})

test_that("significant_cluster_stats works when clusters are found", {
  subj_dir <- test_path("fixtures", "fs7")
  
  med_area = significant_cluster_stats(stat = 'median',
                            ss_dir = file.path(subj_dir, 'ss'),
                            res_dir = file.path(subj_dir, 'results'),
                            term='wisdom',
                            hemi='lh',
                            measure = 'area')
  
  expect_is(result, "data.frame")
  expect_gt(nrow(result), 0)  # there's at least one row
  expect_gt(ncol(result), 1)  # there are multiple columns
})

test_that("significant_cluster_stats works when clusters are not found", {
  subj_dir <- test_path("fixtures", "fs7")
  
  med_area = significant_cluster_stats(stat = 'median',
                            ss_dir = file.path(subj_dir, 'ss'),
                            res_dir = file.path(subj_dir, 'results'),
                            term='sexMale',
                            hemi='lh',
                            measure = 'area')
  
  expect_null(med_area)  # No significant clusters should return NULL
})
