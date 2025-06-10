
subj_dir <- test_path("fixtures","example_data")

pheno <- read.csv(file.path(subj_dir, "phenotype.csv"))
folder_ids <- pheno$folder_id

outp_dir <- subj_dir
measure <- "area"
hemi <- "lh"
n_cores <- 1
fwhmc <- "fwhm10"
target <- "fsaverage3"
error_cutoff <- 2

# Helper function to remove output files
cleanup_output <- function(outp_dir, hemi, measure, target) {
  pattern <- paste0(hemi, "\\.", measure, "\\.", target, ".*")
  files <- list.files(outp_dir, pattern = pattern, full.names = TRUE)
  file.remove(files)
}

test_that("build_supersubject returns FBM with correct dimensions for valid input", {

  # cleanup_output(outp_dir, hemi, measure, target)

  ss <- build_supersubject(
    subj_dir = subj_dir,
    folder_ids = folder_ids,
    outp_dir = outp_dir,
    measure = measure,
    hemi = hemi,
    n_cores = n_cores,
    fwhmc = fwhmc,
    target = target,
    error_cutoff = error_cutoff,
    save_rds = FALSE,
    verbose = TRUE
  )

  expect_s3_class(ss, "FBM")
  expect_equal(dim(ss), c(length(folder_ids), 642))
})

test_that("build_supersubject errors if missing folders exceed error_cutoff", {

  too_many_missing_folders <- folder_ids
  too_many_missing_folders[1:3] <- c("site1/missing_sub1", "site1/missing_sub2", "site1/missing_sub3")

  expect_error(
    build_supersubject(
      subj_dir = subj_dir,
      folder_ids = too_many_missing_folders,
      outp_dir = outp_dir,
      measure = measure,
      hemi = hemi,
      n_cores = n_cores,
      error_cutoff = error_cutoff,
      save_rds = FALSE,
      verbose = FALSE
    ),
    regexp = "observations specified in phenotype were.*not found"
  )
})

test_that("build_supersubject warns and continues if missing folders under cutoff", {

  some_missing_folders <- folder_ids
  some_missing_folders[1] <- "site1/missing_sub1"

  expect_warning(
    build_supersubject(
      subj_dir = subj_dir,
      folder_ids = some_missing_folders,
      outp_dir = outp_dir,
      measure = measure,
      hemi = hemi,
      n_cores = n_cores,
      error_cutoff = error_cutoff,
      save_rds = FALSE,
      verbose = FALSE
    ),
    regexp = "not found in `subj_dir`"
  )
})

test_that("build_supersubject works with parallel processing", {

  # cleanup_output(outp_dir, hemi, measure, target)
  more_cores = 4

  ss <- build_supersubject(
    subj_dir = subj_dir,
    folder_ids = folder_ids,
    outp_dir = outp_dir,
    measure = measure,
    hemi = hemi,
    n_cores = more_cores,
    fwhmc = fwhmc,
    target = target,
    error_cutoff = error_cutoff,
    save_rds = FALSE,
    verbose = FALSE
  )
  expect_s3_class(ss, "FBM")
  expect_equal(dim(ss), c(length(folder_ids), 642))
})

test_that("build_supersubject creates output files when save_rds = TRUE", {

  # cleanup_output(outp_dir, hemi, measure, target)

  ss <- build_supersubject(
    subj_dir = subj_dir,
    folder_ids = folder_ids,
    outp_dir = outp_dir,
    measure = measure,
    hemi = hemi,
    n_cores = n_cores,
    fwhmc = fwhmc,
    target = target,
    error_cutoff = error_cutoff,
    save_rds = TRUE,
    verbose = FALSE
  )

  # Check that the .bk and .rds files exist
  bk_file <- file.path(outp_dir, paste(hemi, measure, target, "supersubject.bk", sep = "."))
  rds_file <- file.path(outp_dir, paste(hemi, measure, target, "supersubject.rds", sep = "."))
  expect_true(file.exists(bk_file))
  expect_true(file.exists(rds_file))
})
