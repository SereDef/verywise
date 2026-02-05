
subj_dir <- test_path("fixtures", "fs3")

pheno <- read.csv(file.path(subj_dir, "phenotype.csv"))
folder_ids <- pheno$folder_id

supsubj_dir <- file.path(subj_dir, "ss")
measure <- "area"
hemi <- "lh"
n_cores <- 1
fwhmc <- "fwhm10"
fs_template <- "fsaverage3"
error_cutoff <- 2

supsubj_file <- paste(hemi, measure, fs_template, "supersubject.rds", sep='.')

# Helper function to remove output files
# cleanup_output <- function(supsubj_dir, hemi, measure, fs_template) {
#   pattern <- paste0(hemi, "\\.", measure, "\\.", fs_template, ".*")
#   files <- list.files(supsubj_dir, pattern = pattern, full.names = TRUE)
#   file.remove(files)
# }

test_that("build_supersubject returns FBM with correct dimensions for valid input", {

  # cleanup_output(supsubj_dir, hemi, measure, fs_template)

  ss <- build_supersubject(
    subj_dir = subj_dir,
    folder_ids = folder_ids,
    supsubj_dir = supsubj_dir,
    measure = measure,
    hemi = hemi,
    n_cores = n_cores,
    fwhmc = fwhmc,
    fs_template = fs_template,
    error_cutoff = error_cutoff,
    save_rds = FALSE,
    verbose = FALSE
  )

  expect_s4_class(ss, "FBM")
  expect_equal(dim(ss), c(length(folder_ids), 642))
})

test_that("build_supersubject errors if missing folders exceed error_cutoff", {

  # skip_on_os("windows")

  too_many_missing_folders <- folder_ids
  too_many_missing_folders[1:3] <- c("site1/missing_sub1",
                                     "site1/missing_sub2",
                                     "site1/missing_sub3")

  expect_error(
    build_supersubject(
      subj_dir = subj_dir,
      folder_ids = too_many_missing_folders,
      supsubj_dir = supsubj_dir,
      measure = measure,
      hemi = hemi,
      fs_template = fs_template,
      n_cores = n_cores,
      error_cutoff = error_cutoff,
      save_rds = FALSE,
      verbose = FALSE
    ),
    regexp = "observations specified in phenotype were.*not found"
  )
})

test_that("build_supersubject warns and continues if missing folders under cutoff", {

  # skip_on_os("windows")

  some_missing_folders <- folder_ids
  some_missing_folders[1] <- "site1/missing_sub1"

  expect_warning(
    build_supersubject(
      subj_dir = subj_dir,
      folder_ids = some_missing_folders,
      supsubj_dir = supsubj_dir,
      measure = measure,
      hemi = hemi,
      fs_template = fs_template,
      n_cores = n_cores,
      error_cutoff = error_cutoff,
      save_rds = FALSE,
      verbose = FALSE
    ),
    regexp = "not found in `subj_dir`"
  )
})

test_that("build_supersubject works with parallel processing", {

  # cleanup_output(supsubj_dir, hemi, measure, fs_template)
  more_cores <- min(2, parallel::detectCores(logical = FALSE))
  if (more_cores < 2) skip("Not enough cores for parallel test")
  
  withr::local_envvar(c(
    OPENBLAS_NUM_THREADS = 1,
    OMP_NUM_THREADS = 1,
    MKL_NUM_THREADS = 1,
    VECLIB_MAXIMUM_THREADS = 1,
    NUMEXPR_NUM_THREADS = 1
  ))

  ss <- build_supersubject(
    subj_dir = subj_dir,
    folder_ids = folder_ids,
    supsubj_dir = supsubj_dir,
    measure = measure,
    hemi = hemi,
    n_cores = more_cores,
    fwhmc = fwhmc,
    fs_template = fs_template,
    error_cutoff = error_cutoff,
    save_rds = FALSE,
    verbose = FALSE
  )
  expect_s4_class(ss, "FBM")
  expect_equal(dim(ss), c(length(folder_ids), 642))
})

test_that("build_supersubject creates output files when save_rds = TRUE", {

  # cleanup_output(supsubj_dir, hemi, measure, fs_template)
  supsubj_dir <- file.path(subj_dir, 'ss')

  ss <- build_supersubject(
    subj_dir = subj_dir,
    folder_ids = folder_ids,
    supsubj_dir = supsubj_dir,
    measure = measure,
    hemi = hemi,
    n_cores = n_cores,
    fwhmc = fwhmc,
    fs_template = fs_template,
    error_cutoff = error_cutoff,
    save_rds = TRUE,
    verbose = FALSE
  )

  # Check that the .bk and .rds files exist
  bk_file <- file.path(supsubj_dir, paste(hemi, measure, fs_template,
                                          "supersubject.bk", sep = "."))
  rds_file <- file.path(supsubj_dir, paste(hemi, measure, fs_template,
                                           "supersubject.rds", sep = "."))
  row_file <- file.path(supsubj_dir, paste(hemi, measure,
                                           "ss.rownames.csv", sep = "."))
  expect_true(file.exists(bk_file))
  expect_true(file.exists(rds_file))
  expect_true(file.exists(row_file))
})

# ==============================================================================

test_that("subset_supersubject works with matching IDs", {

  new_supsubj_dir <- file.path(tempdir(), 'subset')
  id_subset_idx <- c(1, 3, 5, 10)
  id_subset <- folder_ids[id_subset_idx]

  new_ss <- subset_supersubject(
    supsubj_dir = supsubj_dir,
    supsubj_file = supsubj_file,
    folder_ids = id_subset,
    new_supsubj_dir = new_supsubj_dir

  )

  og_ss <- bigstatsr::big_attach(file.path(supsubj_dir, supsubj_file))

  expect_s4_class(new_ss, "FBM")
  expect_equal(og_ss[id_subset_idx, ], new_ss[])

  new_names <- scan(file.path(new_supsubj_dir, "lh.area.ss.rownames.csv"),
                    what = character(), sep = "\n", quiet = TRUE)

  expect_equal(new_names, id_subset)

  expect_true(file.exists(file.path(
    new_supsubj_dir, gsub('.rds$', '.bk', supsubj_file))))
})

test_that("subset_supersubject returns original matrix when no subsetting is required", {

  new_supsubj_dir <- file.path(tempdir(), 'no_subset')

  no_change_ss <- subset_supersubject(
    supsubj_dir = supsubj_dir,
    supsubj_file = supsubj_file,
    folder_ids = folder_ids,
    new_supsubj_dir = new_supsubj_dir

  )

  og_ss <- bigstatsr::big_attach(file.path(supsubj_dir, supsubj_file))

  expect_s4_class(no_change_ss, "FBM")
  expect_equal(og_ss[], no_change_ss[])

  # No ss.rownames or new backing file is created
  expect_false(file.exists(file.path(
    new_supsubj_dir, gsub('.rds$', '.bk', supsubj_file))))
})

test_that("subset_supersubject warns when some IDs missing", {

  new_supsubj_dir <- tempdir()

  expect_warning(
    subset_supersubject(
      supsubj_dir = supsubj_dir,
      supsubj_file = supsubj_file,
      folder_ids = c(folder_ids[1:6], 'not_a_person'), # One missing
      new_supsubj_dir = new_supsubj_dir
    ),
    regexp = "observations specified in phenotype were not found"
  )
})

test_that("subset_supersubject errors when too many IDs are missing", {

  new_supsubj_dir <- tempdir()

  expect_error(
    subset_supersubject(
      supsubj_dir = supsubj_dir,
      supsubj_file = supsubj_file,
      folder_ids = c(folder_ids[1:6], 1:25), # too many missings
      new_supsubj_dir = new_supsubj_dir
    ),
    "observations specified in phenotype were not found"
  )

})

