ds <- list(
  "cohort1" = list(sessions = c("01"), n_subjects = 1)
)

ds_multicohort <- list(
  "cohort1" = list(sessions = c("01", "02"), n_subjects = 5),
  "cohort2" = list(sessions = c("01"), n_subjects = 3)
)

test_that("simulate_long_pheno_data returns correct structure", {
  pheno <- simulate_long_pheno_data(data_structure = ds_multicohort, seed = 123)
  expect_true(is.data.frame(pheno))
  expect_true(all(c("id", "time", "sex", "age", "site", "folder_id") %in% names(pheno)))
  expect_equal(length(unique(pheno$site)), 2)
  expect_equal(nrow(pheno), 5*2 + 3*1)
})

test_that("simulate_dataset creates phenotype file and calls simulate_freesurfer_data", {
  tmp <- tempfile("simdata")
  simulate_dataset(
    path = tmp,
    data_structure = ds,
    overwrite = TRUE,
    seed = 42
  )
  expect_true(file.exists(file.path(tmp, "phenotype.csv")))
  pheno <- utils::read.csv(file.path(tmp, "phenotype.csv"))
  expect_true(is.data.frame(pheno))
  expect_true(nrow(pheno) == 1)
  # Check FreeSurfer folder structure
  expect_true(dir.exists(file.path(tmp, "cohort1", "sub-1_ses-01", "surf")))
})

test_that("simulate_freesurfer_data creates correct files", {
  tmp <- tempfile("fsdata")
  simulate_freesurfer_data(
    path = tmp,
    data_structure = ds,
    vw_resolution = 10,
    measure = "thickness",
    hemi = "lh",
    seed = 1
  )
  surf_dir <- file.path(tmp, "cohort1", "sub-1_ses-01", "surf")
  files <- list.files(surf_dir, pattern = "lh.thickness.*\\.mgh$")
  expect_true(length(files) == 1)
})

test_that("simulate_dataset respects overwrite argument", {
  tmp <- tempfile("simdata2")
  simulate_dataset(
    path = tmp,
    data_structure = ds,
    overwrite = TRUE,
    seed = 1
  )
  old_time <- file.info(file.path(tmp, "phenotype.csv"))$mtime
  Sys.sleep(1)

  simulate_dataset(
    path = tmp,
    data_structure = ds,
    overwrite = FALSE,
    seed = 1
  )
  new_time <- file.info(file.path(tmp, "phenotype.csv"))$mtime
  expect_equal(old_time, new_time)
})

test_that("simulate_freesurfer_data handles simulate_association argument", {
  tmp <- tempfile("fsdata2")
  # Should not error even if simulate_association is a vector
  expect_error(
    simulate_freesurfer_data(
      path = tmp,
      data_structure = ds,
      vw_resolution = 10,
      simulate_association = rep(0.1, 2),
      seed = 1
    ),
    NA
  )
})
