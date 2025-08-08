ds <- list(
  "cohort1" = list(sessions = c("01"), n_subjects = 2)
)
ds_assoc_num <- c(5.5, 4.1)
ds_assoc_chr <- '0.3 * age'

ds_multicohort <- list(
  "cohort1" = list(sessions = c("01", "02"), n_subjects = 5),
  "cohort2" = list(sessions = c("01"), n_subjects = 3)
)

verbosity = FALSE

test_that("simulate_long_pheno_data returns correct structure", {
  pheno <- simulate_long_pheno_data(data_structure = ds_multicohort,
                                    seed = 123, verbose = verbosity)

  expect_s3_class(pheno, "data.frame")
  expect_true(
    all(c("id", "time", "sex", "age", "site", "folder_id") %in% names(pheno)))
  expect_equal(nrow(pheno), 5*2 + 3*1)
  expect_equal(length(unique(pheno$site)), 2)
  # folder_id matches expected pattern
  expect_true(all(grepl("^cohort[0-9]{1}/sub-[0-9]+_ses-[0-9]{2}$",
                        pheno$folder_id)))
})

test_that("simulate_freesurfer_data creates correct file structure", {
  tmp <- file.path(tempdir(), "sym_fs_test")

  simulate_freesurfer_data(
    path = tmp,
    data_structure = ds,
    fs_template = 'fsmicro',
    measure = "thickness",
    hemi = "lh",
    seed = 123,
    verbose = verbosity
  )
  # directory exists
  expect_true(dir.exists(tmp))

  # subject directories exist
  surf_dir <- file.path(tmp, "cohort1", "sub-1_ses-01", "surf")
  files <- list.files(surf_dir, pattern = "lh.thickness.*\\.mgh$")
  expect_true(length(files) == 1)
})

test_that("simulate_dataset creates phenotype file and calls simulate_freesurfer_data", {
  tmp <- file.path(tempdir(), "sym_dset_test")

  simulate_dataset(
    path = tmp,
    data_structure = ds,
    overwrite = TRUE,
    seed = 321,
    verbose = verbosity
  )
  # phenotype file exists
  expect_true(file.exists(file.path(tmp, "phenotype.csv")))
  pheno <- utils::read.csv(file.path(tmp, "phenotype.csv"))
  expect_s3_class(pheno, "data.frame")
  expect_true(nrow(pheno) == 2)

  # FreeSurfer folder structure is correct
  expect_true(dir.exists(file.path(tmp, "cohort1", "sub-1_ses-01", "surf")))
  expect_true(dir.exists(file.path(tmp, "cohort1", "sub-2_ses-01", "surf")))
  # .mgh file exists
  expect_true(any(grepl("\\.mgh$", list.files(tmp, recursive = TRUE))))
})

test_that("simulate_dataset respects overwrite argument", {
  tmp <- file.path(tempdir(), "sym_dset_test2")

  simulate_dataset(
    path = tmp,
    data_structure = ds,
    overwrite = TRUE,
    seed = 123,
    verbose = verbosity
  )
  old_time <- file.info(file.path(tmp, "phenotype.csv"))$mtime
  Sys.sleep(1)

  simulate_dataset(
    path = tmp,
    data_structure = ds,
    overwrite = FALSE,
    seed = 123,
    verbose = verbosity
  )
  new_time <- file.info(file.path(tmp, "phenotype.csv"))$mtime
  expect_equal(old_time, new_time)
})

test_that("simulate_dataset handles simulate_association argument", {
  tmp <- tempdir()

  expect_silent(
    simulate_dataset(
      path = file.path(tmp, 'sym1'),
      data_structure = ds,
      fs_template = "fsmicro",
      simulate_association = ds_assoc_num,
      verbose = verbosity
    )
  )

  expect_silent(
    simulate_dataset(
      path = file.path(tmp, 'sym2'),
      data_structure = ds,
      fs_template = "fsmicro",
      simulate_association = ds_assoc_chr,
      verbose = verbosity
    )
  )

  expect_error(
    simulate_dataset(
      path = file.path(tmp, 'sym3'),
      data_structure = ds,
      fs_template = "fsmicro",
      simulate_association = list("bad"),
      verbose = FALSE
    )
  )

})
