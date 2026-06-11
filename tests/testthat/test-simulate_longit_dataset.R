
ds_single <- list(
  cohort1 = list(sessions = c("01", "02"), n_subjects = 20)
)
ds_multi <- list(
  cohort1  = list(sessions = c("01", "02", "03"), n_subjects = 30),
  cohort2  = list(sessions = c("01", "02"),       n_subjects = 25)
)
bl <- list(
  age    = c(mean = 10, sd = 1),
  sex    = c(levels = c("Male", "Female")),
  wisdom = c(mean = 0, sd = 1)
)
ch <- list(
  age    = c(mean = 4, sd = 0.5),
  wisdom = c(mean = 1, sd = 0.5)
)

# ===========================================================================
# build_minimal_pheno
# ===========================================================================
test_that("build_minimal_pheno returns correct dimensions", {
  out <- build_minimal_pheno(ds_multi)
  # 3 sessions × 30 subj + 2 sessions × 25 subj = 140
  expect_equal(nrow(out), 3 * 30 + 2 * 25)
  expect_named(out, c("id", "site", "time", "folder_id"))
})

test_that("build_minimal_pheno produces correct folder_id format", {
  out <- build_minimal_pheno(ds_single)
  expect_true(all(grepl("^cohort1/sub-\\d+_ses-0[12]$", out$folder_id)))
})

test_that("build_minimal_pheno errors on duplicated sessions", {
  bad <- list(X = list(sessions = c("01", "01"), n_subjects = 5))
  expect_error(build_minimal_pheno(bad))
})

test_that("build_minimal_pheno errors on n_subjects <= 0", {
  bad <- list(X = list(sessions = "01", n_subjects = 0))
  expect_error(build_minimal_pheno(bad))
})

test_that("build_minimal_pheno errors on missing spec keys", {
  expect_error(build_minimal_pheno(list(X = list(sessions = "01"))))
  expect_error(build_minimal_pheno(list(X = list(n_subjects = 5))))
})

# ===========================================================================
# simulate_long_pheno_data
# ===========================================================================
test_that("simulate_long_pheno_data returns correct row count", {
  pheno <- simulate_long_pheno_data(ds_multi, bl, ch, verbose = FALSE)
  expect_equal(nrow(pheno), 3 * 30 + 2 * 25)
})

test_that("simulate_long_pheno_data always contains mandatory columns", {
  pheno <- simulate_long_pheno_data(ds_single, bl, ch, verbose = FALSE)
  expect_true(all(c("site", "id", "time", "folder_id") %in% names(pheno)))
})

test_that("simulate_long_pheno_data produces baseline covariate columns", {
  pheno <- simulate_long_pheno_data(ds_single, bl, ch, verbose = FALSE)
  expect_true(all(c("age", "sex", "wisdom") %in% names(pheno)))
})

test_that("sex is a factor constant within subject across sessions", {
  pheno <- simulate_long_pheno_data(ds_single, bl, ch, verbose = FALSE)
  expect_s3_class(pheno$sex, "factor")
  sex_per_subj <- tapply(as.character(pheno$sex), pheno$id, function(x) length(unique(x)))
  expect_true(all(sex_per_subj == 1))
})

test_that("age increases on average across sessions", {
  pheno <- simulate_long_pheno_data(ds_single, bl, ch, verbose = FALSE)
  mean_age <- tapply(pheno$age, pheno$time, mean)
  expect_true(mean_age["02"] > mean_age["01"])
})

test_that("simulate_long_pheno_data is reproducible with same seed", {
  p1 <- simulate_long_pheno_data(ds_single, bl, ch, seed = 42, verbose = FALSE)
  p2 <- simulate_long_pheno_data(ds_single, bl, ch, seed = 42, verbose = FALSE)
  expect_equal(p1, p2)
})

test_that("different seeds produce different data", {
  p1 <- simulate_long_pheno_data(ds_single, bl, ch, seed = 1, verbose = FALSE)
  p2 <- simulate_long_pheno_data(ds_single, bl, ch, seed = 2, verbose = FALSE)
  expect_false(identical(p1$age, p2$age))
})

test_that("single-session data has no follow-up rows", {
  ds_one <- list(X = list(sessions = "01", n_subjects = 10))
  pheno <- simulate_long_pheno_data(ds_one, bl, ch, verbose = FALSE)
  expect_equal(nrow(pheno), 10)
  expect_equal(unique(pheno$time), "01")
})

test_that("folder_id contains site, subject id, and session", {
  pheno <- simulate_long_pheno_data(ds_single, bl, ch, verbose = FALSE)
  expect_true(all(grepl("^cohort1/sub-\\d+_ses-0[12]$", pheno$folder_id)))
})

test_that("folder_id values are unique per observation", {
  pheno <- simulate_long_pheno_data(ds_multi, bl, ch, verbose = FALSE)
  expect_equal(nrow(pheno), length(unique(pheno$folder_id)))
})

# ===========================================================================
# validate_pheno
# ===========================================================================
test_that("validate_pheno builds minimal pheno when pheno=NULL and no assocs", {
  out <- validate_pheno(NULL, ds_single, list())
  expect_s3_class(out, "data.frame")
  expect_true("folder_id" %in% names(out))
})

test_that("validate_pheno errors when pheno=NULL but associations exist", {
  assoc <- list(temporalpole = c(age = 0.3))
  expect_error(validate_pheno(NULL, ds_single, assoc))
})

test_that("validate_pheno errors when required column is missing", {
  pheno <- data.frame(folder_id = "a/sub-1_ses-01")   # no 'age' column
  assoc <- list(temporalpole = c(age = 0.3))
  expect_error(validate_pheno(pheno, ds_single, assoc))
})

test_that("validate_pheno passes a valid pheno through unchanged", {
  pheno <- data.frame(folder_id = "a/sub-1_ses-01", age = 10)
  assoc <- list(temporalpole = c(age = 0.3))
  out <- validate_pheno(pheno, ds_single, assoc)
  expect_identical(out, pheno)
})

# ===========================================================================
# validate_roi_associations
# ===========================================================================
test_that("validate_roi_associations allows empty list when simulate_other_rois=TRUE", {
  expect_equal(validate_roi_associations(list(), simulate_other_rois = TRUE), list())
})

test_that("validate_roi_associations errors on empty list when simulate_other_rois=FALSE", {
  expect_error(validate_roi_associations(list(), simulate_other_rois = FALSE))
})

test_that("validate_roi_associations errors on duplicated ROI names", {
  assoc <- list(temporalpole = c(age = 0.3), temporalpole = c(age = 0.5))
  expect_error(validate_roi_associations(assoc, simulate_other_rois = FALSE))
})

test_that("validate_roi_associations errors on non-numeric association vector", {
  assoc <- list(temporalpole = c(age = "high"))
  expect_error(validate_roi_associations(assoc, simulate_other_rois = FALSE))
})

# ===========================================================================
# simulate_freesurfer_data
# ===========================================================================

# Helper pheno shared across FS tests
pheno_fs <- simulate_long_pheno_data(
  data_structure = list(GENR = list(sessions = c("01", "02"), n_subjects = 10)),
  baseline = list(age = c(mean = 10, sd = 1), sex = c(levels = c("Male", "Female"))),
  change   = list(age = c(mean = 4, sd = 0.5)),
  verbose  = FALSE
)

test_that("simulate_freesurfer_data writes one .mgh file per observation", {
  tmp <- withr::local_tempdir()
  simulate_freesurfer_data(
    path             = tmp,
    pheno            = pheno_fs,
    roi_associations = list(temporalpole = c(age = 0.3)),
    hemi             = "lh",
    verbose          = FALSE
  )
  mgh_files <- list.files(tmp, pattern = "\\.mgh$", recursive = TRUE)
  expect_equal(length(mgh_files), nrow(pheno_fs))
})

test_that("simulate_freesurfer_data writes files in the correct path layout", {
  tmp <- withr::local_tempdir()
  simulate_freesurfer_data(
    path             = tmp,
    pheno            = pheno_fs,
    roi_associations = list(temporalpole = c(age = 0.3)),
    hemi             = "lh",
    verbose          = FALSE
  )
  # Every folder_id should have a surf/ subdirectory with an .mgh file
  expected_dirs <- file.path(tmp, pheno_fs$folder_id, "surf")
  expect_true(all(dir.exists(expected_dirs)))
})

test_that("simulate_freesurfer_data .mgh filename matches hemi/measure/fwhmc/template", {
  tmp <- withr::local_tempdir()
  simulate_freesurfer_data(
    path             = tmp,
    pheno            = pheno_fs,
    roi_associations = list(temporalpole = c(age = 0.3)),
    hemi             = "lh",
    measure          = "thickness",
    fwhmc            = "fwhm10",
    fs_template      = "fsaverage",
    verbose          = FALSE
  )
  mgh_files <- list.files(tmp, pattern = "\\.mgh$", recursive = TRUE)
  expect_true(all(grepl("lh\\.thickness\\.fwhm10\\.fsaverage\\.mgh$", mgh_files)))
})

test_that("simulate_freesurfer_data returns NULL invisibly", {
  tmp <- withr::local_tempdir()
  out <- simulate_freesurfer_data(
    path             = tmp,
    pheno            = pheno_fs,
    roi_associations = list(temporalpole = c(age = 0.3)),
    hemi             = "lh",
    verbose          = FALSE
  )
  expect_null(out)
})

test_that("simulate_freesurfer_data is reproducible with same seed", {
  tmp1 <- withr::local_tempdir()
  tmp2 <- withr::local_tempdir()
  args <- list(
    pheno            = pheno_fs,
    roi_associations = list(temporalpole = c(age = 0.3)),
    hemi             = "lh",
    seed             = 42,
    verbose          = FALSE
  )
  do.call(simulate_freesurfer_data, c(list(path = tmp1), args))
  do.call(simulate_freesurfer_data, c(list(path = tmp2), args))

  files1 <- sort(list.files(tmp1, pattern = "\\.mgh$", recursive = TRUE))
  files2 <- sort(list.files(tmp2, pattern = "\\.mgh$", recursive = TRUE))
  expect_equal(
    readBin(file.path(tmp1, files1[[1]]), "raw", n = 1e5),
    readBin(file.path(tmp2, files2[[1]]), "raw", n = 1e5)
  )
})

test_that("simulate_freesurfer_data errors when simulate_other_rois=FALSE and roi_associations empty", {
  tmp <- withr::local_tempdir()
  expect_error(
    simulate_freesurfer_data(
      path                = tmp,
      pheno               = pheno_fs,
      roi_associations    = list(),
      simulate_other_rois = FALSE,
      hemi                = "lh",
      verbose             = FALSE
    )
  )
})

test_that("simulate_freesurfer_data errors on negative vw_sd", {
  tmp <- withr::local_tempdir()
  expect_error(
    simulate_freesurfer_data(
      path             = tmp,
      pheno            = pheno_fs,
      roi_associations = list(temporalpole = c(age = 0.3)),
      hemi             = "lh",
      vw_sd            = -0.1,
      verbose          = FALSE
    )
  )
})

test_that("simulate_freesurfer_data vertex values are positive (pmax floor applied)", {
  tmp <- withr::local_tempdir()
  simulate_freesurfer_data(
    path             = tmp,
    pheno            = pheno_fs,
    roi_associations = list(temporalpole = c(age = 0.3)),
    simulate_other_rois = TRUE,
    hemi             = "lh",
    fs_template = "fsaverage3",
    vw_mean          = 2.5,
    verbose          = FALSE
  )
  # Read back one file and check all active vertices are > 0
  first_file <- list.files(tmp, pattern = "\\.mgh$", recursive = TRUE, full.names = TRUE)[[1]]
  vals <- load.mgh(first_file)
  active_vals <- vals[vals != 0]
  expect_true(all(active_vals >= 0.001))
})

test_that("simulate_freesurfer_data writes rh files when hemi='rh'", {
  tmp <- withr::local_tempdir()
  simulate_freesurfer_data(
    path             = tmp,
    pheno            = pheno_fs,
    roi_associations = list(temporalpole = c(age = 0.3)),
    hemi             = "rh",
    verbose          = FALSE
  )
  mgh_files <- list.files(tmp, pattern = "\\.mgh$", recursive = TRUE)
  expect_true(all(grepl("^rh\\.", basename(mgh_files))))
})

test_that("simulate_freesurfer_data accepts pheno=NULL with empty roi_associations", {
  tmp <- withr::local_tempdir()
  expect_no_error(
    simulate_freesurfer_data(
      path                = tmp,
      pheno               = NULL,
      data_structure      = list(GENR = list(sessions = "01", n_subjects = 5)),
      roi_associations    = list(),
      simulate_other_rois = TRUE,
      hemi                = "lh",
      verbose             = FALSE
    )
  )
})

# ===========================================================================
# simulate_longit_dataset
# ===========================================================================

test_that("simulate_longit_dataset returns NULL invisibly", {
  tmp <- withr::local_tempdir()
  out <- simulate_longit_dataset(
    path             = tmp,
    data_structure   = list(GENR = list(sessions = c("01", "02"), n_subjects = 10)),
    roi_associations = list(temporalpole = c(age = 1.3)),
    hemi             = "lh",
    verbose          = FALSE
  )
  expect_null(out)
})

test_that("simulate_longit_dataset writes phenotype.csv", {
  tmp <- withr::local_tempdir()
  simulate_longit_dataset(
    path             = tmp,
    data_structure   = list(GENR = list(sessions = c("01", "02"), n_subjects = 10)),
    roi_associations = list(temporalpole = c(age = 1.3)),
    hemi             = "lh",
    verbose          = FALSE
  )
  expect_true(file.exists(file.path(tmp, "phenotype.csv")))
})

test_that("phenotype.csv has correct number of rows", {
  tmp <- withr::local_tempdir()
  ds  <- list(GENR = list(sessions = c("01", "02"), n_subjects = 10))
  simulate_longit_dataset(
    path             = tmp,
    data_structure   = ds,
    roi_associations = list(temporalpole = c(age = 1.3)),
    hemi             = "lh",
    verbose          = FALSE
  )
  pheno_out <- utils::read.csv(file.path(tmp, "phenotype.csv"))
  expect_equal(nrow(pheno_out), 10 * 2)   # 10 subjects × 2 sessions
})

test_that("simulate_longit_dataset writes .mgh files for both hemispheres", {
  tmp <- withr::local_tempdir()
  simulate_longit_dataset(
    path             = tmp,
    data_structure   = list(GENR = list(sessions = c("01", "02"), n_subjects = 5)),
    roi_associations = list(temporalpole = c(age = 1.3)),
    hemi             = "both",
    verbose          = FALSE
  )
  mgh_files <- list.files(tmp, pattern = "\\.mgh$", recursive = TRUE)
  lh_files  <- grep("^lh\\.", basename(mgh_files), value = TRUE)
  rh_files  <- grep("^rh\\.", basename(mgh_files), value = TRUE)
  expect_gt(length(lh_files), 0)
  expect_gt(length(rh_files), 0)
  expect_equal(length(lh_files), length(rh_files))
})

test_that("simulate_longit_dataset single hemi='lh' writes only lh files", {
  tmp <- withr::local_tempdir()
  simulate_longit_dataset(
    path             = tmp,
    data_structure   = list(GENR = list(sessions = c("01", "02"), n_subjects = 5)),
    roi_associations = list(temporalpole = c(age = 1.3)),
    hemi             = "lh",
    verbose          = FALSE
  )
  mgh_files <- list.files(tmp, pattern = "\\.mgh$", recursive = TRUE)
  expect_true(all(grepl("^lh\\.", basename(mgh_files))))
})

test_that("simulate_longit_dataset is reproducible with same seed", {
  tmp1 <- withr::local_tempdir()
  tmp2 <- withr::local_tempdir()
  args <- list(
    data_structure   = list(GENR = list(sessions = c("01", "02"), n_subjects = 5)),
    roi_associations = list(temporalpole = c(age = 1.3)),
    hemi             = "lh",
    seed             = 99,
    verbose          = FALSE
  )
  do.call(simulate_longit_dataset, c(list(path = tmp1), args))
  do.call(simulate_longit_dataset, c(list(path = tmp2), args))

  csv1 <- utils::read.csv(file.path(tmp1, "phenotype.csv"))
  csv2 <- utils::read.csv(file.path(tmp2, "phenotype.csv"))
  expect_equal(csv1, csv2)
})

test_that("simulate_longit_dataset creates path if it does not exist", {
  tmp  <- withr::local_tempdir()
  new_path <- file.path(tmp, "nested", "sim_out")
  expect_false(dir.exists(new_path))
  simulate_longit_dataset(
    path             = new_path,
    data_structure   = list(GENR = list(sessions = "01", n_subjects = 3)),
    roi_associations = list(),
    simulate_other_rois = TRUE,
    hemi             = "lh",
    verbose          = FALSE
  )
  expect_true(dir.exists(new_path))
})