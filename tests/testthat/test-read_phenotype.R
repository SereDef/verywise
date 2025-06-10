# Test load_pheno_file function
test_that("load_pheno_file works for supported extensions", {

  # Mocking file existence and loading functionality
  pheno_file1 <- test_path("fixtures","example_data",
                          "phenotype.csv")
  pheno_file2 <- test_path("fixtures","example_data",
                           "phenotype_mids.rds")

  data1 <- load_pheno_file(pheno_file1)
  data2 <- load_pheno_file(pheno_file2)

  expect_s3_class(data1, "data.frame")  # Check if it's a data frame
  expect_s3_class(data2, "mids")  # Check if it's a mids object

})

test_that("load_pheno_file throws error for unsupported extensions", {

  pheno_file <- test_path("fixtures","example_data",
                          "site1/sub-1_ses-01",
                          "surf/lh.area.fwhm10.fsaverage.mgh")

  expect_error(load_pheno_file(pheno_file), "Unsupported file extension")

})

test_that("load_pheno_file throws error if file is not found", {
  expect_error(load_pheno_file("non_existent_file.txt"), "file not found.")
})

# Test check_pheno_obj function
test_that("check_pheno_obj works when the object exists", {

  test_obj <- data.frame(A = 1, B = 2)
  assign("test_obj", test_obj, envir = .GlobalEnv)

  obj <- check_pheno_obj("test_obj")
  expect_equal(obj, test_obj)


})

test_that("check_pheno_obj throws error when object does not exist", {
  expect_error(check_pheno_obj("non_existent_obj"),
               "`non_existent_obj` not found in the global environment.")
})
