test_that("imp2list works correctly on mids files", {

  pheno_file <- test_path("fixtures","example_data",
                           "phenotype_mids.rds")

  mids_data <- load_pheno_file(pheno_file)

  data_list <- imp2list(mids_data)

  expect_type(data_list, "list")  # Check if it's a data frame
  expect_length(data_list, 3) # object contains 3 imputed datasets
  expect_s3_class(data_list[[1]], "data.frame") # each element is a data.frame

})
