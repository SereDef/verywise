good_formula <- as.formula("vw_thickness ~ age + sex")

test_that("check_formula works with valid and invalid input", {

  bad_formula <- "not_a_formula"
  bad_metric <- as.formula("charliexcx ~ age + sex")

  # check_formula accepts valid formula and returns measure
  expect_equal(check_formula(good_formula), "thickness")
  # check_formula rejects invalid formula
  expect_error(check_formula(bad_formula),
               regexp = "`formula` must be a valid R formula object")
  # check_formula rejects unknown measures
  expect_error(check_formula(bad_metric),
               regexp = "should be a brain surface metric.")
})

test_that("check_data_list validates structure and contents", {
  df <- data.frame(folder_id = 1:3,
                   age = 20:22,
                   sex = c("M", "F", "M"))
  data_list <- list(df)

  # check_data_list accepts valid data_list
  expect_null(check_data_list(data_list, "folder_id", good_formula))

  # check_data_list errors when data_list is malformed
  bad_data_list <- list(df, df[, -1])
  expect_error(check_data_list(bad_data_list, "folder_id", good_formula),
               regexp = "must have the same dimentions")

  # check_data_list errors on missing folder_id
  df_wrongID <- df; names(df_wrongID)[1] <- 'id'
  bad_data_list <- list(df_wrongID)
  expect_error(check_data_list(bad_data_list, "folder_id", good_formula),
               regexp = "not found in data")

  # check_data_list errors when there are NAs in folder_id
  df_nainID <- df; df_nainID$folder_id[1] <- NA
  bad_data_list <- list(df_nainID)
  expect_error(check_data_list(bad_data_list, "folder_id", good_formula),
               regexp = "contains missing values")

  # check_data_list errors on duplicate folder_id
  df_dupID <- df; df_dupID$folder_id[2] <- df_dupID$folder_id[1]
  bad_data_list <- list(df_dupID)
  expect_error(check_data_list(bad_data_list, "folder_id", good_formula),
               regexp = "contains duplicates")

  # check_data_list errors on missing formula variable
  df_nosex <- df[,c('folder_id', 'age')]
  bad_data_list <- list(df_nosex)
  expect_error(check_data_list(bad_data_list, "folder_id", good_formula),
               regexp = "not present in the data")
})

test_that("check_path creates or checks path correctly", {
  tmp <- tempfile()

  # check_path errors if directory does not exist
  expect_error(check_path(tmp),
               regexp = "does not exist")
  # check_path warns and creates the directory when create_if_not = TRUE
  expect_message(check_path(tmp, create_if_not = TRUE),
                 regexp = "does not exist.*create it")
  expect_true(dir.exists(tmp))

  # check_path returns NULL if file does not exist
  expect_null(check_path(tmp, file_exists = "not_a_file.txt"))
  # check_path returns the file path if file exists
  file.create(file.path(tmp, "a_file.txt"))
  expect_equal(check_path(tmp, file_exists = "a_file.txt"),
               file.path(tmp, "a_file.txt"))
})

test_that("check_stack_file writes or validates stack_names.txt", {
  fixed_terms <- c("age", "sex")

  tmp <- tempfile()
  dir.create(tmp)

  expect_silent(check_stack_file(fixed_terms, tmp))

  # Tamper file
  writeLines("corrupted", file.path(tmp, "stack_names.txt"))
  expect_error(check_stack_file(fixed_terms, tmp),
               regexp = "file already exists.*different than expected")
})

test_that("check_weights validates weight input correctly", {
  df <- data.frame(folder_id = 1:4,
                   age = 20:23,
                   weight_var = c(1, 1.5, 2, 1))

  good_weight_vector <- c(1, 2, 1, 2)
  short_weight_vector <- c(1, 2)
  bad_weight_format1 <- c('1', '2', '1', '2')
  bad_weight_format2 <- list(1, 2, 1, 2)

  # check_weights works when weights are not provided
  expect_silent(check_weights(NULL, df))
  # check_weights works with either a variable in data or a numeric vector
  expect_silent(check_weights("weight_var", df))
  expect_silent(check_weights(good_weight_vector, df))

  # check_weights errors if weight variable is not in data
  expect_error(check_weights("missing_weight", df),
               regexp = "Weights variable 'missing_weight' not found")
  # check_weights errors if weight vector does not have the correct dimensions
  expect_error(check_weights(short_weight_vector, df),
               regexp = "not the same dimensions")
  # check_weights errors if weight vector is not numeric
  expect_error(check_weights(bad_weight_format1, df),
               regexp = "The weights format you provided is not valid")
  # check_weights errors if weight vector is not a vector
  expect_error(check_weights(bad_weight_format2, df),
               regexp = "The weights format you provided is not valid")
})

test_that("check_cores adjusts the number of cores and warns", {
  skip_on_cran()  # system-specific

  # check_cores returns valid core number
  expect_true(check_cores(1) >= 1)
  expect_type(check_cores(1), "integer")
  # handles non-numeric inputs as well
  expect_type(check_cores(1.1), "integer")
  expect_type(check_cores("1"), "integer")

  # check_cores errors on invalid core number
  expect_error(check_cores(0),
               regexp = "should be an integer that is 1 or higher")
})

test_that("check_numeric_param validates numeric parameters", {
  # check_numeric_param accepts valid numeric
  expect_silent(check_numeric_param(5, "param", lower = 1, upper = 10))

  # check_numeric_param errors on non-numeric
  expect_error(check_numeric_param("a", "param"),
               regexp = "`param` must be a single numeric value")

  # check_numeric_param errors on non-integer when integer required
  expect_error(check_numeric_param(1.5, "param", integer = TRUE),
               regexp = "`param` must be an integer")

  # check_numeric_param errors on out-of-bounds
  expect_error(check_numeric_param(0, "param", lower = 1, upper = 10),
               regexp = "must be between 1 and 10")

})

# Note: check_freesurfer_setup modifies environment variables and calls system
# commands, so it's best tested in an integration
