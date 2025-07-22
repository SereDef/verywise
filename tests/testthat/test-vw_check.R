good_formula <- as.formula("thickness ~ age + sex")

test_that("check_formula accepts valid formula and returns measure", {
  expect_equal(check_formula(good_formula), "thickness")
})

test_that("check_formula rejects invalid formula", {
  bad_formula <- "not_a_formula"
  expect_error(check_formula(bad_formula), "must be a valid R formula object")
})

test_that("check_formula rejects unknown measure", {
  bad_metric <- as.formula("charliexcx ~ age + sex")
  expect_error(check_formula(bad_metric),
               "The outcome in `formula` should be a brain surface metric.")
})

test_that("check_data_list accepts valid data_list", {
  df <- data.frame(folder_id = 1:3,
                   age = 20:22,
                   sex = c("M", "F", "M"))
  dl <- list(df)
  expect_null(check_data_list(dl, "folder_id", good_formula))
})

test_that("check_data_list errors on missing folder_id", {
  df <- data.frame(id = 1:3,
                   age = 20:22,
                   sex = c("M", "F", "M"))
  dl <- list(df)
  expect_error(check_data_list(dl, "folder_id", good_formula), "not found in data")
})

test_that("check_data_list errors on duplicate folder_id", {
  df <- data.frame(folder_id = c(1, 1, 2),
                   age = 20:22,
                   sex = c("M", "F", "M"))
  dl <- list(df)
  expect_error(check_data_list(dl, "folder_id", good_formula), "contains duplicates")
})

test_that("check_data_list errors on missing formula variable", {
  df <- data.frame(folder_id = 1:3,
                   age = 20:22)
  dl <- list(df)
  expect_error(check_data_list(dl, "folder_id", good_formula), "not present in the data")
})

test_that("check_path errors if directory does not exist", {
  expect_error(check_path("nonexistent_dir"), "does not exist")
})

test_that("check_path warns if file exists", {
  tmp <- tempfile()
  dir.create(tmp)
  file.create(file.path(tmp, "file.txt"))
  expect_warning(check_path(tmp, file_exists="file.txt"), "already exists")
})

test_that("check_cores returns valid core number", {
  expect_true(check_cores(1) >= 1)
})

test_that("check_cores errors on invalid core number", {
  expect_error(check_cores(0), "1 or higher")
})

test_that("check_numeric_param accepts valid numeric", {
  expect_null(check_numeric_param(5, "param", lower = 1, upper = 10))
})

test_that("check_numeric_param errors on non-numeric", {
  expect_error(check_numeric_param("a", "param"), "must be a single numeric value")
})

test_that("check_numeric_param errors on out-of-bounds", {
  expect_error(check_numeric_param(0, "param", lower = 1, upper = 10), "must be between")
})

test_that("check_numeric_param errors on non-integer when integer required", {
  expect_error(check_numeric_param(1.5, "param", integer = TRUE), "must be an integer")
})

# Note: check_freesurfer_setup modifies environment variables and calls system commands,
# so it's best tested in an integration
