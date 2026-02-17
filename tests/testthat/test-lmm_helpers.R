
# Example data
set.seed(3108)

subj_dir <- test_path("fixtures", "fs7")
pheno <- read.csv(file.path(subj_dir, "phenotype.csv"))

n_obs <- nrow(pheno)
test_y <- rnorm(n_obs, mean = 6.5, sd = 0.2)

pheno['vw_area'] <- test_y

test_that("single_lmm works with model_template", {

  test_formula <- vw_area ~ sex + age + (1 | id)

  test_template <- lmer(test_formula, data = pheno)

  # Use the template for single_lmm
  result <- single_lmm(pheno, test_y, 'vw_area', model_template = test_template)

  expect_type(result, "list")
  expect_named(result, c("stats", "resid", "model_fit", "warning"))
  expect_s3_class(result$stats, "data.frame")
  expect_true(all(c("term", "qhat", "se") %in% names(result$stats)))

  expect_type(result$model_fit, "double")
  expect_length(result$model_fit, 5)

  expect_length(result$resid, n_obs)

  expect_true(all(!is.na(result$stats$qhat)))
})

test_that("single_lmm captures warnings", {
  # Produce a warning with perfectly collinear predictors
  pheno$age_copy <- pheno$age

  test_formula <- vw_area ~ sex + age + age_copy + (1 | id)

  test_template <- lmer(test_formula, data = pheno)

  result <- single_lmm(imp = pheno, model_template = test_template,
                       y = test_y, y_name = 'vw_area')

  expect_true(is.character(result$warning))
  expect_match(result$warning[1],
               regexp = "rank deficient so dropping 1 column")
})
