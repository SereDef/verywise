
# Example data
set.seed(3108)

subj_dir <- test_path("fixtures", "fs7")
pheno <- read.csv(file.path(subj_dir, "phenotype.csv"))

n_obs <- nrow(pheno)
test_y <- rnorm(n_obs, mean = 6.5, sd = 0.2)
test_formula <- vw_area ~ sex + age + (1 | id)

test_that("single_lmm returns correct structure", {
  result <- single_lmm(imp = pheno,
                       y = test_y,
                       formula = test_formula)

  expect_type(result, "list")
  expect_named(result, c("stats", "resid", "warning"))
  expect_s3_class(result$stats, "data.frame")
  expect_true(all(c("term", "qhat", "se") %in% names(result$stats)))

  expect_length(result$resid, n_obs)

  expect_true(all(!is.na(result$stats$qhat)))
})

test_that("single_lmm works with model_template", {
  # Fit a template model
  pheno2 <- pheno
  pheno2[all.vars(test_formula)[1]] <- test_y
  template <- lmer(test_formula, data = pheno2)
  # Use the template for single_lmm
  result <- single_lmm(pheno, test_y, test_formula, model_template = template)

  expect_type(result, "list")
  expect_s3_class(result$stats, "data.frame")
  expect_true(all(c("term", "qhat", "se") %in% names(result$stats)))
})

test_that("single_lmm handles fitting errors gracefully", {
  # Make an invalid formula
  bad_formula <- vw_area ~ non_existent_var + (1 | id)
  result <- single_lmm(imp = pheno,
                       y = test_y,
                       formula = bad_formula)

  expect_true("error" %in% names(result))
  expect_type(result$error, "character")
})

test_that("single_lmm captures warnings", {
  # Produce a warning with perfectly collinear predictors
  pheno$age_copy <- pheno$age

  result <- single_lmm(imp = pheno,
                       y = test_y,
                       formula = vw_area ~ age + age_copy + (1 | id))
  # lmerControl(
  #   optCtrl = list(maxfun = 10),     # low max iterations to induce possible convergence warning
  #   checkConv = list(action = "warning"), # ensure convergence issues produce warnings
  #   check.nobs.vs.nRE = "warning",   # warn if more random effects than observations
  #   check.nlev.gtr.1 = "warning"     # warn if grouping factor has <=1 level
  # )
  expect_true(is.character(result$warning))
  expect_match(result$warning,
               "rank deficient so dropping 1 column / coefficient")
})
