library(lme4)

# Example data
set.seed(3108)

subj_dir <- test_path("fixtures","example_data")
pheno <- read.csv(file.path(subj_dir, "phenotype.csv"))

n_obs <- nrow(pheno)
y <- rnorm(n_obs)
formula <- vw_area ~ sex + age + (1 | id)

test_that("single_lmm returns correct structure", {
  result <- single_lmm(pheno, y, formula)
  expect_type(result, "list")
  expect_named(result, c("stats", "resid"))
  expect_s3_class(result$stats, "data.frame")
  expect_true("term" %in% names(result$stats))
  expect_true("qhat" %in% names(result$stats))
  expect_true("se" %in% names(result$stats))
  expect_true("pval" %in% names(result$stats))
  expect_length(result$resid, n_obs)
})

test_that("single_lmm works with model_template", {
  # Fit a template model
  pheno2 <- pheno
  pheno2[all.vars(formula)[1]] <- y
  template <- lmer(formula, data = pheno2)
  # Use the template for single_lmm
  result <- single_lmm(pheno, y, formula, model_template = template)
  expect_type(result, "list")
  expect_s3_class(result$stats, "data.frame")
})

test_that("single_lmm disables pvalues when requested", {
  result <- single_lmm(pheno, y, formula, pvalues = FALSE)
  expect_false("pval" %in% names(result$stats))
})

test_that("single_lmm errors with incorrect input", {
  bad_pheno <- list(a = 1, b = 2) # not a data.frame
  expect_error(single_lmm(bad_pheno, y, formula))
})

# ==============================================================================
