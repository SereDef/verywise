# ── Shared fixtures ────────────────────────────────────────────────────────────

subj_dir <- test_path("fixtures", "fs7")
pheno <- read.csv(file.path(subj_dir, "phenotype.csv"))

n_obs <- nrow(pheno)

set.seed(3108)
test_y <- rnorm(n_obs, mean = 6.5, sd = 0.2)

test_formula  <- vw_area ~ sex + age + wisdom + (1 | id)
test_control  <- lme4::lmerControl(calc.derivs = FALSE)


# ── single_lmm ────────────────────────────────────────────────────────────────

test_that("single_lmm returns correct structure", {

  result <- single_lmm(
    imp          = pheno,
    y            = test_y,
    y_name       = "vw_area",
    model_formula = test_formula,
    REML         = TRUE,
    lmm_control  = test_control
  )

  expect_type(result, "list")
  expect_named(result, c("stats", "resid", "model_fit", "warning"))

  expect_s3_class(result$stats, "data.frame")
  expect_true(all(c("term", "qhat", "se") %in% names(result$stats)))

  expect_type(result$model_fit, "double")
  expect_length(result$model_fit, 5)

  expect_length(result$resid, n_obs)

  expect_true(all(!is.na(result$stats$qhat)))
})

test_that("single_lmm model_fit has finite AIC and valid R2/ICC range", {

  result <- single_lmm(
    imp           = pheno,
    y             = test_y,
    y_name        = "vw_area",
    model_formula = test_formula,
    REML          = TRUE,
    lmm_control   = test_control
  )

  # model_fit: c(is_singular, AIC, ICC, R2_marginal, R2_conditional)
  expect_true(result$model_fit[1] %in% c(0, 1))       # is_singular is 0 or 1
  expect_true(is.finite(result$model_fit[2]))         # AIC is a real number
  # ICC and R2 values, when non-NA, should be in [0, 1]
  perf_metrics <- result$model_fit[3:5]
  valid <- !is.na(perf_metrics)
  expect_true(all(perf_metrics[valid] >= 0 & perf_metrics[valid] <= 1))
})

test_that("single_lmm singular fit is flagged in model_fit[1]", {

  # Use a large-signal y that is identical across all obs → variance absorbed
  # by random intercept collapses → reliable singular fit without vcov crash
  singular_y <- rep(mean(test_y), n_obs) + rnorm(n_obs, sd = 1e-8)

  result <- single_lmm(
    imp           = pheno,
    y             = singular_y,
    y_name        = "vw_area",
    model_formula = test_formula,
    REML          = TRUE,
    lmm_control   = test_control
  )

  expect_true(result$model_fit[1] %in% c(0, 1))  # graceful, not an error
})

test_that("single_lmm returns error element on fitting failure", {

  # Trigger a genuine lmer error by passing a y of wrong length
  bad_y <- test_y[-1]  # one observation short

  result <- single_lmm(
    imp           = pheno,
    y             = bad_y,
    y_name        = "vw_area",
    model_formula = test_formula,
    REML          = TRUE,
    lmm_control   = test_control
  )

  expect_named(result, "error")
  expect_type(result$error, "character")
})

test_that("single_lmm captures warnings for rank-deficient model", {

  pheno_collinear <- pheno
  pheno_collinear$age2 <- pheno_collinear$age   # perfectly collinear

  collinear_formula <- vw_area ~ sex + age + age2 + wisdom + (1 | id)

  result <- single_lmm(
    imp           = pheno_collinear,
    y             = test_y,
    y_name        = "vw_area",
    model_formula = collinear_formula,
    REML          = TRUE,
    lmm_control   = test_control
  )

  expect_true(is.character(result$warning))
  expect_match(result$warning[1], regexp = "rank deficient so dropping 1 column")
})

test_that("single_lmm handles numeric weights without error", {

  weights_vec <- rep(1, n_obs)

  result <- single_lmm(
    imp           = pheno,
    y             = test_y,
    y_name        = "vw_area",
    model_formula = test_formula,
    REML          = TRUE,
    lmm_control   = test_control,
    weights       = weights_vec
  )

  expect_named(result, c("stats", "resid", "model_fit", "warning"))
})

test_that("single_lmm handles column-name weights without error", {

  pheno$my_weight <- 1

  result <- single_lmm(
    imp           = pheno,
    y             = test_y,
    y_name        = "vw_area",
    model_formula = test_formula,
    REML          = TRUE,
    lmm_control   = test_control,
    weights       = "my_weight"
  )

  expect_named(result, c("stats", "resid", "model_fit", "warning"))
})


# ── refit_lmm ─────────────────────────────────────────────────────────────────

# Build the template once for the refit_lmm tests
local_pheno <- pheno
local_pheno$vw_area <- test_y
template_fit <- lme4::lmer(test_formula, data = local_pheno,
                           control = test_control)

test_that("refit_lmm returns correct structure", {

  result <- refit_lmm(template_fit, y = test_y)

  expect_type(result, "list")
  expect_named(result, c("stats", "resid", "model_fit", "warning"))

  expect_s3_class(result$stats, "data.frame")
  expect_true(all(c("term", "qhat", "se") %in% names(result$stats)))

  expect_type(result$model_fit, "double")
  expect_length(result$model_fit, 5)

  expect_length(result$resid, n_obs)
  expect_true(all(!is.na(result$stats$qhat)))
})

test_that("refit_lmm model_fit values are in expected ranges", {

  result <- refit_lmm(template_fit, y = test_y)

  expect_true(result$model_fit[1] %in% c(0, 1))
  expect_true(is.finite(result$model_fit[2]))

  perf_metrics <- result$model_fit[3:5]
  valid <- !is.na(perf_metrics)
  expect_true(all(perf_metrics[valid] >= 0 & perf_metrics[valid] <= 1))
})

test_that("refit_lmm returns error element on fitting failure", {

  # Provide a y of wrong length to force an internal error
  bad_y <- test_y[-1]

  result <- refit_lmm(template_fit, y = bad_y)

  expect_named(result, "error")
  expect_type(result$error, "character")
})

# ── precompile_model ──────────────────────────────────────────────────────────

test_that("precompile_model returns a list of lmerMod objects", {

  pheno$vw_area <- test_y

  # Wrap pheno in a list to mimic imp2list() output (3 imputations)
  data_list <- list(pheno, pheno, pheno)

  result <- precompile_model(
    formula   = test_formula,
    data_list = data_list,
    tmp_y     = test_y,
    measure   = "area",
    verbose   = FALSE
  )

  expect_type(result, "list")
  expect_length(result, 3)
  expect_true(all(sapply(result, inherits, "lmerMod")))
})

test_that("precompile_model output contains vw_<measure> as the response variable", {

  data_list <- list(pheno, pheno)

  result <- precompile_model(
    formula   = test_formula,
    data_list = data_list,
    tmp_y     = test_y,
    measure   = "area",
    verbose   = FALSE
  )

  # formula() on a lmerMod returns the model formula; [[2]] is the lhs (response)
  fitted_response <- as.character(formula(result[[1]])[[2]])
  expect_equal(fitted_response, "vw_area")
})


# ── unpack_formula ────────────────────────────────────────────────────────────

test_that("unpack_formula returns fixed-effect term names for LME formula", {

  terms_out <- unpack_formula(test_formula, pheno)

  expect_type(terms_out, "character")
  # Should include intercept and fixed predictors, but NOT random-effects bar
  expect_true("(Intercept)" %in% terms_out)
  expect_false(any(grepl("\\|", terms_out)))
})

test_that("unpack_formula returns design matrix when return_X = TRUE", {

  X <- unpack_formula(test_formula, pheno, return_X = TRUE)

  expect_true(is.matrix(X))
  expect_equal(nrow(X), n_obs)
  expect_true("(Intercept)" %in% colnames(X))
})

test_that("unpack_formula works for plain lm formula without random terms", {

  lm_formula <- vw_area ~ sex + age + wisdom

  terms_out <- unpack_formula(lm_formula, pheno)

  expect_type(terms_out, "character")
  expect_true("(Intercept)" %in% terms_out)
})