# Create mock out_stats for 2 imputations
out_stats <- list(
  list(
    stats = data.frame(term = c("(intercept)","x1","x2"),
                       qhat = c(1,2,3),
                       se = rep(0.5,3)),
    resid = c(1, 2, 3),
    model_fit = c(0, 1.5, 0, 0.1, 0.1),
    warning = character(0)
  ),
  list(
    stats = data.frame(term = c("(intercept)","x1","x2"),
                       qhat = c(0.9, 2.1, 3.4),
                       se = rep(0.6,3)),
    resid = c(2, 3, 4),
    model_fit = c(1, 3, 0.1, 0.4, 0.5),
    warning = character(0)
  )
)

out_stats_with_warning <- out_stats
out_stats_with_warning[[1]]$warning <- "Convergence warning"

out_stats_with_error <- out_stats
out_stats_with_error[[3]] <- list(error = 'You suck')

test_that("vw_pool returns pooled stats for valid input", {
  result <- vw_pool(out_stats, m = 2, n_terms = 3)
  expect_named(result, c('coef', 'se', 'p', 'resid', 'mfit', 'cov', 'warning'))
  expect_true(is.numeric(result$coef))
  expect_true(is.numeric(result$se))
  expect_true(is.numeric(result$p))
  expect_true(is.numeric(result$mfit))
  expect_true(is.matrix(result$resid) || is.numeric(result$resid))
  expect_equal(result$warning, "")
})

test_that("vw_pool returns warning string if warnings present", {
  result <- vw_pool(out_stats_with_warning, m = 2, n_terms = 3)
  expect_true(is.character(result$warning))
  expect_match(result$warning, "Convergence warning")
})

test_that("vw_pool returns error info and empty output if errors present", {
  result <- vw_pool(out_stats_with_error, m = 3, n_terms = 3)
  expect_true(is.character(result))
  expect_match(result, "1 / 3 imputations failed. Errors: You suck")
})

test_that("barnard.rubin returns finite degrees of freedom", {
  df <- barnard.rubin(lambda = 0.5, m = 5, dfcom = 100)
  expect_true(is.numeric(df))
  expect_true(df > 0)
})

test_that("barnard.rubin handles infinite dfcom", {
  df <- barnard.rubin(lambda = 0.5, m = 5, dfcom = Inf)
  expect_true(is.numeric(df))
})

# ── vw_pool() covariance pooling (cov_eff) ────────────────────────────────────

make_mock_out_stats_with_cov <- function(cov_vals) {
  list(
    list(
      stats = data.frame(term = c("(intercept)", "x1", "x2"),
                          qhat = c(1, 2, 3), se = rep(0.5, 3)),
      cov = cov_vals[1],
      resid = c(1, 2, 3),
      model_fit = c(0, 1.5, 0, 0.1, 0.1),
      warning = character(0)
    ),
    list(
      stats = data.frame(term = c("(intercept)", "x1", "x2"),
                          qhat = c(0.9, 2.1, 3.4), se = rep(0.6, 3)),
      cov = cov_vals[2],
      resid = c(2, 3, 4),
      model_fit = c(1, 3, 0.1, 0.4, 0.5),
      warning = character(0)
    )
  )
}

test_that("vw_pool passes covariance through unchanged when m = 1 (no pooling)", {

  single_out_stats <- list(
    list(
      stats = data.frame(term = c("(intercept)", "x1", "x2"),
                          qhat = c(1, 2, 3), se = rep(0.5, 3)),
      cov = 0.042,
      resid = c(1, 2, 3),
      model_fit = c(0, 1.5, 0, 0.1, 0.1),
      warning = character(0)
    )
  )

  result <- vw_pool(single_out_stats, m = 1, n_terms = 3, cov_eff = c(2, 3))

  expect_equal(result$cov, 0.042)
})

test_that("vw_pool pools covariance across imputations using Rubin's rules", {

  cov_vals <- c(0.04, 0.06)
  out_stats_cov <- make_mock_out_stats_with_cov(cov_vals)

  result <- vw_pool(out_stats_cov, m = 2, n_terms = 3, cov_eff = c(2, 3))

  # Manually reproduce the pooling formula used inside vw_pool()
  qbar <- matrix(do.call(rbind, lapply(out_stats_cov, `[[`, "stats"))$qhat, nrow = 3)
  ubar_cov <- mean(cov_vals)
  b_cov <- stats::cov(qbar[2, ], qbar[3, ])
  expected_cov <- ubar_cov + (1 + 1 / 2) * b_cov

  expect_equal(result$cov, expected_cov, tolerance = 1e-10)
})

test_that("vw_pool returns NULL covariance when cov_eff is not supplied, even with valid cov data", {

  out_stats_cov <- make_mock_out_stats_with_cov(c(0.04, 0.06))

  result <- vw_pool(out_stats_cov, m = 2, n_terms = 3, cov_eff = NULL)

  expect_null(result$cov)
  # Other pooled quantities should be unaffected by the missing cov_eff
  expect_true(is.numeric(result$coef))
  expect_true(is.numeric(result$se))
})

test_that("vw_pool covariance pooling handles NA covariance values gracefully", {

  out_stats_cov <- make_mock_out_stats_with_cov(c(NA_real_, 0.06))

  result <- vw_pool(out_stats_cov, m = 2, n_terms = 3, cov_eff = c(2, 3))

  # ubar_cov uses na.rm = TRUE internally, so this should not be NA/NaN
  expect_false(is.na(result$cov))
})
