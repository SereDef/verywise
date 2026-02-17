# Create mock out_stats for 2 imputations
out_stats <- list(
  list(
    stats = data.frame(term = "x",
                       qhat = 1,
                       se = 0.5),
    resid = c(1, 2, 3),
    model_fit = c(0, 1.5, 0, 0.1, 0.1),
    warning = character(0)
  ),
  list(
    stats = data.frame(term = "x",
                       qhat = 2,
                       se = 2),
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
  result <- vw_pool(out_stats, m = 2)
  expect_named(result, c("coef", "se", "p", "fitstats", "resid",  "warning"))
  expect_true(is.numeric(result$coef))
  expect_true(is.numeric(result$se))
  expect_true(is.numeric(result$p))
  expect_true(is.numeric(result$fitstats))
  expect_true(is.matrix(result$resid) || is.numeric(result$resid))
  expect_equal(result$warning, "")
})

test_that("vw_pool returns warning string if warnings present", {
  result <- vw_pool(out_stats_with_warning, m = 2)
  expect_true(is.character(result$warning))
  expect_match(result$warning, "Convergence warning")
})

test_that("vw_pool returns error info and empty output if errors present", {
  result <- vw_pool(out_stats_with_error, m = 3)
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
