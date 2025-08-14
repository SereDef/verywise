# Create mock out_stats for 2 imputations
out_stats <- list(
  list(
    stats = data.frame(term = "x",
                       qhat = 1,
                       se = 0.5,
                       tval = 1.5,
                       pval = 0.5),
    resid = c(1, 2, 3),
    warning = NULL,
    error = NULL
  ),
  list(
    stats = data.frame(term = "x",
                       qhat = 2,
                       se = 2,
                       tval = 1.5,
                       pval = 0.5),
    resid = c(2, 3, 4),
    warning = NULL,
    error = NULL
  )
)

out_stats_with_warning <- out_stats
out_stats_with_warning[[1]]$warning <- "Convergence warning"

out_stats_with_error <- out_stats
out_stats_with_error[[2]]$error <- "Singular fit"

test_that("vw_pool returns pooled stats for valid input", {
  result <- vw_pool(out_stats, m = 2)
  expect_named(result, c("coef", "se", "p", "resid", "warning"))
  expect_true(is.numeric(result$coef))
  expect_true(is.numeric(result$se))
  expect_true(is.numeric(result$p))
  expect_true(is.matrix(result$resid) || is.numeric(result$resid))
  expect_equal(result$warning, "")
})

test_that("vw_pool returns warning string if warnings present", {
  result <- vw_pool(out_stats_with_warning, m = 2)
  expect_true(is.character(result$warning))
  expect_match(result$warning, "Convergence warning")
})

test_that("vw_pool returns error info and empty output if errors present", {
  result <- vw_pool(out_stats_with_error, m = 2)
  expect_true(is.character(result))
  expect_match(result, "Singular fit")
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
