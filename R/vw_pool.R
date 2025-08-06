#' @title
#' Pool \code{lme4::lmer} model output across imputed datasets
#'
#' @description
#' This function combines estimates, standard errors and p-values across imputed
#' datasets, for a single model (i.g., one vertex).
#' The function was largely taken from \code{mice::pool()} and
#' \code{mice::summary.mipo()} code. It averages the estimates of the
#' complete data model, and computes relevant statistics, following Rubin's
#' rules (Rubin, 1987, p. 76).
#'
# #' The degrees of freedom calculation for the pooled estimates uses the
# #' Barnard-Rubin adjustment for "small" samples (<= 100.000), see
# #' \code{\link{barnard.rubin}} helper.
#' P-values are estimated using the *t-as-z* approach at the moment. This is known
#' to the anti-conservative for small sample sizes. However we preferred a
#' relatively lenient (and computationally inexpensive) solution at this stage.
#' We will be addressing Type I error mores strictly at the cluster forming stage.
#'
#' The residuals of the model are currently simply averaged across imputed
#' datasets, for lack of a better idea of how to combine them.
#'
#' @param out_stats : Output of \code{\link{single_lmm}}, i.e.: a list with two
#' elements:
#' \enumerate{
#' \item \code{"stats"}: a dataframe with estimates, SEs and p-values for each
#' fixed effect term)
#' \item \code{"resid"}: a vector of residuals for the given model.
# #' \item \code{"df"}: residual degrees of freedom for the give model.
#' }
#' @param m : Number of imputed dataset (to avoid recomputing it)
#'
#' @note
#' Used inside \code{\link{run_vw_lmm}}.
#'
#' @return A list containing the pooled coefficients, SEs, t- and p- values.
#'
#' @importFrom rlang .data
#' @importFrom stats pt
#'
#' @author Serena Defina, 2024.
#'
#' @export
#'
#' @references
#' Rubin, D.B. (1987). \emph{Multiple Imputation for Nonresponse in Surveys}.
#' New York: John Wiley and Sons.
#'
vw_pool <- function(out_stats, m) {

  fails <- vapply(out_stats,
                  function(i) if (!is.null(i$error)) i$error else NA_character_,
                  character(1))
  num_failed <- sum(!is.na(fails))

  if (num_failed > 0) {
    unique_errors <- paste(num_failed, "of", length(out_stats), " imputations failed. Errors: ",
                           paste(unique(stats::na.omit(fails)), collapse = " | "))
    return(unique_errors)
  }

  if (m == 1) {
    # Just reformat output (no pooling needed)
    s <- out_stats[[1]]$stats
    stats <- list(
      coef = s$qhat,
      se = s$se,
      t = s$tval,
      p = s$pval,
      resid = as.vector(out_stats[[1]]$resid),
      warning = out_stats[[1]]$warning
    )
    return(stats)
  }

  # Extract warnings (if any)
  warnings <- vapply(out_stats,
                     function(i) if (!is.null(i$warning)) i$warning else NA_character_,
                     character(1))
  num_warned <- sum(!is.na(warnings))

  if (num_warned > 0) {
    warning_msg <- paste(num_warned, "of", length(out_stats), " imputations gave Warnings: ",
                         paste(unique(stats::na.omit(warnings)), collapse = " | "))
  } else {
    warning_msg <- NULL
  }

  # Extract estimates, standard errors and p-values
  model_output <- do.call(rbind, lapply(out_stats, `[[`, 1)) # stats data.frame

  # Residual degrees of freedom (assumed equal across imputations)
  # dfcom <- out_stats[[1]][[3]]
  # If sample is large enough, do not perform Barnard-Rubin adjustment
  # dfcom <- ifelse(dfcom > 1e+05, Inf, dfcom)

  # Pool model output
  pooled_stats <- model_output %>%
    dplyr::group_by(.data$term) %>%
    dplyr::summarize(
      # Number of imputations
      m = m,
      # Mean coefficient
      qbar = mean(.data$qhat),
      # Mean squared SE
      ubar = mean(.data$se^2),
      b = stats::var(.data$qhat),
      # Calculate the total variance
      t = .data$ubar + (1 + 1 / .data$m) * .data$b,
      # Proportion of total variance due to missingness
      lambda = (1 + 1 / .data$m) * .data$b / .data$t,
      # Model degrees of freedom
      # dfcom = dfcom,
      # Barnard-Rubin adjusted degrees of freedom
      # df = barnard.rubin(.data$lambda, .data$m, .data$dfcom),
      # Relative increase in variance due to non-response
      # riv = (1 + 1 / .data$m) * .data$b / .data$ubar,
      # Fraction of missing information
      # fmi = (.data$riv + 2 / (.data$df + 3)) / (.data$riv + 1),
    )

  # Pooled estimates
  coef = pooled_stats$qbar
  # Pooled standard errors
  se = sqrt(pooled_stats$t)
  # Pooled test statistics
  tval = coef / se
  # Pooled p-values
  pval <- 2 * (1 - pnorm(abs(tval))) # t-as-z approach
  # pval = 2 * pt(-abs(tval), df = pooled_stats$df)

  # Average residuals across imputed datasets
  resid <- as.matrix(colMeans(do.call(rbind, lapply(out_stats, `[[`, 2))))

  list(
    "coef" = coef,
    "se" = se,
    "t" = tval,
    "p" = pval,
    "resid" = resid,
    "warning" = warning_msg
  )
}

#' @title
#' Pooled degrees of freedom calculation
#'
#' @description
#' Edited from \code{mice::barnard.rubin}. This function computes the residual
#' degrees of freedom for hypothesis testing, as proposed by Barnard & Rubin (1999).
#' It is used by the \code{\link{vw_pool}} function.
#'
#' @param lambda : Proportion of total variance due to missingness
#' @param m : Number of imputed datasets
#' @param dfcom : Residual degrees of freedom, estimated using \code{stats::df.residual()}
#'
#' @return  Value for the residual degrees of freedom
#'
#' @references
#' Barnard, J. and Rubin, D.B. (1999). Small sample degrees of
#' freedom with multiple imputation. \emph{Biometrika}, 86, 948-955.
#'
barnard.rubin <- function(lambda, m, dfcom = Inf) {
  lambda[lambda < 1e-04] <- 1e-04
  dfold <- (m - 1) / lambda^2
  if (is.infinite(dfcom)) {
    df <- dfold
  } else {
    dfobs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lambda)
    df <- dfold * dfobs / (dfold + dfobs)
  }
  df
}
