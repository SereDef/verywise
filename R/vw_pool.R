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
#' \strong{P-values estimation}
#' Default:
#' P-values are estimated using the *t-as-z* approach. This is known
#' to the anti-conservative for small sample sizes but provides a
#' computationally efficient solution. Type I error control is addressed
#' more rigorously at the cluster-forming stage.
#' Wald Chi-square test
#' Satterthwaite Approximation
#'
#' The residuals of the model are currently simply averaged across imputed
#' datasets, for lack of a better idea of how to combine them.
#'
#' @param out_stats Output of \code{\link{single_lmm}}, i.e.: a list containing:
#'  \itemize{
#'    \item \code{"stats"}: a dataframe with estimates and SEs for each fixed
#'                          effect term)
#'    \item \code{"resid"}: a vector of residuals for the given model.
#'    \item \code{"warning"}: a character vector with warning messages (if any)
#'  }
#'  Note that if a model failed for that dataset \code{out_stats} has form
#'  \code{list("error"="Error message")}
#' @param m Integer indicating the number of (imputed) dataset
#' @param pvalue_method String indicating which approximation methods should be
#'   used to compute p-values. Options: \code{"t-as-z"}, \code{"wald-chi2"},
#'   \code{"satterthwaite"} (see Details). Default: \code{"t-as-z"}.
#' @param min_pvalue Float, used to avoid pvalues == 0L for which log10 is Inf.
#'   Set this to 0L if no trimming should be applied.
#'   Default: \code{2^-149} ( == 1.401298e-45) the smallest positive subnormal
#'   float.
#' @note
#' Used inside \code{\link{run_vw_lmm}}.
#'
#' @return A list containing the pooled coefficients, SEs, t- and p- values.
#'
#' @importFrom rlang .data
#' @importFrom stats pnorm
#'
#' @author Serena Defina, 2024.
#'
#' @export
#'
#' @references
#' Rubin, D.B. (1987). \emph{Multiple Imputation for Nonresponse in Surveys}.
#' New York: John Wiley and Sons.
#'
vw_pool <- function(out_stats, m,
                    pvalue_method = "t-as-z",
                    min_pvalue =  2^-149
                    ) {

  fails <- vapply(out_stats,
                  function(i) if (!is.null(i$error)) i$error else NA_character_,
                  character(1))
  num_failed <- sum(!is.na(fails))

  if (num_failed > 0) {
    unique_errors <- paste(num_failed, " of ", m, " imputations failed. Errors: ",
                           paste(unique(stats::na.omit(fails)), collapse = " || "))
    return(unique_errors)
  }

  #

  if (m == 1) {
    # No pooling needed, just reformat output (no pooling needed)
    s <- out_stats[[1]]$stats
    tval <- s$qhat / s$se
    # TODO: # https://www.r-bloggers.com/2014/02/three-ways-to-get-parameter-specific-p-values-from-lmer/
    if (pvalue_method == "t-as-z") {
      # t-as-z method
      # normal approximation
      s$pval <- 2 * (1 - pnorm(abs(tval)))
      s$pval[s$pval < min_pvalue] <- min_pvalue
      # t-distributed (but df are not defined)
      # s$pval <- 2 * (1 - stats::pt(abs(tval), df = resid_df))

    # } else if (pvalue_method == "wald-chi2") {
    #   # Wald chi-square tests
    #   # anova_tab <- Anova(fit, type = 3) # gives one pvalues per term
    #   #
    #   ## Wald chi-square per coefficient (1 df)
    #   chisq <- tval^2
    #   outp_stats$pval <- stats::pchisq(chisq, df = 1, lower.tail = FALSE)
    } else {
      warning("The specified p-value method is currently not supported, a t-as-z approach was used.")
      s$pval <- 2 * (1 - pnorm(abs(tval)))
      s$pval[s$pval < min_pvalue] <- min_pvalue
    }

    stats <- list(
      coef = s$qhat,
      se = s$se, # t = s$tval,
      p = s$pval,
      resid = as.vector(out_stats[[1]]$resid),
      warning = paste(out_stats[[1]]$warning, collapse = " || ")
    )
    return(stats)
  }

  # Extract warnings (if any)
  warnings <- lapply(out_stats, `[[`, "warning")

  warning_count <- sum(lengths(warnings) > 0)

  if (warning_count > 0) {
    # Replace only floating point numbers with <NUM>
    normalized_warnings <- gsub("\\b\\d+\\.\\d+(e[+-]?\\d+)?\\b", "<NUM>",
                                unlist(warnings, use.names = FALSE))
    warning_msg <- paste(warning_count, " of ", m, " imputations gave Warnings: ",
                         paste(unique(normalized_warnings), collapse = " || "))
  } else {
    warning_msg <- ""
  }

  # Extract estimates and standard errors (i.e. the "stats" data.frame)
  model_output <- do.call(rbind, lapply(out_stats, `[[`, "stats"))

  # Residual degrees of freedom (assumed equal across imputations)
  # dfcom <- out_stats[[1]][["df"]]
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
      # lambda = (1 + 1 / .data$m) * .data$b / .data$t,
      # Model degrees of freedom
      # dfcom = dfcom,
      # Barnard-Rubin adjusted degrees of freedom
      # df = barnard.rubin(.data$lambda, .data$m, dfcom = Inf), #.data$dfcom),
      # Df with Rubin's rules
      # df <- (m - 1) * (1 + (u_bar / ((1 + 1/m) * b)))^2
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
  pval[pval < min_pvalue] <- min_pvalue

  # pval = 2 * pt(-abs(tval), df = pooled_stats$df)

  # Average residuals across imputed datasets
  resid <- as.matrix(colMeans(do.call(rbind, lapply(out_stats, `[[`, 2))))

  list(
    "coef" = coef,
    "se" = se, # "t" = tval,
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
