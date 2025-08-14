#' @title
#' Run a single linear mixed model and extract statistics
#'
#' @description
#' Fits a linear mixed model to a single vertex outcome using
#' \code{lme4::lmer()} and extracts fixed effects statistics. This function
#' is called repeatedly across all cortical vertices during vertex-wise
#' analysis.
#'
#' @param imp A data.frame containing the phenotype dataset
#'   (in verywise format).
#' @param y A numeric vector of outcome values representing a single
#'   vertex from the super-subject matrix.
#' @param formula An R formula object describing the linear mixed model using
#'   \code{lme4} notation.
#' @param model_template Optional pre-compiled model object for faster
#'   estimation. When provided, \code{single_lmm} will use an "update"-based
#'   workflow instead of refitting the model from scratch. This minimizes
#'   repeated parsing and model construction overhead, significantly reducing
#'   computation time for large-scale vertex-wise analyses.
#'   Default: \code{NULL}.
#' @inheritParams run_vw_lmm
#'
#' @details
#' Additional parameters are currently passed to the \code{lme4::lmer} call using
#' the \code{lmm_control} argument.
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item \code{stats} - A data.frame with columns:
#'       \itemize{
#'         \item \code{term} - Fixed effect term names
#'         \item \code{qhat} - Parameter estimates
#'         \item \code{se} - Standard errors
#'       }
#'     \item \code{resid} - A numeric vector of model residuals
#'     \item \code{warning} - Warning message(s) if any
#'   }
#'
#' @seealso \code{\link{run_vw_lmm}} for the main interface.
#'
#' @importFrom lme4 lmer
#' @importFrom lme4 fixef
#' @importFrom stats residuals
#'
#' @export
#'
single_lmm <- function(
    imp, y, formula,
    model_template = NULL,
    weights = NULL,
    lmm_control = lme4::lmerControl())
{

  if (!is.null(weights) && !is.numeric(weights)) {
    weights <- imp[, weights]
  }

  # Add (vertex) outcome to (single) dataset
  imp[all.vars(formula)[1]] <- y


  error_msg <- NULL
  warning_msg <- character(0)
  message_msg <- character(0)

  # Fit linear mixed model using `lme4::lmer`
  fit <- withCallingHandlers(
    tryCatch(
      {
        if (!is.null(model_template)) {
          stats::update(model_template, data = imp, weights = weights)
        } else {
          lmer(formula = formula, data = imp, weights = weights,
               control = lmm_control)
        }
      },
      error = function(e) {
        error_msg <<- conditionMessage(e)
        NULL
      }
    ),
    warning = function(w) {
      warning_msg <<- c(warning_msg, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      message_msg <<- c(message_msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  if (!is.null(error_msg)) {
    result <- list("error" = error_msg)
    return(result)
  }

  # coef(fit)$id # A matrix of effects by random variable
  # lme4::ranef(fit) # extract random effects (these should sum to 0)
  # lme4::VarCorr(fit) # estimated variances, SDs, and correlations between the random-effects terms

  # Extract estimates, standard errors and p-values
  outp_stats <- data.frame(
    "qhat" = as.matrix(fixef(fit)), # Fixed effects estimates
    "se" = as.matrix(summary(fit)$coefficients[, "Std. Error"]) # standard error
  )

  # TODO: implement "stack of interest"?

  # Extract residual degrees of freedom for Barnard-Rubin adjustment
  # Normally, this would be the number of independent observation minus the
  # number of fitted parameters, but not exactly what is done here.
  # Following the `broom.mixed` package approach, which `mice::pool` relies on
  # resid_df <- df.residual(fit)

  # Save row.names (i.e. terms) as a column so these can be grouped later
  outp_stats <- data.frame("term" = row.names(outp_stats), outp_stats,
                           row.names = NULL)

  # Also extract model residuals for smoothness estimation
  resid <- residuals(fit)

  result <- list("stats" = outp_stats,
                 "resid" = resid,
                 "warning" = c(warning_msg, message_msg))

  return(result)
}


#' @title
#' Unpack \code{lme4} formula
#'
#' @description
#' Get the terms for the fixed effect for the given formula.
#'
#' @param formula : model formula object (this should specify a LME model)
#' @param dset : the data.frame, for example the first element in the list
#' output of \code{\link{imp2list}}.
#'
#' @return A character vector of fixed terms.
#'
#' @author Serena Defina, 2024.
#'
get_terms <- function(formula, dset) {
  # Add a placeholder for vw_* outcome
  dset[all.vars(formula)[1]] <- 999
  # Unpack the formula
  lf <- lme4::lFormula(formula = formula, data = dset)
  fixed_terms <- colnames(lf$X) # fixed-effects design matrix

  return(fixed_terms)
}

precompile_model <- function(use_model_template,
                             formula,
                             tmp_data,
                             tmp_y,
                             measure,
                             lmm_control,
                             verbose) {

  if (!use_model_template) return(NULL)

  vw_message(" * construct model template", verbose = verbose)
  tmp_data[paste0("vw_", measure)] <- tmp_y  # dummy outcome

  # Fit model once
  model_template <- suppressMessages(
    lme4::lmer(formula = formula,
               data = tmp_data,
               control = lmm_control)
  )

  return(model_template)
}
