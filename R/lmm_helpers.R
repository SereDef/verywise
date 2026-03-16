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
#' @param y_name String with the outcome name (as specified in the formula)
#' @param model_template Pre-compiled model object for faster
#'   estimation. \code{single_lmm} uses an "update"-based
#'   workflow instead of refitting the model from scratch. This minimizes
#'   repeated parsing and model construction overhead, significantly reducing
#'   computation time for large-scale vertex-wise analyses.
#' @param weights Optional string or numeric vector of weights for the linear mixed model.
#'   You can use this argument to specify inverse-probability weights. If this
#'   is a string, the function look for a column with that name in the phenotype
#'   data. Note that these are not normalized or standardized in any way.
#'   Default: \code{NULL} (no weights).
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
#' @importFrom lme4 lmer fixef isSingular
#' @importFrom stats update vcov residuals AIC
#' @importFrom insight get_variance
#' 
#' @export
#'
single_lmm <- function(
    imp, y, y_name,
    model_template = NULL,
    weights = NULL)
{

  if (!is.null(weights) && !is.numeric(weights)) {
    weights <- imp[, weights]
  }

  # Add (vertex) outcome to (single) dataset
  imp[y_name] <- y


  error_msg <- NULL
  warning_msg <- character(0)

  # Fit linear mixed model using `lme4::lmer`
  fit <- withCallingHandlers(
    tryCatch(
      update(model_template, data = imp, weights = weights),
      error = function(e) { 
        error_msg <<- conditionMessage(e)
        NULL
      }),
    warning = function(w) {
      warning_msg <<- c(warning_msg, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      warning_msg <<- c(warning_msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  if (!is.null(error_msg)) {
    result <- list("error" = error_msg)
    return(result)
  }

  coefs <- fixef(fit) # Fixed effects estimates
  ses   <- sqrt(diag(as.matrix(vcov(fit)))) # Their standard errors
  
  fixed_stats <- data.frame(
    term = names(coefs),
    qhat = as.numeric(coefs),
    se   = as.numeric(ses),
    row.names = NULL,
    check.names = FALSE)
  
  # Also extract model residuals for smoothness estimation
  resid <- residuals(fit)

  # Model performance ------------------------------------
  
  is_singular <- as.numeric(isSingular(fit))

  # Add singularity to performance grid, but do not print the warning
  if (is_singular) {
    warning_msg <- warning_msg[!grepl("boundary (singular) fit", warning_msg, 
                                      fixed = TRUE)]
  }

  # Add model performance metrics 
  aic <- AIC(fit)
  # icc <- icc(fit)[['ICC_adjusted']]
  # r2 <- t(r2_nakagawa(fit))
  
  varpart <- get_variance(fit, tolerance = 1e-12) # more lenient than default

  icc <- safe_calc(varpart$var.random / (varpart$var.random + varpart$var.residual))

  r2_margin <- safe_calc(varpart$var.fixed / 
    (varpart$var.fixed + varpart$var.random + varpart$var.residual))
  
  r2_condit <- safe_calc((varpart$var.fixed + varpart$var.random) /
        (varpart$var.fixed + varpart$var.random + varpart$var.residual))
  
  # dropping names: singularity, aic, icc, r2_marginal, r2_conditional
  perf <- c(is_singular, aic, icc, r2_margin, r2_condit)
  
  # Extract residual degrees of freedom for Barnard-Rubin adjustment
  # Normally, this would be the number of independent observation minus the
  # number of fitted parameters, but not exactly what is done here.
  # Following the `broom.mixed` package approach, which `mice::pool` relies on
  # resid_df <- df.residual(fit)

  # TODO: implement "stack of interest"?

  # coef(fit)$id # A matrix of effects by random variable
  # lme4::ranef(fit) # extract random effects (these should sum to 0)
  # lme4::VarCorr(fit) # estimated variances, SDs, and correlations between the random-effects terms

  list("stats" = fixed_stats, "resid" = resid, "model_fit" = perf,
       "warning" = warning_msg)

}


#' @title
#' Unpack R formula
#'
#' @description
#' Get the terms for the (fixed) effect for the given formula.
#' Wors for `lme4` as well as regular `lm` style formulas.
#'
#' @param formula A formula object, for example specifying a LME model.
#' @param df A data.frame, for example the first element in the list
#'   output of \code{\link{imp2list}}.
#' @param return_X Logical. Whether to return the whole design matrix or
#'   just the (fixed) term names.
#'
#' @return A design matrix or a character vector of (fixed) term names.
#'
#' @author Serena Defina, 2026.
#'
unpack_formula <- function(formula, df, return_X = FALSE) {

  # 1. Drop random terms (if present)
  # width.cutoff to ensure long formulas are returned as a single 
  # character string rather than split after 60 char.
  if (any(grepl("\\|", deparse(formula, width.cutoff = 500L)))) {
    formula <- lme4::nobars(formula)
  }

  # 2. Drop response (lhs)
  rhs <- stats::delete.response(stats::terms(formula))

  # 3. Build design matrix 
  X <- stats::model.matrix(rhs, data = df)

  if (return_X) {
    return(X)
  } else {
    # Extract term names
    fixed_terms <- colnames(X)
    return(fixed_terms)
  }
}


precompile_model <- function(formula,
                             tmp_data,
                             tmp_y,
                             measure,
                             REML = TRUE,
                             lmm_control = lme4::lmerControl(calc.derivs=FALSE),
                             verbose) {

  vw_message(" * construct model template", verbose = verbose)
  tmp_data[paste0("vw_", measure)] <- tmp_y  # dummy outcome

  # Fit model once
  model_template <- suppressMessages(
    lme4::lmer(formula = formula,
               data = tmp_data,
               REML = REML,
               control = lmm_control)
  )

  # Embed the control settings into the call so update() can access them    
  model_template@call$control <- lmm_control
  model_template@call$REML <- REML

  return(model_template)
}

# Safe calculation helper - returns NA if any inputs are NULL/invalid
safe_calc <- function(expr) {
  result <- tryCatch(expr, error = function(e) NA_real_)
  if (is.null(result) || length(result) == 0) NA_real_ else result
}
