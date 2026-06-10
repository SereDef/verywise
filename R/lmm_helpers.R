#' @title Unpack R formula
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
    formula <- drop_random_terms(formula)
  }

  # 2. Drop response (lhs)
  rhs <- stats::delete.response(stats::terms(formula))

  # 3. Build design matrix 
  X <- stats::model.matrix(rhs, data = df)

  if (return_X) return(X) else {
    # Extract term names
    fixed_terms <- colnames(X)
    return(fixed_terms)
  }
}

#' @title Pre-compile a linear mixed model template
#'
#' @description 
#' Fits a `lme4::lmer` model (on each imputed dataset) using a dummy outcome
#' variable. It produces a list of fitted model templates that can be updated
#' with vertex outcome data, avoiding the overhead of re-specifying the 
#' random-effects structure for every vertex/surface measure.
#'
#' @param formula A two-sided `formula` object passed to [lme4::lmer()].
#'   The left-hand side is replaced by a dummy outcome column named
#'   `vw_<measure>` before fitting.
#' @param data_list A `list` of `data.frame` objects, typically constructed
#'   by [imp2list()] Each element must contain all variables referenced by 
#'   `formula`.
#' @param tmp_y A numeric vector used as the dummy outcome. Must have the
#'   same length as the number of rows in each element of `data_list`.
#' @param measure A character string naming the brain surface measure
#'   (e.g. `"thickness"`, `"area"`). Used to construct the temporary
#'   outcome column `vw_<measure>` injected into each imputed dataset.
#' @param weights Either `NULL` (no observation weights), a numeric vector
#'   of weights of the same length as the rows in each imputed dataset, or
#'   a character string giving the name of a column in each dataset that
#'   contains the weights.
#' @param REML Logical. Whether to use restricted maximum likelihood
#'   estimation. Passed directly to [lme4::lmer()]. Defaults to `TRUE`.
#' @param lmm_control A `lmerControl` object produced by
#'   [lme4::lmerControl()], used to fine-tune the optimiser behaviour.
#'   Defaults to `lme4::lmerControl(calc.derivs = FALSE)`.
#' @param verbose Logical. If `TRUE`, emits a `cli` progress step during
#'   model construction.
#'
#' @returns A `list` of the same length as `data_list`, where each element
#'   is a fitted `lmerMod` object estimated on the corresponding imputed
#'   dataset with the dummy outcome. These objects are intended for use
#'   with [lme4::refit()] or equivalent downstream functions in the
#'   `verywise` pipeline.
#' 
#' @author Serena Defina, 2026
#' 
precompile_model <- function(formula,
                             data_list,
                             tmp_y,
                             measure,
                             weights = NULL,
                             REML = TRUE,
                             lmm_control = lme4::lmerControl(calc.derivs=FALSE),
                             verbose) {
  
  if (verbose) cli::cli_progress_step('Construct model template', spinner=TRUE)

  # Fit model template (once per imputed dataset)
  model_template <- lapply(data_list, function(imp) {

    imp[paste0("vw_", measure)] <- tmp_y  # dummy outcome
    
    weight_vec <- if (is.character(weights)) imp[[weights]] else weights
   
    tmp_fit <- suppressMessages(do.call(lme4::lmer, 
      list(formula = formula, data = imp, REML = REML, 
           control = lmm_control, weights = weight_vec)))

    # Embed the control settings into the call so update() can access them ?  
    # tmp_fit@call$control <- lmm_control
    # mp_fit@call$REML <- REML

    return(tmp_fit)
  })

  if (verbose) cli::cli_progress_done()
  model_template
}

# precompile_model <- function(formula,
#                              tmp_data,
#                              tmp_y,
#                              measure,
#                              REML = TRUE,
#                              lmm_control = lme4::lmerControl(calc.derivs=FALSE),
#                              verbose) {
  
#   if (verbose) cli::cli_progress_step('Construct model template', spinner=TRUE)

#   tmp_data[paste0("vw_", measure)] <- tmp_y  # dummy outcome

#   # Fit model once
#   model_template <- suppressMessages(
#     lme4::lmer(formula = formula,
#                data = tmp_data,
#                REML = REML,
#                control = lmm_control)
#   )

#   # Embed the control settings into the call so update() can access them    
#   model_template@call$control <- lmm_control
#   model_template@call$REML <- REML

#   return(model_template)
# }

#' @title Refit a pre-compiled linear mixed model with a new outcome and extract statistics
#'
#' @description
#' Refits a linear mixed model with data from a single vertex outcome and extracts fixed
#' effect and model fit statistics. Errors and warnings are caught and returned in the 
#' output list rather than signalled to the caller, because the function is designed to be 
#' used inside parallel loops (over vertices).
#'
#' @param model_template_i Pre-compiled model object for (MI) dataset `i`, generated
#'   using [precompile_model()]. The random-effects structure is reused as-is; only 
#'   the response vector is replaced.
#' @param y A numeric vector of outcome values representing a single vertex from the
#'   super-subject matrix.
#'
#' @returns A named `list` with one of two shapes:
#' **On error:**
#' \describe{
#'   \item{`error`}{`character(1)`. The error message produced by [lme4::refit()].}
#' }
#' **On success:**
#' \describe{
#'   \item{`stats`}{A `data.frame` with columns `term` (character), `qhat` (numeric, 
#'     fixed-effect estimate), and `se` (numeric, standard error), one row per 
#'     fixed-effect term.}
#'   \item{`resid`}{Named numeric vector of model residuals from [stats::residuals()], 
#'     used downstream for smoothness estimation.}
#'   \item{`model_fit`}{Unnamed numeric vector of length 5 containing, in order: 
#'    `is_singular` (0/1), `AIC`, `ICC`, `R2_marginal`, `R2_conditional`.}
#'   \item{`warning`}{`character` vector of any warning or message strings captured 
#'     during refitting. Zero-length if none occurred. Singular-fit boundary warnings
#'     are suppressed from this vector when singularity is already flagged in `model_fit`.}
#' }
#'
#' @details
#' `verywise` uses an "update" (or rather refit)-based workflow instead of fitting
#' the model from scratch at each vertex. This minimizes repeated parsing and
#' model construction overhead, significantly reducing computation time for
#' large-scale vertex-wise analyses.
#'
#' Variance components used for ICC and R\eqn{^2} are computed directly from
#' [lme4::VarCorr()] with no additional package dependencies:
#' \itemize{
#'   \item \strong{Var(random)} — sum of diagonal elements across all
#'     random-effect covariance matrices.
#'   \item \strong{Var(residual)} — residual variance
#'     \eqn{\hat{\sigma}^2}{sigma^2}.
#'   \item \strong{Var(fixed)} — variance of the marginal linear predictor
#'     \eqn{\mathrm{Var}(\mathbf{X}\hat{\beta})}{Var(X * beta)}.
#' }
#'
#' The three R\eqn{^2}/ICC quantities follow Nakagawa & Schielzeth (2013):
#'
#' \deqn{
#'   \mathrm{ICC} = \frac{\sigma^2_{\mathrm{rand}}}{\sigma^2_{\mathrm{rand}} + \sigma^2_{\varepsilon}}
#' }{ICC = var_rand / (var_rand + var_resid)}
#'
#' \deqn{
#'   R^2_{\mathrm{marginal}} = \frac{\sigma^2_{\mathrm{fix}}}{\sigma^2_{\mathrm{fix}} + \sigma^2_{\mathrm{rand}} + \sigma^2_{\varepsilon}}
#' }{R2_marginal = var_fix / (var_fix + var_rand + var_resid)}
#'
#' \deqn{
#'   R^2_{\mathrm{conditional}} = \frac{\sigma^2_{\mathrm{fix}} + \sigma^2_{\mathrm{rand}}}{\sigma^2_{\mathrm{fix}} + \sigma^2_{\mathrm{rand}} + \sigma^2_{\varepsilon}}
#' }{R2_conditional = (var_fix + var_rand) / (var_fix + var_rand + var_resid)}
#'
#' Division-by-zero or other numeric edge cases in ICC and R\eqn{^2} are handled
#' by the internal helper [safe_calc()].
#' 
#' @references
#' Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for
#' obtaining R² from generalized linear mixed-effects models.
#' *Methods in Ecology and Evolution*, 4(2), 133–142.
#' \doi{10.1111/j.2041-210x.2012.00261.x}
#'
#' @seealso [precompile_model()] and [run_vw_lmm()] for the main interface.
#'
#' @importFrom lme4 refit fixef isSingular VarCorr
#' @importFrom stats vcov residuals AIC var predict
#' 
#' @export
#'
refit_lmm <- function(model_template_i, y) {

  error_msg <- NULL
  warning_msg <- character(0)

  # Refit linear mixed model using `lme4::refit`
  fit <- withCallingHandlers(
    tryCatch(
      refit(model_template_i, newresp = y),
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
  
  # Use lme4 directly – no extra dependencies, no interactive prompts:
  vc <- VarCorr(fit)
  var_rand <- sum(sapply(vc, function(x) sum(diag(as.matrix(x)))))
  var_resid <- attr(vc, "sc")^2
  var_fix <- var(predict(fit, re.form = NA)) # var(as.vector(X %*% beta))

  icc <- safe_calc(var_rand / (var_rand + var_resid))
  r2_margin <- safe_calc(var_fix  / (var_fix + var_rand + var_resid))
  r2_condit <- safe_calc((var_fix + var_rand) / (var_fix + var_rand + var_resid))
  
  # dropping names: singularity, aic, icc, r2_marginal, r2_conditional
  perf <- c(is_singular, aic, icc, r2_margin, r2_condit)
  
  list("stats" = fixed_stats, "resid" = resid, "model_fit" = perf,
       "warning" = warning_msg)

}

#' @title Estimate a single outcome linear mixed model and extract statistics
#'
#' @description
#' Fits a linear mixed model to a single vertex outcome using [lme4::lmer()] and 
#' extracts fixed effects and model fit statistics. This function
#' is called repeatedly across all cortical vertices during vertex-wise
#' analysis in the [run_vw_lmm2()] pipeline.
#'
#' @param imp A data.frame containing the phenotype dataset
#'  (in verywise format).
#' @param y A numeric vector of outcome values representing a single
#'  vertex from the super-subject matrix.
#' @param y_name String with the outcome name (as specified in the formula)
#' @param model_formula A model formula object. This should specify a linear mixed
#'  model `lme4` syntax. Example: `vw_thickness ~ age * sex + site + (1|participant_id)`.
#' @param REML Logical specifying whether to optimize the REML criterion (as opposed to
#'  the log-likelihood). Default: TRUE. Use `REML = FALSE` if you intend to do model
#'  comparison (using AIC output).
#' @param lmm_control Optional list (of correct class, resulting from [lme4::lmerControl()] 
#'  containing control parameters to be passed to [lme4::lmer()] (e.g. optimizer choice, 
#'  convergence criteria, see the `?lmerControl` documentation for details.
#'  Default: (uses default settings).
#' @param weights Optional string or numeric vector of weights for the linear mixed model.
#'   You can use this argument to specify inverse-probability weights. If this
#'   is a string, the function look for a column with that name in the phenotype
#'   data. Note that these are not normalized or standardized in any way.
#'   Default: \code{NULL} (no weights).
#'
#' @returns A named `list` with one of two shapes:
#' **On error:**
#' \describe{
#'   \item{`error`}{`character(1)`. The error message produced by [lme4::refit()].}
#' }
#' **On success:**
#' \describe{
#'   \item{`stats`}{A `data.frame` with columns `term` (character), `qhat` (numeric, 
#'     fixed-effect estimate), and `se` (numeric, standard error), one row per 
#'     fixed-effect term.}
#'   \item{`resid`}{Named numeric vector of model residuals from [stats::residuals()], 
#'     used downstream for smoothness estimation.}
#'   \item{`model_fit`}{Unnamed numeric vector of length 5 containing, in order: 
#'    `is_singular` (0/1), `AIC`, `ICC`, `R2_marginal`, `R2_conditional`.}
#'   \item{`warning`}{`character` vector of any warning or message strings captured 
#'     during refitting. Zero-length if none occurred. Singular-fit boundary warnings
#'     are suppressed from this vector when singularity is already flagged in `model_fit`.}
#' }
#'
#' @seealso [run_vw_lmm2()] for the main interface.
#'
#' @importFrom lme4 lmer fixef isSingular VarCorr
#' @importFrom stats vcov residuals AIC var predict
#' 
#' @export
#'
single_lmm <- function(imp, y, y_name, model_formula = NULL,
    REML, lmm_control, weights = NULL){

  if (!is.null(weights) && !is.numeric(weights)) {
    weights <- imp[[weights]]
  }

  error_msg <- NULL
  warning_msg <- character(0)

  # Fit linear mixed model using `lme4::lmer`
  fit <- withCallingHandlers(
    tryCatch({

      # Add (vertex) outcome to (single) dataset
      imp[y_name] <- y

      lmer(formula = model_formula, data = imp,
           weights = weights, REML = REML, control = lmm_control)},
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
  

  # Use lme4 directly – no extra dependencies, no interactive prompts:
  vc <- VarCorr(fit)
  var_rand <- sum(sapply(vc, function(x) sum(diag(as.matrix(x)))))
  var_resid <- attr(vc, "sc")^2
  var_fix <- var(predict(fit, re.form = NA)) # var(as.vector(X %*% beta))

  icc <- safe_calc(var_rand / (var_rand + var_resid))
  r2_margin <- safe_calc(var_fix  / (var_fix + var_rand + var_resid))
  r2_condit <- safe_calc((var_fix + var_rand) / (var_fix + var_rand + var_resid))
  
  # varpart <- get_variance(fit, tolerance = 1e-12) # more lenient than default

  # icc <- safe_calc(varpart$var.random / (varpart$var.random + varpart$var.residual))

  # r2_margin <- safe_calc(varpart$var.fixed / 
  #   (varpart$var.fixed + varpart$var.random + varpart$var.residual))
  
  # r2_condit <- safe_calc((varpart$var.fixed + varpart$var.random) /
  #       (varpart$var.fixed + varpart$var.random + varpart$var.residual))
  
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


# #' @title Run a single outcome linear mixed model and extract statistics
# #'
# #' @description
# #' Fits a linear mixed model to a single vertex outcome using
# #' \code{lme4::lmer()} and extracts fixed effects statistics. This function
# #' is called repeatedly across all cortical vertices during vertex-wise
# #' analysis.
# #'
# #' @param imp A data.frame containing the phenotype dataset
# #'   (in verywise format).
# #' @param y A numeric vector of outcome values representing a single
# #'   vertex from the super-subject matrix.
# #' @param y_name String with the outcome name (as specified in the formula)
# #' @param model_template Pre-compiled model object for faster
# #'   estimation. \code{single_lmm} uses an "update"-based
# #'   workflow instead of refitting the model from scratch. This minimizes
# #'   repeated parsing and model construction overhead, significantly reducing
# #'   computation time for large-scale vertex-wise analyses.
# #' @param weights Optional string or numeric vector of weights for the linear mixed model.
# #'   You can use this argument to specify inverse-probability weights. If this
# #'   is a string, the function look for a column with that name in the phenotype
# #'   data. Note that these are not normalized or standardized in any way.
# #'   Default: \code{NULL} (no weights).
# #'
# #' @return A list with two elements:
# #'   \itemize{
# #'     \item \code{stats} - A data.frame with columns:
# #'       \itemize{
# #'         \item \code{term} - Fixed effect term names
# #'         \item \code{qhat} - Parameter estimates
# #'         \item \code{se} - Standard errors
# #'       }
# #'     \item \code{resid} - A numeric vector of model residuals
# #'     \item \code{warning} - Warning message(s) if any
# #'   }
# #'
# #' @seealso \code{\link{run_vw_lmm}} for the main interface.
# #'
# #' @importFrom lme4 lmer fixef isSingular VarCorr
# #' @importFrom stats update vcov residuals AIC var predict
# # #' @importFrom insight get_variance
# #' 
# #' @export
# #'
# update_lmm <- function(imp, y, y_name, model_template = NULL, weights = NULL) {

#   if (!is.null(weights) && !is.numeric(weights)) {
#     weights <- imp[[weights]]
#   }

#   # Add (vertex) outcome to (single) dataset
#   imp[y_name] <- y

#   error_msg <- NULL
#   warning_msg <- character(0)

#   # Fit linear mixed model using `lme4::lmer`
#   fit <- withCallingHandlers(
#     tryCatch(
#       update(model_template, data = imp, weights = weights),
#       error = function(e) { 
#         error_msg <<- conditionMessage(e)
#         NULL
#       }),
#     warning = function(w) {
#       warning_msg <<- c(warning_msg, conditionMessage(w))
#       invokeRestart("muffleWarning")
#     },
#     message = function(m) {
#       warning_msg <<- c(warning_msg, conditionMessage(m))
#       invokeRestart("muffleMessage")
#     }
#   )

#   if (!is.null(error_msg)) {
#     result <- list("error" = error_msg)
#     return(result)
#   }

#   coefs <- fixef(fit) # Fixed effects estimates
#   ses   <- sqrt(diag(as.matrix(vcov(fit)))) # Their standard errors
  
#   fixed_stats <- data.frame(
#     term = names(coefs),
#     qhat = as.numeric(coefs),
#     se   = as.numeric(ses),
#     row.names = NULL,
#     check.names = FALSE)
  
#   # Also extract model residuals for smoothness estimation
#   resid <- residuals(fit)

#   # Model performance ------------------------------------
  
#   is_singular <- as.numeric(isSingular(fit))

#   # Add singularity to performance grid, but do not print the warning
#   if (is_singular) {
#     warning_msg <- warning_msg[!grepl("boundary (singular) fit", warning_msg, 
#                                       fixed = TRUE)]
#   }

#   # Add model performance metrics 
#   aic <- AIC(fit)

#   # Use lme4 directly – no extra dependencies, no interactive prompts:
#   vc <- VarCorr(fit)
#   var_rand <- sum(sapply(vc, function(x) sum(diag(as.matrix(x)))))
#   var_resid <- attr(vc, "sc")^2
#   var_fix <- var(predict(fit, re.form = NA)) # var(as.vector(X %*% beta))

#   icc <- safe_calc(var_rand / (var_rand + var_resid))
#   r2_margin <- safe_calc(var_fix  / (var_fix + var_rand + var_resid))
#   r2_condit <- safe_calc((var_fix + var_rand) / (var_fix + var_rand + var_resid))
  
#   # dropping names: singularity, aic, icc, r2_marginal, r2_conditional
#   perf <- c(is_singular, aic, icc, r2_margin, r2_condit)

#   list("stats" = fixed_stats, "resid" = resid, "model_fit" = perf,
#        "warning" = warning_msg)

# }

#' @title Safe calculation helper
#' 
#' @description returns NA if any inputs are NULL/invalid
#' 
#' @keywords internal 
safe_calc <- function(expr) {
  result <- tryCatch(expr, error = function(e) NA_real_)
  if (is.null(result) || length(result) == 0) NA_real_ else result
}

#' @title Remove random terms from formula
#' 
#' @keywords internal 
drop_random_terms <- function(formula) {
  # Remove (... | ...) bar terms 
  rhs_str <- deparse(formula[[length(formula)]], width.cutoff = 500L)
  
  # Remove each (expr | group) block, including surrounding +
  rhs_str <- gsub("\\s*\\+\\s*\\([^)]*\\|[^)]*\\)", "", rhs_str)
  rhs_str <- gsub("\\([^)]*\\|[^)]*\\)\\s*\\+\\s*", "", rhs_str)
  rhs_str <- gsub("^\\([^)]*\\|[^)]*\\)$", "1", rhs_str)
  rhs_str <- trimws(rhs_str)
  if (nchar(rhs_str) == 0L) rhs_str <- "1"

  response <- if (length(formula) == 3L) deparse(formula[[2L]]) else NULL
  f <- stats::reformulate(rhs_str, response = response)
  environment(f) <- environment(formula)
  f
}
