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
