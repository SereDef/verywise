#' @title
#' Convert imputation object to a list of dataframes
#'
#' @description
#' Converts known imputation objects (e.g., \code{amelia}, \code{aregImpute},
#' \code{mi}, \code{mids}, \code{missForest}) to a \code{list} of dataframes.
#'
#' @details
#' This function attempts to unify imputation object formats in supplying
#' datasets. This is done by extracting the imputed datasets from the imputation
#' object and assembling them into a \code{list}. If an unknown object is
#' supplied, the function will throw an error with the request to supply a list
#' of data frames instead.
#'
#' This is a generic function: methods can be defined for it directly,
#' see \code{methods(imp2list)} for a list of all available methods.
#'
#' @param obj The imputation object that contains the imputed datasets
#'
#' @return A list of imputed datasets
#'
#' @author Sander Lamballais, 2018.
#' @export

imp2list <- function(obj) UseMethod("imp2list", obj)

#' @export
imp2list.amelia <- function(obj) obj$imputations

#' @export
imp2list.aregImpute <- function(obj) {
  stop(
    "The aregImpute format is not supported, ",
    "as it does not include the raw dataset. ",
    "Please create your own design matrices."
  )
}

#' @export
imp2list.data.frame <- function(obj) list(obj)

#' @export
imp2list.default <- function(obj) {
  x <- class(obj)[1]
  stop(
    "The provided data seems of class `", x, "`. Please create a list object of ",
    "your imputed datasets and supply that instead."
  )
}

#' @export
imp2list.list <- function(obj) {
  if (length(obj) > 1) {
    s <- t(sapply(obj, names))
    if (all(is.null(s))) stop("The provided datasets (in the list) do not have column names")
    if (any(is.null(s))) stop("Some provided datasets (in the list) do not have column names")
    count <- sum(!duplicated(s))
    if (count == 1) s else stop("The supplied datasets (in the list) do not share column names.")
  }
  obj
}

#' @export
imp2list.matrix <- function(obj) imp2list(as.data.frame(obj))

#' @export
imp2list.mi <- function(obj) {
  if (!requireNamespace("mi", quietly = TRUE)) {
    stop("You are trying to load in an `mi` output object, but `mi` is not installed.")
  }
  mi::complete(obj)
}

#' @export
imp2list.mids <- function(obj) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("You are trying to load in an `mice` output object (class `mids`), but `mice` is not installed.")
  }
  lapply(seq_len(obj$m), function(y) mice::complete(obj, y))
}

#' @export
imp2list.missForest <- function(obj) list(obj$objimp)
