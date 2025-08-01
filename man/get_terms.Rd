% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_terms}
\alias{get_terms}
\title{Unpack \code{lme4} formula}
\usage{
get_terms(formula, data_list)
}
\arguments{
\item{formula}{: model formula object (this should specify a LME model)}

\item{data_list}{: the data, formatted as a list of datasets, as in the
output of \code{\link{imp2list}}.}
}
\value{
A character vector of fixed terms.
}
\description{
Get the terms for the fixed effect for the given formula.
}
\author{
Serena Defina, 2024.
}
