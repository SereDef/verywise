
#' @title
#' Clean out vertices that are not on the cortex
#'
#' @description
#' Clean the super-subject matrix of the vertices that have 0 value or are not
#' located on the cortex, according to FreeSurfer \code{*cortical.mask.mgh} file
#' distributed in \code{inst/extdata}. This function is used in
#' \code{\link{build_supersubject}}.
#'
#' @param hemi : (default = "lh") hemisphere ("lh" or "rh")
#' @param target : (default = "fsaverage")
#'
#' @return A (large) logical vector for cortical (or technical, i.e. vertex = 0) mask.
#'
mask_cortex <- function(hemi = "lh", target = "fsaverage") {
  # Load mask file
  mask_file <- file.path(system.file("extdata", package = "verywise"), "cortex_mask",
                         paste0(hemi, ".", target, ".cortex.mask.mgh"))
  # Turn into logical vector
  mask <- as.logical(load.mgh(mask_file)$x) # TRUE where cortex

}
