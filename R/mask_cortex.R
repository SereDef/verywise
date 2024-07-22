
#' @title
#' Clean out vertices that are not on the cortex
#'
#' @description
#' Clean the super-subject matrix of the vertices that have 0 value or are not
#' located on the cortex, according to FreeSurfer \code{*cortical.mask.mgh} file
#' distributed in \code{inst/extdata}.
#'
#' @param ss : the super-subject file-backed matrix (FBM) object
#' @param hemi : (default = "lh") hemisphere ("lh" or "rh")
#' @param target : (default = "fsaverage")
#' @param n_cores : (default = 1) number of cores for this operation.
#'
#' @return A masked super-subject file-backed matrix (FBM).
#'
mask_cortex <- function(ss, hemi = "lh", target = "fsaverage", n_cores = 1) {
  # Load mask file
  mask_path <- file.path(system.file("extdata", package = "verywise"), "cortex_mask")
  mask_file <- file.path(mask_path, paste0(hemi, ".", target, ".cortex.mask.mgh"))
  mask <- load.mgh(mask_file)$x

  # mask rows with all 0 and cortex
  mgh_is_0 <- fbm_row_is_0(ss, n_cores = n_cores)
  masked_ss <- !as.logical(mask) | mgh_is_0

  masked_ss
}
