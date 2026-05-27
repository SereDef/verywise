


# ============================= Folder navigators =============================
#' @title
#' List sub-directories till depth \code{n}
#'
#' @description
#' Small recursive utility to list directories given complex folder structure.
#' This is used in \code{\link{build_supersubject}} for listing all sub-directories
#' until the correct level.
#'
#' @param path : the path where to begin search
#' @param n : the depth of the recursive process (i.e., how many levels of sub-directories to list)
#'
#' @author Serena Defina, 2024.
#' @return A list of directory/sub-directory names.
#'
list.dirs.till <- function(path, n) {
  res <- list.dirs(path, recursive = FALSE)

  if (n > 1) {
    add <- list.dirs.till(res, n - 1)
    return(add)
  } else {
    return(res)
  }
}

# ==============================================================================
# Read and save annotation files for internal use
# fs_home = "/Applications/freesurfer/7.4.1"
# hemis <- c("lh", "rh")
# aparc.annot <- lapply(hemis, function(hemi){
#   annot_file <- file.path(fs_home, "subjects/fsaverage/label",
#                           paste0(hemi, ".aparc.annot"))
#   freesurferformats::read.fs.annot(annot_file)})
# names(aparc.annot) <- hemis
# usethis::use_data(aparc.annot, internal = TRUE, overwrite = TRUE)
# ==============================================================================


`%||%` <- function(x, y) if (is.null(x)) y else x


#' @importFrom cli cli_progress_step cli_progress_done
NULL
