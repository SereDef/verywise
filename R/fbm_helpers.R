# ============================== FBM operations ===============================
#' @title
#' Check if all row elements are 0 in FBM
#'
#' @description
#' Check if all row elements in a FBM are 0. This is used to clean up the
#' super-subject matrix in \code{\link{mask_cortex}}.
#'
#' @param X : the file-backed matrix (FBM) object
#' @param n_cores : (default = 1) number of cores for parellalization
#' @param row.ind : indicator for rows
#' @param col.ind : indicator for columns
#' @param row.mask : (default = NULL) specify a subset of rows
#' @param col.mask : (default = NULL) specify a subset of columns
#' @param verbose : (default = TRUE)
#'
#' @importFrom bigstatsr rows_along
#' @importFrom bigstatsr cols_along
#' @importFrom bigstatsr big_apply
#'
#' @author Serena Defina, 2024.
#'
#' @return A (large) logical vector for where rows are all 0.
#'
fbm_col_has_0 <- function(X, n_cores = 1,
                          row.ind = bigstatsr::rows_along(X),
                          col.ind = bigstatsr::cols_along(X),
                          row.mask = NULL, col.mask = NULL,
                          verbose = TRUE) {
  # Any sub-selection?
  if (is.numeric(row.mask)) row.mask <- as.logical(row.mask)
  if (is.numeric(col.mask)) col.mask <- as.logical(col.mask)
  if (!is.null(row.mask)) row.ind <- row.ind[row.mask]
  if (!is.null(col.mask)) col.ind <- col.ind[col.mask]

  problem_verts <- bigstatsr::big_apply(X,
                                        a.FUN = function(X, ind) {
                                          apply(
                                            X[row.ind, col.ind[ind]], 2,
                                            function(q) any(q == 0e-5)
                                          )
                                        },
                                        a.combine = "c",
                                        ncores = n_cores,
                                        ind = seq_along(col.ind)
  )

  if (sum(problem_verts) > 0) {
    vw_message(" * Ignoring ", sum(problem_verts),
               " vertices that contained 0 values.\n   These may be located",
               " at the edge of the cortical map and\n   are potentially",
               " problematic.", verbose = verbose)
  }
  return(problem_verts)
}


build_output_bks <- function(result_path, res_bk_names, verbose = TRUE) {

  vw_message(" * generate file-backed output containers", verbose = verbose)

  res_bk_paths <- paste(result_path, res_bk_names, sep = ".")

  # Note: always remove backing files if they already exist
  res_bk_files <- paste0(res_bk_paths, ".bk")

  if (any(file.exists(res_bk_files))) {
    vw_message("   ! WARNING: overwriting existing results backing files.",
               verbose = verbose)
    file.remove(res_bk_files[file.exists(res_bk_files)])
  }
  names(res_bk_paths) <- res_bk_names

  return(res_bk_paths)
}
