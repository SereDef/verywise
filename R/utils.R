# ============================== Load balancing ================================

#' @title
#' Define chunks of vertices for analyses
#'
#' @description
#' This function takes a vector of values representing the dimension of the
#' vertex-wise data to use and splits it into chunks for memory-efficient
#' parallel processing.
#'
#' @param iv : indices of the vertex-wise brain data to use.
#' @param chunk_size : (default = 1000) how big are the chunks
#'
#' @author Serena Defina, 2024.
#'
#' @return A list with n_chunks elements. Each element is a vector of vertex
#'   positions.
#'
make_chunk_sequence <- function(iv, chunk_size = 1000) {

  chunk_seq <- split(iv, ceiling(seq_along(iv) / chunk_size))

  # Add chunk index as attribute (so it can be accessed for progress updates)
  # chunk_seq <- Map(function(chunk, i) {
  #   attr(chunk, "chunk_idx") <- i; chunk
  #   }, chunk_seq, seq_along(chunk_seq))

  for (i in seq_along(chunk_seq)) {
    attr(chunk_seq[[i]], "chunk_idx") <- i
  }

  return(chunk_seq)
}

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

# ======================== simulation helpers ==================================

locate_roi <- function(rois = NULL, n_verts = 163842, verbose = TRUE) {
  # using Desikan-Killiany 36 regions

  all_rois <- aparc.annot$vd_label[1:n_verts]
  roi_counts <- as.data.frame(table(all_rois))
  names(roi_counts) <- c('roi_id', 'vw_count')
  roi_counts['vw_prop'] <- round(roi_counts[,'vw_count']/n_verts, 3)

  notations <- aparc.annot$LUT[, c('LUT_labelname', 'LUT_value')]
  names(notations) <- c('roi_label', 'roi_id')

  roi_lookup <- merge(notations, roi_counts, by = 'roi_id',
                      all.x = FALSE, all.y = TRUE)
  roi_lookup <- roi_lookup[order(-roi_lookup$vw_prop), ]

  # If no input ROIs are provided just return the entire lookup
  if (is.null(rois)) {
    return(roi_lookup)
  }

  # If only one ROI is provided make sure this is a vector
  rois <- as.vector(rois, mode = "character")

  selected_rois <- roi_lookup[roi_lookup[,'roi_label'] %in% rois, ]
  vw_message(' * selected ', sum(selected_rois$vw_count), ' vertices (',
             sum(selected_rois$vw_prop)*100, '%)', verbose = verbose)

  roi_ids <- selected_rois[, 'roi_id']

  roi_locs <- all_rois %in% as.integer(roi_ids)

  return(roi_locs)
}

# ???? =========================================================================
# qdecr_decon <- function(x, y = environment()) {
#   y <- as.list(substitute(x, env = y))
#   y[-1] <- lapply(y[-1], eval)
#   y
# }
#
# # Get function from string =====================================================
# get_function <- function(x, ...) {
#   if(is.character(x)){
#     fn <- strsplit(x, "::")[[1]] # example > "stats" "lm"
#     x <- if(length(fn) == 1) {
#       get(fn[[1]], envir = parent.frame(), mode = "function", ...)
#     } else {
#       get(fn[[2]], envir = asNamespace(fn[[1]]), mode = "function", ...)
#     }
#   }
#   x
# }
#
# do.call2 <- function(what, args, ...) {
#   what <- get_function(what)
#   do.call(what, as.list(args), ...)
# }
#
#
# # ==============================================================================\
# # Set all arguments for given function call (model) and list of user defined
# # arguments (margs)
# set_arguments <- function(margs, model){
#   m <- names(margs)
#   if(is.null(m)) m <- rep("", length(margs))
#
#   f <- methods::formalArgs(get_function(model))
#   f2 <- formals(get_function(model))
#   # Rename arguments if necessary
#   b <- which(m[-1] == "")
#   m[b+1] <- f[!f %in% m[-1]][length(b)]
#   names(margs) <- m
#   # Add default arguments...?
#   margs <- c(margs, f2[!names(f2) %in% m])
#   if(is.symbol(margs$`...`)) margs$`...` <- NULL
#   margs
# }
