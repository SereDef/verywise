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
  progress_marks <- ceiling(chunk_size * c(0.25, 0.50, 0.75))

  for (i in seq_along(chunk_seq)) {

    attr(chunk_seq[[i]], "chunk_idx") <- i

    if (length(chunk_seq[[i]]) != chunk_size) { # e.g. last chunk
      this_chunk_size <- length(chunk_seq[[i]])
      progress_marks <- ceiling(this_chunk_size * c(0.25, 0.50, 0.75))
    }

    attr(chunk_seq[[i]], "chunk_25%") <- chunk_seq[[i]][progress_marks[1]]
    attr(chunk_seq[[i]], "chunk_50%") <- chunk_seq[[i]][progress_marks[2]]
    attr(chunk_seq[[i]], "chunk_75%") <- chunk_seq[[i]][progress_marks[3]]

  }

  return(chunk_seq)
}


# ========================== External environment ==============================

with_tmp_sysenv <- function(tmp_sysenv = c(OMP_NUM_THREADS = 1,
                                           MKL_NUM_THREADS = 1,
                                           OPENBLAS_NUM_THREADS = 1,
                                           VECLIB_MAXIMUM_THREADS = 1,
                                           NUMEXPR_NUM_THREADS = 1),
                            code) {
  sysenv <- Sys.getenv(names(tmp_sysenv), unset = NA, names = TRUE)

  for (nm in names(tmp_sysenv)) Sys.setenv(nm = tmp_sysenv[[nm]])

  on.exit({
    for (nm in names(sysenv)) {
      if (is.na(sysenv[[nm]])) Sys.unsetenv(nm)
      else Sys.setenv(nm = sysenv[[nm]])
    }
  }, add = TRUE)

  force(code)
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

# ======================== Simulation helpers ==================================

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

# ==============================================================================

count_vertices <- function(fs_template) {

  n_verts <- switch(fs_template,
                    fsmicro = 100,
                    fsaverage = 163842,
                    fsaverage6 = 40962,
                    fsaverage5 = 10242,
                    fsaverage4 = 2562,
                    fsaverage3 = 642)
  return(n_verts)
}
