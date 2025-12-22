#' @title
#' Calculate Significant Cluster Statistics
#'
#' @description
#' This function computes statistics (mean or median) for significant clusters based on 
#' the supersubject data and results directories. It retrieves the clusters from an MGH file, 
#' extracts relevant measurements for each cluster, and computes either the mean or median 
#' for each subject in the supersubject data.
#'
#' @param stat A character string specifying the statistic to compute. Either "mean" or "median".
#' @param ss_dir A character string indicating the directory containing the super-subject data matrix.
#' @param res_dir A character string indicating the directory containing the result data, 
#'   including the stack information and the significant clusters.
#' @param term A character string indicating the term or "stack name" to find in the stack file.
#' @param measure A character string indicating the measure used (e.g., "thickness", "volume").
#' @param hemi A character string indicating the hemisphere. One of "lh" for left hemisphere 
#'   or "rh" for right hemisphere.
#'
#' @return A data frame with the computed statistics for each cluster and subject. The data frame 
#'   contains one column for each computed statistic, indexed by the subject's folder ID.
#'
#' @examplesIf exists(ss_dir)
#' result <- significant_cluster_stats(
#'   stat = "mean", 
#'   ss_dir = "path/to/ss_dir", 
#'   res_dir = "path/to/res_dir", 
#'   term = "age", 
#'   measure = "thickness", 
#'   hemi = "lh"
#' )
#'
#' @export
#' 
significant_cluster_stats <- function(stat, ss_dir, res_dir, term, measure, hemi) {
  # Match mesh arguments
  hemi <- match.arg(hemi, c("lh", "rh"))
  stat <- match.arg(stat, c('mean','median'))

  # Validate paths
  if (!dir.exists(res_dir))
    stop("Results directory does not exist: ", res_dir)
  if (!dir.exists(ss_dir))
    stop("Super-subject directory does not exist: ", res_dir)

  # Find term stack index
  stack_file <- file.path(res_dir, "stack_names.txt")
  if (!file.exists(stack_file)) {
    stop(sprintf("Cannot find the `stack_names.txt` in `%s`.", res_dir),
           " Results folder is incorrect or corrupted.")
  }
  # Read existing file
  stack_ids <- utils::read.table(stack_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (!term %in% stack_ids$stack_name) {
    stop("Term '", term, "' not found in stack_names.txt")
  }
  # Extract stack
  stack <- paste0('stack', stack_ids[ stack_ids$stack_name == term, 'stack_number'])

  ocn_map_file <- file.path(res_dir, paste(hemi, measure, stack, "cache.th30.abs.sig.ocn.mgh", sep = "."))
  
  clusters <- load.mgh(ocn_map_file)$x

  if (max(clusters) == 0) {
    message('No significant clusters found')
    return(NULL)
  }

  vw_message(table(clusters))

  ss_file <- file.path(ss_dir, paste(hemi, measure, 'fsaverage.supersubject.rds', sep='.'))
  ss <- bigstatsr::big_attach(ss_file)

  rownames_file <- file.path(ss_dir, paste(hemi, measure, 'ss.rownames.csv', sep='.'))
  ss_rownames <- scan(file = rownames_file, what = character(), sep = "\n", quiet = TRUE)

  output = data.frame('folder_id' = ss_rownames)

  for (cluster in 1:max(clusters)) {
    mask <- clusters == cluster
    ss_subset = ss[, mask]
    if (stat == 'mean') {
      med_meas = rowMeans(ss_subset) 
    } else {
      med_meas <- apply(ss_subset, 1, median)
    }
    output[paste(stat, hemi, measure, paste0('stack', stack), paste0('cluster', cluster))] <- med_meas
  }
  
  return(output)
}
