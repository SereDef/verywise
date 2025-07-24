#' @title Estimate full-width half maximum (FWHM)
#'
#' @description
#' This function is a wrapper for the FreeSurfer command \code{mris_fwhm}. It
#' estimates the full-width half maximum (FWHM) of the brain data based on the
#' residuals of your model.
#'
#' @param result_path : path where model residuals are stored and where results
#'    are generated: this includes two files:
#'    "hemi.measure.fwhm.dat" and "hemi.measure.finamMask.mgh".
#' @param hemi : hemisphere.
#' @param mask : apply a mask (default = cortical mask + excluding problematic vertices)
#' @param fs_template : (default = "fsaverage") data template.
#' @param verbose : (default = FALSE) verbosity.
#'
#' @return An integer FWHM value
#' @export
#'
estimate_fwhm <- function(result_path,
                          hemi,
                          mask = NULL,
                          fs_template = "fsaverage",
                          verbose = FALSE) {

  fwhm_estim_path <- paste0(result_path, ".fwhm.dat")
  final_mask_path <- paste0(result_path, ".finalMask.mgh")

  resid_mgh_path <- paste0(result_path, ".residuals.mgh")

  cmd_str <- paste("mris_fwhm",
                   "--i", resid_mgh_path,
                   "--hemi", hemi,
                   "--subject", fs_template,
                   "--prune",
                   "--cortex",
                   "--dat", fwhm_estim_path,
                   "--out-mask", final_mask_path)

  if (!is.null(mask)) paste(cmd_str, "--mask", mask)

  # message2(verbose = verbose, cmdStr)
  system(cmd_str, ignore.stdout = !verbose)

  tryCatch({fwhm <- round(utils::read.table(fwhm_estim_path))},
           error = function(e) {
             message("mris_fwhm command not running correctly")
             # Choose a return value in case of error
             NA} )

  return(fwhm)
}

#' @title Compute significant clusters of vertices
#'
#' @description
#' This function is a wrapper for the FreeSurfer command \code{mri_surfcluster}.
#' It computes the clusters of significant vertices.
#'
#' @param stack_path : path where this term's p values are stored and where
#'    results are generated.
#' @param hemi : hemisphere.
#' @param fwhm : full-width half maximum estimate of data smoothness.
#' @param FS_HOME : FreeSurfer directory, i.e. \code{$FREESURFER_HOME}
#' @param mcz_thr : (default = 0.001) numeric value for the Monte Carlo simulation threshold.
#' Any of the following are accepted (equivalent values separate by `/`):
#'  * 13 / 1.3 / 0.05,
#'  * 20 / 2.0 / 0.01,
#'. * 23 / 2.3 / 0.005,
#'  * 30 / 3.0 / 0.001, \* default
#'  * 33 / 3.3 / 0.0005,
#'  * 40 / 4.0 / 0.0001.
#' @param cwp_thr : (default = 0.025, when both hemispheres are ran, else 0.05)
#' the cluster-wise p-value threshold on top of all corrections.
#' @param csd_sign : (default = "abs")
#' @param mask : apply a mask (default = cortical mask + excluding problematic vertices)
#' @param verbose : (default = FALSE) verbosity.
#'
#' @export
#'
compute_clusters <- function(stack_path,
                             hemi,
                             fwhm,
                             FS_HOME,
                             cwp_thr = 0.025,
                             mcz_thr = 30,
                             csd_sign = "abs",
                             mask = NULL,
                             verbose = FALSE) {

  pval_mgh_file <- paste0(stack_path, ".p.mgh")

  # Format montecarlo simulation threshold
  mcz_thr_str <- paste0("th", mcz_thr)

  # Format FWHM (esure a leading 0 if < 10)
  fwhm_str <- paste0("fwhm", sprintf("%02d", fwhm))

  csd_file <- file.path(FS_HOME, "average", "mult-comp-cor", "fsaverage", hemi,
                        "cortex", fwhm_str, csd_sign, mcz_thr_str, "mc-z.csd")

  other_files <- paste0(stack_path, ".cache.th", mcz_thr, ".", csd_sign, ".sig.",
                        c("cluster.mgh", # map of cluster-wise significances
                          "voxel.mgh", # map of corrected voxel-wise significances
                          "cluster.summary", # text summary file
                          "ocn.mgh", # value is cluster number
                          "ocn.annot", # output clusters as an annotation
                          "masked.mgh") # input with non-clusters set to 0
                        )

  cmd_str <- paste("mri_surfcluster",
                   "--in", pval_mgh_file,
                   "--csd", csd_file,
                   "--cwsig", other_files[1],
                   "--vwsig", other_files[2],
                   "--sum", other_files[3],
                   "--ocn", other_files[4],
                   "--oannot", other_files[5],
                   "--annot aparc", # report annotation for max vertex
                   "--cwpvalthresh", cwp_thr, # clusterwise threshold
                   "--o", other_files[6],
                   "--no-fixmni", # <do not> fix MNI talairach coordinates
                   "--surf", "white") # get coorindates from surface (white)

  cmd_str <- if (!is.null(mask)) {
    paste(cmd_str, "--mask", mask)
  } else {
    paste(cmd_str, "--cortex")
  }

  system(cmd_str, ignore.stdout = !verbose)

  NULL
}


