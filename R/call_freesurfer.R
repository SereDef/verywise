#' @title Estimate full-width half maximum (FWHM)
#'
#' @description
#' This function is a wrapper for the FreeSurfer command \code{mris_fwhm}. It
#' estimates the full-width half maximum (FWHM) of the brain data based on the
#' residuals of your model.
#'
#' @param outp_dir : output path, where results are stored. The function generates
#' two files: "fwhm.dat" and "finamMask.mgh".
#' @param hemi : hemisphere.
#' @param resid_file : path the .mgh file containing the residuals of the model.
#' @param mask : apply a mask (default = cortical mask + excluding problematic vertices)
#' @param target : (default = "fsaverage") data template.
#' @param verbose : (default = FALSE) verbosity.
#'
#' @return An integer FWHM value
#' @export
#'
estimate_fwhm <- function(outp_dir,
                          hemi,
                          resid_file,
                          mask = NULL,
                          target = "fsaverage",
                          verbose = FALSE) {

  fwhm_estim_path <- file.path(outp_dir, "fwhm.dat")
  final_mask_path <- file.path(outp_dir, "finalMask.mgh")

  cmd_str <- paste("mris_fwhm",
                   "--i", resid_file,
                   "--hemi", hemi,
                   "--subject", target,
                   "--prune",
                   "--cortex",
                   "--dat", fwhm_estim_path,
                   "--out-mask", final_mask_path)

  if (!is.null(mask)) paste(cmd_str, "--mask", mask)

  # message2(verbose = verbose, cmdStr)
  system(cmd_str, ignore.stdout = !verbose)

  fwhm <- round(utils::read.table(fwhm_estim_path))

  return(fwhm)
}

#' @title Compute significant clusters of vertices
#'
#' @description
#' This function is a wrapper for the FreeSurfer command \code{mri_surfcluster}.
#' It computes the clusters of significant vertices.
#'
#' @param outp_dir : output path, where results are stored.
#' @param hemi : hemisphere.
#' @param term_number : which (fixed) term are the clusters computed for.
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
compute_clusters <- function(outp_dir,
                             hemi,
                             term_number,
                             fwhm,
                             FS_HOME,
                             cwp_thr = 0.025,
                             mcz_thr = 30,
                             csd_sign = "abs",
                             mask = NULL,
                             verbose = FALSE) {

  pval_mgh_file <- file.path(
    outp_dir,
    paste(hemi, paste0("stack",term_number), "p.mgh", sep = ".")
  )

  # Montecarlo simulation threshold
  mcz_thr_str <- paste0("th", mcz_thr)

  if (fwhm < 10) fwhm <- paste0("0", fwhm)

  csd <- file.path(FS_HOME, "average", "mult-comp-cor", "fsaverage", hemi, "cortex",
                   paste0("fwhm", fwhm), csd_sign, mcz_thr_str, "mc-z.csd")

  other_files <- paste0("stack", term_number, ".cache.th", mcz_thr, ".", csd_sign, ".sig.",
                        c("cluster.mgh",
                          "voxel.mgh",
                          "cluster.summary",
                          "ocn.mgh",
                          "ocn.annot",
                          "masked.mgh"))

  cmd_str <- paste("mri_surfcluster",
                   "--in", pval_mgh_file,
                   "--csd", csd,
                   "--cwsig", other_files[1],
                   "--vwsig", other_files[2],
                   "--sum", other_files[3],
                   "--ocn", other_files[4],
                   "--oannot", other_files[5],
                   "--annot aparc",
                   "--cwpvalthresh", cwp_thr,
                   "--o", other_files[6],
                   "--no-fixmni",
                   "--surf", "white")

  cmd_str <- if (!is.null(mask)) {
    paste(cmd_str, "--mask", mask)
  } else {
    paste(cmd_str, "--cortex")
  }

  system(cmd_str, ignore.stdout = !verbose)

  NULL
}
