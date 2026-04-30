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
  
  if (verbose) cli::cli_progress_step("Estimate data smoothness for multiple testing correction")

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

  if (!is.null(mask)) {
    mask_mgh_path <- paste0(result_path, ".inputMask.mgh")
    write_mask_mgh(mask, fs_template = fs_template, file_name = mask_mgh_path)
    # on.exit(if (file.exists(mask_mgh_path)) file.remove(mask_mgh_path), add = TRUE)
    cmd_str <- paste(cmd_str, "--mask", mask_mgh_path)
  }

  # message2(verbose = verbose, cmdStr)
  system(cmd_str, ignore.stdout = !verbose)

  tryCatch({fwhm <- round(utils::read.table(fwhm_estim_path))},
           error = function(e) {
             message("mris_fwhm command not running correctly")
             # TODO: choose a return value in case of error?
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
#' @param fs_template : (default = "fsaverage") data template.
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
#' @param full_surfcluster_output : (default = FALSE) whether to save additional
#'   files: \code{.cluster.mgh}: a map of cluster-wise significances
#'   files: \code{.voxel.mgh}: a map of corrected voxel-wise significance
#'   files: \code{.ocn.annot}: output clusters as an annotation
#'   files: \code{.masked.mgh}: input with non-clusters set to 0
#' @param mask : apply a mask (default = cortical mask + excluding problematic vertices)
#' @param verbose : (default = FALSE) verbosity.
#'
#' @export
#'
compute_clusters <- function(stack_path,
                             hemi,
                             fwhm,
                             FS_HOME,
                             fs_template = "fsaverage",
                             cwp_thr = 0.025,
                             mcz_thr = 30,
                             csd_sign = "abs",
                             full_surfcluster_output = FALSE,
                             mask = NULL,
                             verbose = FALSE) {

  # Compute cluster-wise pvalues
  # (i.e. pvalue of the cluster corrected for multiple comparisons

  # Format monte carlo simulation threshold
  mcz_thr_str <- paste0("th", mcz_thr)

  # Format FWHM (esure a leading 0 if < 10)
  fwhm_str <- paste0("fwhm", sprintf("%02d", as.integer(fwhm)))

  # Cluster Simulation Data. This file is produced by running mri_glmfit with --sim.
  csd_file <- file.path(FS_HOME, "average", "mult-comp-cor", fs_template, hemi,
                        "cortex", fwhm_str, csd_sign, mcz_thr_str, "mc-z.csd")
  
  if (!file.exists(csd_file)) {
    vw_message(c(
      "!" = "Cluster simulation data (CSD) are required for cluster-wise p-value estimation, but no CSD file was found for {.field {fs_template}} template.",
      "i" = "FreeSurfer home: {.file {FS_HOME}}",
      ">" = "Use the default `fsaverage` template, or run {.code mri_glmfit --sim} with
            {.field {fs_template}} to generate the CSD."
    ), verbose = TRUE)
    return(NULL)
  }
  
  pval_mgh_file <- paste0(stack_path, ".-log10p.mgh")

  outp_prefix <- paste0(stack_path, ".cache.th", mcz_thr, ".", csd_sign, ".sig.")

  outp_files <- list("ocn" = paste0(outp_prefix, "ocn.mgh"),
                     "sum" = paste0(outp_prefix, "cluster.summary"))

  # ── Build command ───────────────────────────────────────────────────────────
  args <- c(
    "--in",  pval_mgh_file,
    "--csd", csd_file,
    "--ocn", outp_files[["ocn"]], # cluster number
    "--sum", outp_files[["sum"]], # text summary file
    "--annot", "aparc", # report annotation for max vertex
    "--cwpvalthresh", cwp_thr, # clusterwise threshold
    "--no-fixmni", # <do not> fix MNI talairach coordinates
    "--surf", "white" # get coorindates from surface (white)
  )

  if (full_surfcluster_output) {
    other_files <- list(
      "cwsig"  = paste0(outp_prefix, "cluster.mgh"),
      "vwsig"  = paste0(outp_prefix, "voxel.mgh"),
      "oannot" = paste0(outp_prefix, "ocn.annot"),
      "o"      = paste0(outp_prefix, "masked.mgh")
    )
    args <- c(args,
      "--cwsig",  other_files[["cwsig"]],
      "--vwsig",  other_files[["vwsig"]],
      "--oannot", other_files[["oannot"]],
      "--o",      other_files[["o"]]
    )
  }

  args <- c(args, if (!is.null(mask)) c("--mask", mask) else "--cortex")

  # ── Run and capture output ──────────────────────────────────────────────────
  # stderr is captured separately: suppress the known "seed not set" warning
  # but surface any real errors to the user
  stderr_file <- tempfile("mri_surfcluster_stderr_")
  on.exit(unlink(stderr_file), add = TRUE)

  exit_code <- system2(
    command = "mri_surfcluster",
    args = args,
    stdout = if (verbose) "" else FALSE,  # "" → print, FALSE → discard
    stderr = stderr_file                  # always capture, never blindly discard
  )

  # ── Handle stderr ───────────────────────────────────────────────────────────
  stderr_lines <- readLines(stderr_file, warn = FALSE)

  # Filter the known harmless warning
  noise_pattern <- "supposed to be reproducible but seed not set"
  real_errors   <- stderr_lines[!grepl(noise_pattern, stderr_lines, fixed = FALSE)]
  real_errors   <- real_errors[nzchar(trimws(real_errors))]  # drop blank lines

  # ── Check exit code ─────────────────────────────────────────────────────────
  if (exit_code != 0) {
    abort_msg <- c("mri_surfcluster failed with exit code {exit_code}.",
                   "i" = "Stack: {.file {basename(stack_path)}}")
    
    if (length(real_errors) > 0) {
      abort_msg <- c(abort_msg, "x" = paste(real_errors, collapse = "\n"))
    } else {
      abort_msg <- c(abort_msg, "i" = "No stderr captured. Try re-running {.code compute_clusters}.")
    }

    cli::cli_abort(abort_msg)
  }

  if (length(real_errors) > 0) {
    vw_message(c("!" = "mri_surfcluster stderr:",
      stats::setNames(real_errors, rep("x", length(real_errors)))), verbose = TRUE)
  }
  
  ocn_map <- load.mgh(outp_files[["ocn"]])$x

  return(ocn_map)
}
