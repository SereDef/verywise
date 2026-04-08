#' @title Run vertex-wise linear model on local site data
#'
#' @description
#' This the "local" function for conducting a distributed linear model
#' analyses on brain surface metrics. It will first check use inputs, prepare
#' the phenotype data and compute vertex-wise sufficient statistics for the
#' specified hemisphere. These can be then compressed, shared and finally
#' aggregated using using the \code{\link{run_vw_fed_aggr}} function.
#'
#' @param site_name A character string indicating the site name or ID.
#' @param formula A model formula object. This should specify a linear
#'   model. The outcome variable should be one of the
#'   supported brain surface metrics (see Details). Example:
#'   \code{vw_thickness ~ age * sex + ethnicity}.
#' @param pheno Either a \code{data.frame}/\code{tibble} containing the
#'   "phenotype" data (i.e., already loaded in the global environment), or a
#'   string specifying the file path to phenotype data. Supported file formats:
#'   .rds, .csv, .txt, .sav (SPSS).
#'   The data contain all the variables specified in the left-hand side of the 
#'   \code{formula} (i.e., after the `~`) plus the \code{folder_id} column.
#' @param subj_dir Character string specifying the path to FreeSurfer data
#'   directory. Must follow the verywise directory structure (see package
#'   vignette for details).
#' @param outp_dir Character string specifying the output directory for results.
#'   If \code{NULL} (default), creates a "verywise_results" sub-directory in the
#'   current working directory (not recommended).
#' @param hemi Character string specifying which hemisphere to analyze.
#'   Options: `"lh"` (left hemisphere: default), `"rh"` (right hemisphere).
#' @param fs_template Character string specifying the FreeSurfer template for
#'   vertex registration. Options:
#'   \itemize{
#'   \item \code{"fsaverage"} (default) = 163842 vertices (highest resolution),
#'   \item \code{"fsaverage6"} = 40962 vertices,
#'   \item \code{"fsaverage5"} = 10242 vertices,
#'   \item \code{"fsaverage4"} = 2562 vertices,
#'   \item \code{"fsaverage3"} = 642 vertices
#'   }
#'   Note that lower resolutions should be only used to downsample the brain
#'   map, for faster model tuning. The final analyses should also run using
#'   \code{fs_template = "fsaverage"} to avoid (small) imprecisions in vertex
#'   registration and smoothing.
#' @param apply_cortical_mask Logical indicating whether to exclude non-cortical
#'   vertices from analysis. Default: \code{TRUE} (recommended).
#' @param folder_id Character string specifying the column name in \code{pheno}
#'   that contains subject directory names of the input neuroimaging data
#'   (e.g. "sub-001_ses-baseline" or "site1/sub-010_ses-F1"). These are expected
#'   to be nested inside \code{subj_dir}.
#'   Default: \code{"folder_id"}.
#' @param tolerate_surf_not_found Integer indicating how many brain surface
#'   files listed in \code{folder_id} can be missing from \code{subj_dir}. If
#'   the number of missing or corrupted files is
#'   \code{ > tolerate_surf_not_found } execution will stop.
#'   Default: \code{20L}.
#' @param seed Integer specifying the random seed for reproducibility
#'   Default: 3108.
#' @param n_cores Integer specifying the number of CPU cores for parallel
#'   processing.
#'   Default: 1.
#' @param chunk_size Integer specifying the number of vertices processed per
#'   chunk in parallel operations. Larger values use more memory but may be
#'   faster.
#'   Default: 1000.
#' @param fwhm Numeric value specifying the full-width half-maximum for
#'   smoothing kernel. This is used to read FreeSurfer files. Default: 10.
#' @param save_ss Logical indicating whether to save the super-subject matrix
#'  ("ss") as an .rds file that can be then re-used in future analyses. This can
#'  also be a character string specifying the directory where ss should be saved.
#'  When \code{TRUE}, the ss matrix will be saved in \code{<outp_dir>/ss} by
#'  default. Default: \code{FALSE}.
#' @param verbose Logical indicating whether to display progress messages.
#'   Default: \code{TRUE}.
#'
#' @details
#' The function does not currently support multiple imputed datasets or IPW weights
#' (this is for future development)
#' 
#' \strong{Supported Brain Surface Metrics:}
#' The outcome specified in \code{formula} should be a brain surface metric
#' among:
#' \itemize{
#'   \item \code{vw_thickness} - Cortical thickness
#'   \item \code{vw_area} - Cortical surface area (white  surface)
#'   \item \code{vw_area.pial} - Cortical surface area (pial surface)
#'   \item \code{vw_curv} - Mean curvature
#'   \item \code{vw_jacobian_white} - Jacobian determinant (white surface)
#'   \item \code{vw_pial} - Pial surface coordinates
#'   \item \code{vw_pial_lgi} - Local gyrification index (pial surface)
#'   \item \code{vw_sulc} - Sulcal depth
#'   \item \code{vw_volume} - Gray matter volume
#'   \item \code{vw_w_g.pct} - White/gray matter intensity ratio
#'   \item \code{vw_white.H} - Mean curvature (white surface)
#'   \item \code{vw_white.K} - Gaussian curvature (white surface)
#' }
#'
#' \strong{Parallel processing:}
#' The \code{verywise} package employs a carefully designed parallelization
#' strategy to maximize computational efficiency while avoiding the
#' performance penalties associated with nested parallelization.
#' Left and right cortical hemispheres are processed sequentially by default.
#' Parallel processing of the two hemispheres (and/or different metrics, models)
#' should be handled by the user (e.g., using SLURM job arrays or similar,
#' see vignette on parallelization).
#' Within each hemisphere, vertices are divided into chunks of size
#' \code{chunk_size} and processed in parallel across \code{n_cores} workers
#' (when \code{n_cores > 1}). When multiple imputed datasets are present,
#' these are processed sequentially within each vertex.
#'
#' Note that, on some systems, implicit parallelism in low-level matrix algebra
#' libraries (BLAS/LAPACK) can interfere with explicit parallelization. If you
#' feel like processing is taking too long, I recommend disabling these implicit
#' threading libraries before starting R.
#' For example:
#' \preformatted{
#' export OPENBLAS_NUM_THREADS=1
#' export OMP_NUM_THREADS=1
#' export MKL_NUM_THREADS=1
#' export VECLIB_MAXIMUM_THREADS=1
#' export NUMEXPR_NUM_THREADS=1
#' }
#'
#' Also note that using a very large number of cores (e.g. >120) may sometimes
#' cause worker initialization or other issues (e.g. R parallel processes
#' limits)
#'
#' @return A list of site-specific information summary matrices, of which
#' some are (\code{bigstatsr::FBM} objects). These should be compressed
#' before sending them to the aggregation center.
#'
#' @note
#' \itemize{
#'   \item Large datasets may require substantial memory. Consider adjusting
#'     \code{chunk_size} and \code{n_cores} based on your system specifications.
#'   \item For reproducibility, always specify a \code{seed}.
#' }
#'
#' @seealso
#' \code{\link{chunk_Ymats}} for vertex-chunk modeling,
# #' \code{vignette("03-run-vw-lmm", package = "verywise")} for detailed
# #' usage examples.
#'
#'
#' @author Serena Defina, 2026.
#'
#' @importFrom foreach %dopar%
#' 
#' @export
#'
run_vw_fed_local <- function( 
  # Basic settings
  site_name,
  formula,
  pheno,
  subj_dir,
  outp_dir = NULL,
  # Brain data processing
  hemi = c("lh", "rh"),
  fs_template = "fsaverage",
  apply_cortical_mask = TRUE,
  folder_id = "folder_id",
  tolerate_surf_not_found = 20,
  fwhm = 10,
  # Reproducibility and parallel processing
  seed = 3108,
  n_cores = 1,
  chunk_size = 1000,
  save_ss = FALSE,
  verbose = TRUE) {
  
  # Check user input ===========================================================

  hemi <- match.arg(hemi)
  measure <- check_formula(formula)

  # TODO enforce site_name; define term names to ensure harmonisarion

  vw_pretty_message(outcome_name(hemi, measure), verbose = verbose)
  vw_message('* Model: ', deparse(formula), verbose = verbose)
  # vw_pretty_message('', fill = '-', verbose = verbose)

  vw_message("Checking user inputs...", verbose = verbose)

  ss_file <- paste(hemi, measure, fs_template, "supersubject.rds", sep = ".")

  subj_dir <- check_path(subj_dir)
  ss_exists <- check_ss_exists(subj_dir, ss_file)

  outp_dir <- check_path(outp_dir, create_if_not = TRUE)
  
  n_cores <- check_cores(n_cores)
  
  # Other set-up stuff =========================================================

  # Avoid bigstatsr wanrining about lost precision (float vs. double)
  old_opts <- options(bigstatsr.downcast.warning = FALSE)
  on.exit(options(old_opts), add = TRUE)

  # Esure reproducible seeds in parallel settings
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  # Read phenotype data (if not already loaded) ================================
  vw_message("Checking and preparing phenotype dataset...", verbose = verbose)

  if (is.character(pheno)) pheno <- load_pheno_file(pheno)

  # Check that "folder_id" and all variables specified in the formula are
  # present in the data
  check_data_frame(pheno, folder_id, formula)

  vw_message(" * Dataset of dimention: ", nrow(pheno), " x ", ncol(pheno), verbose = verbose)
  
  folder_ids <- pheno[, folder_id, drop=TRUE] # ensure this is always a character vector 

  # Read and clean vertex data =================================================

  vw_message("Checking and preparing brain surface data...", verbose = verbose)

  if (is.character(save_ss)) {
    ss_dir <- check_path(save_ss, create_if_not = TRUE)
    save_ss <- TRUE
  } else {
    ss_dir <- file.path(outp_dir, 'ss')
    if (!save_ss) on.exit(unlink(ss_dir, recursive = TRUE), add = TRUE)
  }

  if (ss_exists) {
    vw_message(" * reading super-subject file from: ", subj_dir,
               verbose = verbose)

    ss <- subset_supersubject(
      supsubj_dir = subj_dir,
      supsubj_file = ss_file,
      folder_ids = folder_ids,
      new_supsubj_dir = ss_dir,
      n_cores = n_cores,
      save_rds = save_ss,
      error_cutoff = tolerate_surf_not_found,
      verbose = verbose)

  } else {

    vw_message("Building ", outcome_name(hemi, measure), " super-subject matrix...",
               verbose = verbose)

    ss <- build_supersubject(
      subj_dir = subj_dir,
      folder_ids = folder_ids,
      supsubj_dir = ss_dir,
      measure = measure,
      hemi = hemi,
      fs_template = fs_template,
      n_cores = n_cores,
      fwhmc = paste0("fwhm", fwhm),
      save_rds = save_ss,
      error_cutoff = tolerate_surf_not_found,
      verbose = verbose
    )
  }

  vw_message(" * cleaning super-subject matrix...", verbose = verbose)

  # Cortical mask
  is_cortex <- mask_cortex(hemi = hemi, fs_template = fs_template)

  # Additionally check that there are no vertices that contain any 0s
  problem_verts <- fbm_col_has_0(ss, n_cores = 1L)

  good_verts <- which(!problem_verts & is_cortex); rm(problem_verts)

  # Ensure phenotype and ss row order matches ==================================
  pheno <- check_row_match(ss_file = ss$bk, pheno = pheno, folder_ids = folder_ids)

  # Unpack model ===============================================================
  vw_message("Statistical model preparation...",
             "\n * call: ", deparse(formula), verbose = verbose)
  
  # Get design matrix (p x p)
  X <- unpack_formula(formula, pheno, return_X = TRUE)

  # Number of terms
  p <- ncol(X)
  # Number of observations
  n_obs = nrow(pheno)

  # Collect "static" i.e. not vertex dependent site info
  site_info <- list(
    # Number of participants
    n_obs = n_obs,
    # Model
    terms = colnames(X), # dim: p
    XtX = crossprod(X),  # dim: p x p
    # Note, only in models without intercept this is != XtX[,1]
    X1 = drop(crossprod(X, rep(1, n_obs))), # dim: p 
    fs_template = fs_template,
    # Effective number of vertices that were computed
    n_good_vx = length(good_verts),
    date_created = format(Sys.time(), "%d-%m-%Y %H:%M"),
    verywise_version = as.character(utils::packageVersion("verywise"))
  )

  # Prepare FBM output =========================================================
  
  # Number of vertices
  vw_n <- length(is_cortex); rm(is_cortex)

  result_path <- file.path(outp_dir, paste(site_name, hemi, measure, sep = "."))

  # Temporary output matrices
  res_bk_names <- c("XtY", "psumsY") # 1tY and YtY
  res_bk_paths <- build_output_bks(result_path, res_bk_names = res_bk_names,
                                   verbose = verbose)

  fbm_precision <- "float" # single precision – 32 bits

  XtY <- bigstatsr::FBM(p, vw_n, init = 0, type = fbm_precision,
                        backingfile = res_bk_paths["XtY"]) 
  psumsY <- bigstatsr::FBM(2, vw_n, init = 0, type = fbm_precision,
                        backingfile = res_bk_paths["psumsY"])
  
  # log_file <- paste0(result_path, ".issues.log") # Log model fitting issues

  # Prepare chunk sequence =====================================================
  vw_message(" * chunk dataset", verbose = verbose)
  chunk_seq <- make_chunk_sequence(good_verts, chunk_size = chunk_size)

  # Parallel analyses ==========================================================
  vw_message("Running analyses...\n",
             " * dimentions: ", site_info$n_obs, " observations x ", site_info$n_good_vx,
             " (of ", vw_n, " total) vertices.", verbose = verbose)

  # progress_file <- paste0(result_path, ".progress.log")
  # on.exit(if (file.exists(progress_file)) file.remove(progress_file), add = TRUE)

  with_parallel(n_cores = n_cores, 
    seed = seed,
    verbose = verbose, 
    expr = {
      foreach::foreach(chunk = chunk_seq, 
        .packages = c("bigstatsr"), 
        .export = c("chunk_Ymats", "X")
                   # "init_progress_tracker", "update_progress_tracker")
    ) %dopar% { # Only parallel if n_cores > 1

      # Progress updates
      # progress_tracker <- init_progress_tracker(chunk, chunk_seq, verbose=TRUE)

      ss_chunk <- ss[, chunk]

      chunk_out <- chunk_Ymats(X, ss_chunk)
      
      # Write results to their respective FBM
      XtY[, chunk] <- chunk_out$XtY
      psumsY[, chunk] <- chunk_out$psumsY
    } 
  })
  
  # Write output to outp_dir
  saveRDS(site_info, paste0(result_path,'.static.rds'))
  XtY$save()
  psumsY$save()

  return(c(site_info, XtY = XtY, psumsY = psumsY))
}