#' @title
#' Run vertex-wise linear mixed model using \code{lme4::lmer()}
#'
#' @description
#' This is is the main function for conducting vertex-wise linear mixed model
#' analyses on brain surface metrics. It will first check use inputs, prepare
#' the phenotype data(list) and run a linear mixed model at each vertex of the
#' specified hemisphere using the \code{\link{single_lmm}} function.
#'
#' The function supports analysis of both single and multiple imputed datasets.
#' It also automatically handles cortical masking, and provides cluster-wise
#' correction for multiple testing using FreeSurfer's Monte Carlo simulation
#' approach.
#'
#' @param formula A model formula object. This should specify a linear mixed
#'   model \code{lme4} syntax. The outcome variable should be one of the
#'   supported brain surface metrics (see Details). Example:
#'   \code{vw_thickness ~ age * sex + site + (1|participant_id)}.
#' @param pheno Either a \code{data.frame}/\code{tibble} containing the
#'   "phenotype" data (i.e., already loaded in the global environment), or a
#'   string specifying the file path to phenotype data. Supported file formats:
#'   .rds, .csv, .txt, .sav (SPSS).
#'   The data should be in \strong{long} format and it should contain all the
#'   variables specified in the \code{formula} plus the \code{folder_id} column.
#' @param subj_dir Character string specifying the path to FreeSurfer data
#'   directory. Must follow the verywise directory structure (see package
#'   vignette for details).
#' @param outp_dir Character string specifying the output directory for results.
#'   If \code{NULL} (default), creates a "verywise_results" subdirectory within
#'   \code{subj_dir}.
#' @param hemi Character string specifying which hemisphere(s) to analyze.
#'   Options: \code{"both"} (default), \code{"lh"} (left only), \code{"rh"}
#'   (right only).
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
#' @param use_model_template Logical indicating whether to pre-compile the model
#'   template for faster estimation.
#'   Default: \code{TRUE} (recommended).
#' @param weights Optional string or numeric vector of weights for the linear mixed model.
#'   You can use this argument to specify inverse-probability weights. If this
#'   is a string, the function look for a column with that name in the phenotype
#'   data. Note that these are not normalized or standardized in any way.
#'   Default: \code{NULL} (no weights).
#' @param lmm_control Optional list (of correct class, resulting from
#'   \code{lmerControl()} containing control parameters to be passed to
#'   \code{lme4::lmer()} (e.g. optimizer choice, convergence criteria,
#'   see the \code{*lmerControl} documentation for details.
#'   Default: (uses default settings).
#' @param seed Integer specifying the random seed for reproducibility
#'   Default: 3108.
#' @param n_cores Integer specifying the number of CPU cores for parallel
#'   processing.
#'   Default: 1.
#' @param chunk_size Integer specifying the number of vertices processed per
#'   chunk in parallel operations. Larger values use more memory but may be
#'   faster.
#'   Default: 1000.
#' @param FS_HOME Character string specifying the FreeSurfer home directory.
#'   Defaults to \code{FREESURFER_HOME} environment variable.
#' @param fwhm Numeric value specifying the full-width half-maximum for
#'   smoothing kernel. Default: 10.
#' @param mcz_thr Numeric value for Monte Carlo simulation threshold.
#'   Any of the following are accepted (equivalent values separate by \code{/}):
#'   \itemize{
#'   \item 13 / 1.3 / 0.05,
#'   \item 20 / 2.0 / 0.01,
#'   \item 23 / 2.3 / 0.005,
#'   \item 30 / 3.0 / 0.001, \* default
#'   \item 33 / 3.3 / 0.0005,
#'   \item 40 / 4.0 / 0.0001.
#' }
#' @param cwp_thr Numeric value for cluster-wise p-value threshold (on top of
#'   all corrections). Set this to 0.025 when both hemispheres are analyzed,
#'   0.05 for single hemisphere.
#'   Default: 0.025.
#' @param save_optional_cluster_info Logical indicating whether to save additional
#'  output form \code{mri_surfcluster} call. See \code{\link{compute_clusters}}
#'  for details. Default: \code{FALSE}.
#' @param save_ss Logical indicating whether to save the super-subject matrix as
#'   an .rds file for faster future processing.
#'   Default: \code{TRUE} (recommended for repeated analyses).
#' @param save_residuals Logical indicating whether to save the residuals.mgh
#'   file. Default: \code{FALSE}.
#' @param verbose Logical indicating whether to display progress messages.
#'   Default: \code{TRUE}.
#'
#' @details
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
#' \strong{Statistical Approach:}
#' The function uses \code{lme4::lmer()} for mixed-effects modeling, enabling
#' analysis of longitudinal and hierarchical data. P-values are computed using
#' the t-as-z approximation, with cluster-wise correction applied using
#' FreeSurfer's Monte Carlo simulation approach.
#'
#' \strong{Multiple Imputation:}
#' The function automatically detects and handles multiple imputed datasets
#' (created with \code{mice} or similar packages), pooling results according
#' to Rubin's rules.
#'
#' \strong{Parallel processing:}
#' The \code{verywise} package employs a carefully designed parallelization
#' strategy to maximize computational efficiency while avoiding the
#' performance penalties associated with nested parallelization.
#' Left and right cortical hemispheres are processed sequentially by default.
#' Parallel processing of the two hemispheres (and/or different metrics, models)
#' should be handled by the user (e.g., using SLURM job arrays or similar,
#' see vignette on parallelisation [COMING UP]).
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
#' \strong{Output Files:}
#' Results are saved in FreeSurfer-compatible .mgh format for visualization
#' with [verywiseWIZard](https://github.com/SereDef/verywise-wizard),
#' FreeView or other neuroimaging software.
#'
#' @return A list of file-backed matrices (\code{bigstatsr::FBM} objects)
#' containing pooled coefficients, SEs, t- and p- values and residuals.
#' Results are also automatically saved to disk in .mgh format.
#'
#' @note
#' \itemize{
#'   \item Ensure FreeSurfer is properly installed and the
#'     \code{FREESURFER_HOME} environment variable is set.
#'   \item Large datasets may require substantial memory. Consider adjusting
#'     \code{chunk_size} and \code{n_cores} based on your system specifications.
#'   \item For reproducibility, always specify a \code{seed}.
#' }
#'
#' @seealso
#' \code{\link{single_lmm}} for single-vertex modeling,
#' \code{vignette("03-run-vw-lmm", package = "verywise")} for detailed
#' usage examples.
#'
#' @examples
#' # Basic cortical thickness analysis
#' \dontrun{
#' results <- run_vw_lmm(
#'   formula = vw_thickness ~ age + sex + site + (1|participant_id),
#'   pheno = "path/to/phentype/data", # or data.frame object
#'   subj_dir = "/path/to/freesurfer/subjects",
#'   outp_dir = "/path/to/output",
#'   hemi = "lh",
#'   n_cores = 4)
#' }
#'
#' @author Serena Defina, 2024.
#'
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#'
#' @export
#'
run_vw_lmm <- function(
  # Basic settings
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
  # Modeling settings
  use_model_template = TRUE,
  weights = NULL,
  lmm_control = lme4::lmerControl(),
  # Reproducibility and parallel processing
  seed = 3108,
  n_cores = 1,
  chunk_size = 1000,
  # Cluster estimation
  FS_HOME = Sys.getenv("FREESURFER_HOME"),
  fwhm = 10,
  mcz_thr = 30,
  cwp_thr = 0.025,
  # Output control
  save_optional_cluster_info = FALSE,
  save_ss = TRUE,
  save_residuals = FALSE,
  verbose = TRUE
) {

  hemi <- match.arg(hemi)

  hemi_name <- if (hemi == "lh") "left" else "right"
  vw_message(pretty_message(paste(hemi_name, "hemisphere")), verbose = verbose)

  # Check user input ===========================================================
  vw_message("Checking user inputs...", verbose = verbose)

  measure <- check_formula(formula)

  ss_file <- paste(hemi, measure, fs_template, "supersubject.rds", sep = ".")

  ss_in_subj_dir <- check_path(subj_dir,
                               file_exists = c(ss_file,
                                               file.path('ss', ss_file)))

  if (is.null(outp_dir)) {
    outp_dir <- file.path(subj_dir, "verywise_results")
    vw_message(" ! WARNING: outpur directory unspecified, which is not recommended.",
               " You can find the results at ", outp_dir)
    dir.create(outp_dir, showWarnings = FALSE)
  }

  ss_in_outp_dir <- check_path(outp_dir, create_if_not = TRUE,
                               file_exists = file.path('ss', ss_file))

  check_freesurfer_setup(FS_HOME, verbose = verbose)

  n_cores <- check_cores(n_cores)

  check_numeric_param(seed, integer = TRUE, lower = 0)
  check_numeric_param(chunk_size, integer = TRUE,
                      lower = 1, upper = 5000) # for memory safety
  check_numeric_param(fwhm, lower = 1, upper = 30)
  # check_numeric_param(mcz_thr, lower = 0)
  # check_numeric_param(cwp_thr, lower = 0)

  # Read phenotype data (if not already loaded) ================================
  vw_message("Checking and preparing phenotype dataset...", verbose = verbose)

  if (is.character(pheno)) pheno <- load_pheno_file(pheno)

  # Transform to list of dataframes (imputed and single datasets alike)
  data_list <- imp2list(pheno); rm(pheno)

  # Check that the data_list is not empty, it contains data.frames of the same
  # dims, and that "folder_id" and all variables specified in the formula are
  # present in the data
  check_data_list(data_list, folder_id, formula)

  # Extract first dataset
  data1 <- data_list[[1]]

  vw_message(" * ", length(data_list), " datasets of dimention: ", nrow(data1),
             " x ", ncol(data1), verbose = verbose)

  check_weights(weights, data1)

  fixed_terms <- get_terms(formula, data1)

  # Check that the stacks are not overwritten by mistake and
  # Save the stack names (i.e. fixed terms) to a lookup file
  check_stack_file(fixed_terms, outp_dir)

  # Build supersubject =========================================================
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  # Read and clean vertex data =================================================

  vw_message("Checking and preparing brain surface data...", verbose = verbose)

  ss_file_locs <- c(ss_in_subj_dir, ss_in_outp_dir)
  ss_exists <- ss_file_locs[!is.null(ss_file_locs)]

  if (!is.null(ss_exists)) {
    vw_message(" * reading super-subject file from: ", ss_exists,
               verbose = verbose)

    ss <- bigstatsr::big_attach(ss_exists)

    ss_folder <- dirname(ss_exists)

  } else {

    ss_folder <- file.path(outp_dir, "ss")

    ss <- build_supersubject(
      subj_dir = subj_dir,
      folder_ids = data1[, folder_id],
      supsubj_dir = ss_folder,
      measure = measure,
      hemi = hemi,
      n_cores = n_cores,
      fwhmc = paste0("fwhm", fwhm),
      fs_template = fs_template,
      save_rds = save_ss,
      error_cutoff = tolerate_surf_not_found
    )
  }

  ss_rownames <- scan(file = file.path(ss_folder,
                                       paste(hemi, measure, 'ss.rownames.csv', sep = '.')),
                      what = character(), sep = "\n", quiet = TRUE)

  if (!identical(ss_rownames, data1[, folder_id])) {
    # Assume all data.frames in data_list have the same order... one hopes
    vw_message(' * matching phenotype with available brain surface data')
    data_list <- lapply(data_list,
                        function(df) {  # Match the row names in ss
                          return(df[match(ss_rownames, df[, folder_id]), ])})
  }

  vw_message(" * cleaning super-subject matrix...", verbose = verbose)

  # Cortical mask
  is_cortex <- mask_cortex(hemi = hemi, fs_template = fs_template)

  # Additionally check that there are no vertices that contain any 0s
  problem_verts <- fbm_col_has_0(ss)

  good_verts <- which(!problem_verts & is_cortex); rm(problem_verts)

  # Unpack model ===============================================================
  vw_message("Statistical model preparation...",
             "\n * call: ", deparse(formula), verbose = verbose)

  # Number of vertices
  vw_n <- length(is_cortex); rm(is_cortex)
  # Number of terms (excluding random terms)
  fe_n <- length(fixed_terms)
  # Number of participants*timepoint (long format)
  n_obs <- nrow(data1)
  # Number of (imputed) datasets
  m <- length(data_list)

  # "Pre-compile the model"
  # cache the model frame to avoid re-generating it them each time
  # single_lmm can leverage an "update"-based workflow to minimize
  # repeated parsing and model construction overhead
  model_template <- precompile_model(use_model_template,
    formula = formula, tmp_data = data1, tmp_y = ss[,1],
    measure = measure, lmm_control = lmm_control, verbose = verbose)

  # Prepare FBM output =========================================================

  result_path <- file.path(outp_dir, paste(hemi, measure, sep = "."))

  # Temporary output matrices
  res_bk_names <- c("coef", "se", "t", "p", "resid")
  res_bk_paths <- build_output_bks(result_path, res_bk_names = res_bk_names,
                                   verbose = verbose)
  on.exit(file.remove(paste0(res_bk_paths, ".bk")), add = TRUE)

  c_vw <- bigstatsr::FBM(fe_n, vw_n, init = 0,
                         backingfile = res_bk_paths["coef"])  # Coefficients
  s_vw <- bigstatsr::FBM(fe_n, vw_n, init = 0,
                         backingfile = res_bk_paths["se"])    # Standard errors
  t_vw <- bigstatsr::FBM(fe_n, vw_n, init = 0,
                         backingfile = res_bk_paths["t"])     # t values
  p_vw <- bigstatsr::FBM(fe_n, vw_n, init = 1,
                         backingfile = res_bk_paths["p"])     # P values
  r_vw <- bigstatsr::FBM(n_obs, vw_n, init = 0,
                         backingfile = res_bk_paths["resid"]) # Residuals

  log_file <- paste0(result_path, ".issues.log") # Log model fitting issues

  # Prepare chunk sequence =====================================================
  vw_message(" * chunk dataset", verbose = verbose)
  chunk_seq <- make_chunk_sequence(good_verts, chunk_size = chunk_size)

  # Parallel analyses ==========================================================
  vw_message("Running analyses...", verbose = verbose)
  vw_message("* dimentions: ", n_obs, " observations x ",
             length(good_verts), " (of ", vw_n, " total) vertices.",
             verbose = verbose)

  # log_file <- file.path(outp_dir, paste0(hemi,".", measure,"model.log"))

  if (n_cores > 1) vw_message("* preparing cluster of ", n_cores, " workers...",
                              verbose = verbose)
  # Set up parallel processing
  cluster <- parallel::makeCluster(n_cores, type = "FORK",
                                   outfile = paste0(result_path, ".progress.log"))
  # NOTE: would not work on Windows anyway till freesurfer dependency is needed
  # type = ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK"))

  doParallel::registerDoParallel(cluster)
  # Perform %dopar% as %dorng% loops, for reproducible random numbers
  doRNG::registerDoRNG(seed)

  # Sync library paths: ensures all workers look for packages in the same place
  # as the main session
  # parallel::clusterCall(cluster, function(x) .libPaths(x), .libPaths())

  vw_message("* fitting linear mixed models...")

  # utils::capture.output(pb <- utils::txtProgressBar(0, length(chunk_seq),
  #                                                   style = 3), file = "/dev/null")

  # Progress bar setup # note progressr only works with doFuture not doParallel
  chunk <- NULL
  foreach::foreach(chunk = chunk_seq,
                   .packages = c("bigstatsr"),
                   .export = c("single_lmm", "vw_pool", "vw_message")
                   ) %dopar% {

    # Progress updates
    if (verbose) {
      chunk_idx <- as.integer(attr(chunk, 'chunk_idx'))
      worker_id <- Sys.getpid() # useful for debugging

      vw_message(sprintf("   - Processing chunk %d/%d (worker: %s)",
                         chunk_idx, length(chunk_seq), worker_id))
      utils::flush.console()
    }

    for (v in chunk) {

      # NOTE: ss does not need to be copied to each worker with doParallel
      vertex <- ss[, v]

      # Loop through imputed datasets and run analyses
      out_stats <- lapply(data_list, single_lmm,
                          y = vertex,
                          formula = formula,
                          model_template = model_template,
                          weights = weights,
                          lmm_control = lmm_control,
                          pvalues = (m == 1))

      # Pool results
      pooled_stats <- vw_pool(out_stats, m = m)

      # Log errors (if any)
      if (is.character(pooled_stats)) {
        cat(paste0(v, "\t", pooled_stats, "\n"), file = log_file,
            append = TRUE)
        # & skip to the next value of v
        next
      }

      # Log warnings (if any)
      if (!is.null(pooled_stats$warning)) {
        cat(paste0(v, "\t", pooled_stats$warning, "\n"), file = log_file,
            append = TRUE)
      }

      # Write results to their respective FBM
      c_vw[, v] <- pooled_stats$coef
      s_vw[, v] <- pooled_stats$se
      t_vw[, v] <- pooled_stats$t
      p_vw[, v] <- -1 * log10(pooled_stats$p)
      r_vw[, v] <- pooled_stats$resid

    }
  }

  parallel::stopCluster(cluster)

  out <- list(c_vw, s_vw, t_vw, p_vw, r_vw)
  # "coefficients", "standard_errors", "t_values", "p_values", "residuals"
  names(out) <- res_bk_names

  # Post-processing ============================================================

  # Split all coefficients into separate .mgh files
  vw_message("Post-processing\n * converting coefficients, SEs, t- and p-values",
             " to .mgh format", verbose = verbose)

  results_grid <- expand.grid(seq_along(fixed_terms), # how many terms
                              # all statistics except the residuals
                              res_bk_names[-length(res_bk_names)])

  stats_mgh_paths <- paste(result_path, paste0("stack", results_grid[, 1]),
                           results_grid[, 2], "mgh", sep = ".")

  lapply(seq_len(nrow(results_grid)), function(grid_row) {
    fbm2mgh(fbm = out[[results_grid[grid_row, 2]]],
            fname = stats_mgh_paths[grid_row],
            filter = as.integer(results_grid[grid_row, 1]))
  })

  vw_message(" * converting residuals to .mgh format...", verbose = verbose)
  resid_mgh_path <- paste(result_path, "residuals.mgh", sep = ".")
  fbm2mgh(fbm = out[[res_bk_names[length(res_bk_names)]]],
          fname = resid_mgh_path)
  if (!save_residuals) on.exit(file.remove(resid_mgh_path), add = TRUE)

  # Estimate full-width half maximum (using FreeSurfer) ========================
  vw_message("Estimating data smoothness for multiple testing correction...",
             verbose = verbose)

  fwhm <- estimate_fwhm(result_path = result_path,
                        hemi = hemi,
                        mask = good_verts,
                        fs_template = fs_template)

  # Clamp fwhm to [1, 30]
  fwhm_clamped <- min(max(fwhm, 1), 30)

  if (fwhm != fwhm_clamped) {
    direction <- if (fwhm > 30) "high. Reduced to 30." else "low. Increased to 1."
    vw_message(sprintf(" * estimated smoothness is %s, which is really %s",
                       fwhm, direction), verbose = verbose)
    fwhm <- fwhm_clamped
  } else {
    vw_message(sprintf(" * estimated smoothness = %s", fwhm), verbose = verbose)
  }

  # Apply cluster-wise correction (using FreeSurfer) ===========================
  vw_message("Clusterwise correction...", verbose = verbose)

  for (stack_n in seq_along(fixed_terms)){
    stack_path <- paste0(result_path, ".stack", stack_n)
    compute_clusters(stack_path = stack_path,
                     hemi = hemi,
                     fwhm = fwhm,
                     FS_HOME = FS_HOME,
                     mcz_thr = mcz_thr,
                     cwp_thr = cwp_thr,
                     full_surfcluster_output = save_optional_cluster_info,
                     mask = paste0(result_path, ".finalMask.mgh"))
  }

  vw_message(pretty_message("All done! :)"))

  return(out)
}


#' @title
#' Run a single linear mixed model and extract statistics
#'
#' @description
#' Fits a linear mixed model to a single vertex outcome using
#' \code{lme4::lmer()} and extracts fixed effects statistics. This function
#' is called repeatedly across all cortical vertices during vertex-wise
#' analysis.
#'
#' @param imp A data.frame containing the phenotype dataset
#'   (in verywise format).
#' @param y A numeric vector of outcome values representing a single
#'   vertex from the super-subject matrix.
#' @param formula An R formula object describing the linear mixed model using
#'   \code{lme4} notation.
#' @param model_template Optional pre-compiled model object for faster
#'   estimation. When provided, \code{single_lmm} will use an "update"-based
#'   workflow instead of refitting the model from scratch. This minimizes
#'   repeated parsing and model construction overhead, significantly reducing
#'   computation time for large-scale vertex-wise analyses.
#'   Default: \code{NULL}.
#' @param pvalues Logical indicating whether to compute p-values using the
#'   *t-as-z* approximation (see Details). Default: \code{TRUE}.
#' @inheritParams run_vw_lmm
#'
#' @details
#' No additional parameters are currently passed to the \code{lme4::lmer} call.
#' P-values are estimated using the *t-as-z* approach. This is known
#' to the anti-conservative for small sample sizes.sizes but provides a
#' computationally efficient solution. Type I error control is addressed
#' more rigorously at the cluster-forming stage.
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item \code{stats} - A data.frame with columns:
#'       \itemize{
#'         \item \code{term} - Fixed effect term names
#'         \item \code{qhat} - Parameter estimates
#'         \item \code{se} - Standard errors
#'         \item \code{tval} - t-statistics (if \code{pvalues = TRUE})
#'         \item \code{pval} - p-values (if \code{pvalues = TRUE})
#'       }
#'     \item \code{resid} - A numeric vector of model residuals
#'   }
#'
#' @seealso \code{\link{run_vw_lmm}} for the main interface.
#'
#' @importFrom lme4 lmer
#' @importFrom lme4 fixef
#' @importFrom stats vcov
#' @importFrom stats pnorm
#' @importFrom stats residuals
# #' @importFrom stats df.residual
# #' @importFrom car Anova
#'
#' @export
#'
single_lmm <- function(
  imp, y, formula,
  model_template = NULL,
  pvalues = TRUE,
  weights = NULL,
  lmm_control = lme4::lmerControl())
  {

  if (!is.null(weights) && !is.numeric(weights)) {
    weights <- imp[, weights]
  }

  # Add (vertex) outcome to (single) dataset
  imp[all.vars(formula)[1]] <- y

  warning_msg <- NULL
  error_msg <- NULL

  # Fit linear mixed model using `lme4::lmer`
  fit <- withCallingHandlers(
    tryCatch(
      {
        if (!is.null(model_template)) {
          stats::update(model_template, data = imp, weights = weights)
        } else {
          suppressMessages(lmer(formula = formula, data = imp,
                                weights = weights, control = lmm_control))
        }
      },
      error = function(e) {
        error_msg <<- conditionMessage(e)
        NULL
      }
    ),
    warning = function(w) {
      warning_msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )

  if (!is.null(error_msg)) {
    result <- list("stats" = NA, "resid" = NA, "error" = error_msg)
    return(result)
  }

  # coef(fit)$id # A matrix of effects by random variable
  # lme4::ranef(fit) # extract random effects (these should sum to 0)
  # lme4::VarCorr(fit) # estimated variances, SDs, and correlations between the random-effects terms

  # Extract estimates, standard errors and p-values
  stats <- data.frame(
    "qhat" = as.matrix(fixef(fit)), # Fixed effects estimates
    "se" = as.matrix(summary(fit)$coefficients[, "Std. Error"]) # standard errors
    # microbench
    # sqrt(diag(summary(fit)$vcov))
    # sqrt(diag(as.matrix(vcov(fit)))))
  )

  # "pval" = as.matrix(Anova(fit, type = 3)[, "Pr(>Chisq)"]) # Wald chi-square tests
  # TODO: check pbkrtest:: Kenward-Roger and Satterthwaite Based estimation
  # TODO: ask Bing: why type 3 and not 2
  # -log10() ?
  # TODO: implement stack of interest?

  # Calculate p-values using t-as-z method
  if (pvalues) {
    stats$tval <- stats$qhat / stats$se
    stats$pval <- 2 * (1 - pnorm(abs(stats$tval)))
  }

  # Save row.names (i.e. terms) as a column so these can be grouped later
  stats <- data.frame("term" = row.names(stats), stats, row.names = NULL)

  # Also extract model residuals for smoothness estimation
  resid <- residuals(fit)

  # Finally, extract residual degrees of freedom for Barnard-Rubin adjustment
  # Normally, this would be the number of independent observation minus the
  # number of fitted parameters, but not exactly what is done here.
  # Following the `broom.mixed` package approach, which `mice::pool` relies on
  # df <- df.residual(fit)

  result <- list("stats" = stats, "resid" = resid, "warning" = warning_msg)

  return(result)
}
