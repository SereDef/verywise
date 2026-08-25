#' @title Run vertex-wise linear mixed model using [lme4::lmer()]
#'
#' @description
#' This is the main function for conducting vertex-wise linear mixed model
#' analyses on brain surface metrics. It will first check use inputs, prepare
#' the phenotype data(list) and the brain (outcome) data, then it will fit a 
#' linear mixed model at each vertex of the specified hemisphere using the 
#' [refit_lmm()] function, extract key statistics and perform multiple testing
#' correction.
#'
#' The function supports analysis of both single and multiple imputed datasets.
#' It also automatically handles cortical masking, and provides two multiple 
#' testing strategies: FDR or cluster-wise correction using FreeSurfer's Monte 
#' Carlo simulation approach.
#'
#' @param formula A model formula object. This should specify a linear mixed
#'   model `lme4` syntax. The outcome variable should be one of the
#'   supported brain surface metrics (see Details). Example:
#'   `vw_thickness ~ age * sex + site + (1|participant_id)`.
#' @param pheno Either a `data.frame`/`tibble` containing the
#'   "phenotype" data (i.e., already loaded in the global environment), or a
#'   string specifying the file path to phenotype data. Supported file formats:
#'   .rds, .csv, .txt, .sav (SPSS).
#'   The data should be in **long** format and it should contain all the
#'   variables specified in the left-hand side of the `formula` (i.e., after the `~`) 
#'   plus the `folder_id` column.
#' @param subj_dir Character string specifying the path to FreeSurfer data
#'   directory. Must follow the verywise directory structure (see package
#'   vignette for details).
#' @param outp_dir Character string specifying the output directory for results.
#'   If `NULL` (default), creates a "verywise_results" sub-directory in the
#'   current working directory (not recommended).
#' @param hemi Character string specifying which hemisphere to analyze.
#'   Options: `"lh"` (left hemisphere: default), `"rh"` (right hemisphere).
#' @param fs_template Character string specifying the FreeSurfer template for
#'   vertex registration. Options:
#'   \itemize{
#'   \item `"fsaverage"` (default) = 163842 vertices (highest resolution),
#'   \item `"fsaverage6"` = 40962 vertices,
#'   \item `"fsaverage5"` = 10242 vertices,
#'   \item `"fsaverage4"` = 2562 vertices,
#'   \item `"fsaverage3"` = 642 vertices
#'   }
#'   Note that lower resolutions should be only used to downsample the brain
#'   map, for faster model tuning. The final analyses should also run using
#'   `fs_template = "fsaverage"` to avoid (small) imprecisions in vertex
#'   registration and smoothing.
#' @param apply_cortical_mask Logical indicating whether to exclude non-cortical
#'   vertices from analysis. Default: `TRUE` (recommended).
#' @param folder_id Character string specifying the column name in `pheno`
#'   that contains subject directory names of the input neuroimaging data
#'   (e.g. "sub-001_ses-baseline" or "site1/sub-010_ses-F1"). These are expected
#'   to be nested inside `subj_dir`.
#'   Default: `"folder_id"`.
#' @param tolerate_surf_not_found Integer indicating how many brain surface
#'   files listed in `folder_id` can be missing from `subj_dir`. If
#'   the number of missing or corrupted files is ` > tolerate_surf_not_found ` 
#'   execution will stop.
#'   Default: `20L`.
#' @param weights Optional string or numeric vector of weights for the linear mixed model.
#'   You can use this argument to specify inverse-probability weights. If this
#'   is a string, the function look for a column with that name in the phenotype
#'   data. Note that these are not normalized or standardized in any way.
#'   Default: `NULL` (no weights).
#' @param lmm_control Optional list (of correct class, resulting from
#'   [lme4::lmerControl()] containing control parameters to be passed to
#'   [lme4::lmer()] (e.g. optimizer choice, convergence criteria,
#'   see the `?lmerControl` documentation for details.
#'   Default: (uses default settings).
#' @param REML Logical specifying whether to optimize the REML criterion (as opposed to
#'   the log-likelihood). Default: TRUE. Use `REML = FALSE` if you intend to do model
#'   comparison (using AIC output).
#' @param seed Integer specifying the random seed for reproducibility
#'   Default: 3108.
#' @param n_cores Integer specifying the number of CPU cores for parallel
#'   processing.
#'   Default: 1.
#' @param chunk_size Integer specifying the number of vertices processed per
#'   chunk in parallel operations. Larger values use more memory but may be
#'   faster.
#'   Default: 1000.
#' @param mtc Character string: multiple testing correction strategy. 
#'   Options: `"mcz"` (FreeSurfer MonteCarlo-based cluster correction: default), `"fdr"` (False Discovery Rate).
#' @param FS_HOME Character string specifying the FreeSurfer home directory.
#'   Defaults to `FREESURFER_HOME` environment variable.
#' @param fwhm Numeric value specifying the full-width half-maximum for
#'   smoothing kernel. Default: 10.
#' @param mcz_thr Numeric value for the Monte Carlo simulation threshold. 
#' Any of the following are accepted (equivalent values are separated by `/`):
#'   \itemize{
#'   \item 13 / 1.3 / 0.05,
#'   \item 20 / 2.0 / 0.01,
#'   \item 23 / 2.3 / 0.005,
#'   \item 30 / 3.0 / 0.001, \* default
#'   \item 33 / 3.3 / 0.0005,
#'   \item 40 / 4.0 / 0.0001.
#' }
#' @param cwp_thr Numeric value for cluster-wise p-value threshold (on top of all
#'  corrections). Set this should be set to `0.025` when both hemispheres are analyzed,
#'  and `0.05` for single hemisphere analyses.
#'  Default: `0.025`.
#' @param save_optional_cluster_info Logical indicating whether to save additional
#'  output form `mri_surfcluster` call. See [compute_clusters()]
#'  for details. Default: `FALSE`.
#' @param save_ss Logical indicating whether to save the super-subject matrix
#'  ("ss") as an .rds file that can be then re-used in future analyses. This can
#'  also be a character string specifying the directory where ss should be saved.
#'  When `TRUE`, the ss matrix will be saved in `<outp_dir>/ss` by
#'  default. Default: `FALSE`.
#' @param save_residuals Logical indicating whether to save the residuals.mgh
#'   file. Default: `FALSE`.
#' @param save_cov Character vector of two fixed effect term names of `NULL`. 
#'   When specified, the covariance between the two terms is extracted and saved
#'   for later analysis (e.g. simple slopes). 
#' @param verbose Logical indicating whether to display progress messages.
#'   Default: `TRUE`.
#'
#' @details
#' \strong{Supported Brain Surface Metrics:}
#' The outcome specified in `formula` should be a brain surface metric
#' among:
#' \itemize{
#'   \item `vw_thickness` - Cortical thickness
#'   \item `vw_area` - Cortical surface area (white  surface)
#'   \item `vw_area.pial` - Cortical surface area (pial surface)
#'   \item `vw_curv` - Mean curvature
#'   \item `vw_jacobian_white` - Jacobian determinant (white surface)
#'   \item `vw_pial` - Pial surface coordinates
#'   \item `vw_pial_lgi` - Local gyrification index (pial surface)
#'   \item `vw_sulc` - Sulcal depth
#'   \item `vw_volume` - Gray matter volume
#'   \item `vw_w_g.pct` - White/gray matter intensity ratio
#'   \item `vw_white.H` - Mean curvature (white surface)
#'   \item `vw_white.K` - Gaussian curvature (white surface)
#' }
#'
#' \strong{Statistical Approach:}
#' The function uses [lme4::lmer()] for mixed-effects modeling, enabling
#' analysis of longitudinal and hierarchical data. P-values are computed using
#' the t-as-z approximation, with cluster-wise correction applied using
#' FreeSurfer's Monte Carlo simulation approach.
#'
#' \strong{Multiple Imputation:}
#' The function automatically detects and handles multiple imputed datasets
#' (created with `mice` or similar packages), pooling results according
#' to Rubin's rules.
#'
#' \strong{Parallel processing:}
#' The `verywise` package employs a carefully designed parallelization
#' strategy to maximize computational efficiency while avoiding the
#' performance penalties associated with nested parallelization.
#' Left and right cortical hemispheres are processed sequentially by default.
#' Parallel processing of the two hemispheres (and/or different metrics, models)
#' should be handled by the user (e.g., using SLURM job arrays or similar,
#' see vignette on parallelization).
#' Within each hemisphere, vertices are divided into chunks of size
#' `chunk_size` and processed in parallel across `n_cores` workers
#' (when `n_cores > 1`). When multiple imputed datasets are present,
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
#' using `verywise` plotting functions ([plot_vw_map()],[plot_vw_diff()]...), 
#' the [verywiseWIZard](https://github.com/SereDef/verywise-wizard),
#' FreeView or other neuroimaging software.
#'
#' @return A list of file-backed matrices ([bigstatsr::FBM] objects)
#' containing pooled coefficients, SEs, t- and p- values and residuals.
#' Results are also automatically saved to disk in .mgh format.
#'
#' @note
#' \itemize{
#'   \item Ensure FreeSurfer is properly installed and the
#'     `FREESURFER_HOME` environment variable is set.
#'   \item Large datasets may require substantial memory. Consider adjusting
#'     `chunk_size` and `n_cores` based on your system specifications.
#'   \item For reproducibility, always specify a `seed`.
#' }
#'
#' @seealso [single_lmm()] for single-vertex modeling,
#' `vignette("03-run-vw-lmm", package = "verywise")` for detailed
#' usage examples.
#'
# #' @examplesIf file.exists("path/to/phentype/data")
# #' # Basic cortical thickness analysis
# #' results <- run_vw_lmm(
# #'   formula = vw_thickness ~ age + sex + site + (1|participant_id),
# #'   pheno = "path/to/phentype/data", # or data.frame object
# #'   subj_dir = "/path/to/freesurfer/subjects",
# #'   outp_dir = "/path/to/output",
# #'   hemi = "lh",
# #'   n_cores = 4)
#'
#' @author Serena Defina, 2024.
#'
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
  weights = NULL,
  lmm_control = lme4::lmerControl(calc.derivs=FALSE),
  REML = TRUE,
  # Reproducibility and parallel processing
  seed = 3108,
  n_cores = 1,
  chunk_size = 1000,
  # Cluster estimation
  mtc = c('mcz', 'fdr'),
  FS_HOME = Sys.getenv("FREESURFER_HOME"),
  fwhm = 10,
  mcz_thr = 30,
  cwp_thr = 0.025,
  # Output control
  save_optional_cluster_info = FALSE,
  save_ss = FALSE,
  save_residuals = FALSE,
  save_cov = NULL,
  verbose = TRUE) {
  
  vw_init_message('Linear mixed model', verbose = verbose)

  hemi <- match.arg(hemi)
  measure <- check_formula(formula)
  model_desc <- paste(as.character(formula)[c(1,3)], collapse = ' ') # Only lhs

  vw_message('* Outcome: {.val2 {outcome_name(hemi, measure)}}', verbose = verbose)
  vw_message('* Model:   {.val2 {model_desc}}', verbose = verbose)

  # Check user input ===========================================================
  vw_message("User input validation and set-up", type='step', verbose = verbose)

  # Path specifications
  if (verbose) cli::cli_progress_step('Input and output paths', spinner=TRUE)
  
  subj_dir <- check_path(subj_dir)
  outp_dir <- check_path(outp_dir, create_if_not = TRUE)

  ss_file <- check_ss_exists(subj_dir, hemi, measure, fs_template)

  # Numeric input
  if (verbose) cli::cli_progress_step('Check settings and prepare environment', spinner=TRUE)

  check_numeric_param(seed, integer = TRUE, lower = 0)
  check_numeric_param(chunk_size, integer = TRUE, lower = 1, upper = 5000) # for memory safety

  mtc <- match.arg(mtc) 

  if (mtc == 'mcz') {
    check_numeric_param(fwhm, lower = 1, upper = 30)
    check_numeric_param(mcz_thr, set=c(13, 20, 23, 30, 33, 40))
    check_numeric_param(cwp_thr, set=c(0.025, 0.05))
    check_freesurfer_setup(FS_HOME, verbose = verbose)
  }
  
  n_cores <- check_cores(n_cores)

  if (verbose) cli::cli_progress_done()

  # Avoid bigstatsr warning about lost precision (float vs. double)
  old_opts <- options(bigstatsr.downcast.warning = FALSE)
  on.exit(options(old_opts), add = TRUE)

  # Esure reproducible seeds in parallel settings
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  start.time <- Sys.time()
  
  # Read phenotype data (if not already loaded) ================================

  if (verbose) cli::cli_progress_step('Load and transform phenotype', spinner=TRUE)

  if (is.character(pheno)) pheno <- load_pheno_file(pheno)

  # Transform to list of dataframes (imputed and single datasets alike)
  data_list <- imp2list(pheno); rm(pheno)

  # Check that the data_list is not empty, it contains data.frames of the same
  # dims, and that "folder_id" and all variables specified in the formula are
  # present in the data
  check_data_list(data_list, folder_id, formula)

  # Extract first dataset
  data1 <- data_list[[1]]

  if (verbose) cli::cli_progress_done()

  vw_message(" * Phenotype: {.val {length(data_list)}} dataset{?s} of dimensions
             {.val2 { nrow(data1) }} x {.val2 { ncol(data1) }}.", verbose = verbose)
  
  # Unpack model ===============================================================
  check_weights(weights, data1)

  fixed_terms <- unpack_formula(formula, data1)

  # Check that the stacks are not overwritten by mistake and
  # Save the stack names (i.e. fixed terms) to a lookup file
  check_stack_file(fixed_terms, outp_dir)

  folder_ids <- data1[, folder_id, drop=TRUE] # ensure this is always a character vector 

  cov_terms <- check_cov_spec(save_cov, fixed_terms)

  # Read and clean vertex data =================================================
  
  vw_message("Brain data processing", type='step', verbose = verbose)
  
  if (is.character(save_ss)) {
    ss_dir <- check_path(save_ss, create_if_not = TRUE)
    save_ss <- TRUE
  } else {
    ss_dir <- file.path(outp_dir, 'ss')
    if (!save_ss) on.exit(unlink(ss_dir, recursive = TRUE), add = TRUE)
  }

  if (is.null(ss_file)) {

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
      verbose = verbose)

  } else {

    ss <- subset_supersubject(
      supsubj_dir = subj_dir,
      supsubj_file = ss_file,
      folder_ids = folder_ids,
      new_supsubj_dir = ss_dir,
      fs_template = fs_template,
      n_cores = n_cores,
      save_rds = save_ss,
      error_cutoff = tolerate_surf_not_found,
      verbose = verbose)
  }

  if (verbose) cli::cli_progress_step('Clean and chunk super-subject matrix', spinner=TRUE)

  # Cortical mask
  is_cortex <- mask_cortex(hemi = hemi, fs_template = fs_template)

  # Additionally check that there are no vertices that contain any 0s
  problem_verts <- fbm_col_has_0(ss, n_cores = 1L, verbose = verbose)

  good_verts <- which(!problem_verts & is_cortex); rm(problem_verts)

  # Ensure phenotype and ss row order matches ==================================
  data_list <- check_row_match(ss_file = ss$bk, pheno = data_list, 
                               folder_ids = folder_ids)

  # Prepare chunk sequence =====================================================
  chunk_seq <- make_chunk_sequence(good_verts, chunk_size = chunk_size)

  # Number of vertices
  vw_n <- length(is_cortex); rm(is_cortex)
  # Number of terms (excluding random terms)
  fe_n <- length(fixed_terms)
  # Number of participants*timepoint (long format)
  n_obs <- nrow(data_list[[1]])
  # Number of (imputed) datasets
  m <- length(data_list)

  if (verbose) cli::cli_progress_done()

  vw_message(c(">" = "Ready to run {.val {length(good_verts)}} models 
                     (split in {.val2 {length(chunk_seq)}} chunks)."),
             verbose = verbose)

  vw_message("Statistical model fitting", type='step', verbose = verbose)
  
  # Cache the model frame: `refit_lmm` uses an "update"-based workflow to minimize
  # repeated parsing and model construction overhead
  model_template <- precompile_model(formula = formula, data_list = data_list, 
    tmp_y = ss[, good_verts[1]], measure = measure, weights = weights,
    lmm_control = lmm_control, REML = REML, verbose = verbose)
  
  n_random_groups <- summary(model_template[[1]])$ngrps; storage.mode(n_random_groups) <- "integer"
  n_obs_effective <- sapply(model_template, stats::nobs)  
  
  if (length(unique(n_obs_effective)) > 1) {
    vw_error('Missing pattern is inconsistent across imputations.')
  } else {
    n_obs_effective <- n_obs_effective[1]
  }
  
  vw_message(c("i"= "model includes {.val2 {fe_n}} fixed parameters and {.val2 {n_random_groups}} groups"), 
    verbose = verbose)

  # Prepare FBM output =========================================================

  result_path <- file.path(outp_dir, paste(hemi, measure, sep = "."))

  # Temporary output matrices # note: default single precision (32 bits)
  # Coefficients, SE and P-values 
  c_vw <- build_output_fbm(result_path, "coef", nrow = fe_n, ncol = vw_n, verbose = verbose) 
  s_vw <- build_output_fbm(result_path, "se",   nrow = fe_n, ncol = vw_n) 
  p_vw <- build_output_fbm(result_path, "p",    nrow = fe_n, ncol = vw_n)
  # Residuals
  r_vw <- build_output_fbm(result_path, "resid", nrow = n_obs_effective, ncol = vw_n)
  # Fit statistics: singular_fits, aic, r2_conditional, r2_marginal, icc by group
  f_vw <- build_output_fbm(result_path, "mfit", nrow = (4L + length(n_random_groups)), ncol = vw_n)
  # Covariance between two terms
  if (!is.null(cov_terms)) {
    cov_vw <-  build_output_fbm(result_path, "cov", nrow = 1, ncol = vw_n)
  }
  
  log_file <- paste0(result_path, ".issues.log") # Log model fitting issues

  # Parallel analyses ==========================================================

  progress_file <- paste0(result_path, ".progress.log")
  on.exit(if (file.exists(progress_file)) file.remove(progress_file), add = TRUE)

  # Progress bar setup # note progressr only works with doFuture not doParallel
  if (verbose) {
    cli::cli_progress_step("Fitting linear mixed models... this may take some time", spinner=TRUE)
    vw_message(c("i"="Check the {.file {basename(progress_file)}} file for updates."))
  }
  with_parallel(n_cores = n_cores, 
    seed = seed,
    verbose = verbose, 
    expr = {
      foreach::foreach(chunk = chunk_seq, 
        .packages = c("bigstatsr"), 
        .export = c("refit_lmm", "vw_pool",
                    "init_progress_tracker", "update_progress_tracker")
    ) %dopar% { # Only parallel if n_cores > 1

      # Progress updates
      progress_tracker <- init_progress_tracker(chunk, chunk_seq, 
        progress_file = progress_file, verbose=verbose)

      for (v in chunk) {

        update_progress_tracker(v, progress_tracker, 
          progress_file = progress_file, verbose = verbose)

        # NOTE: ss does not need to be copied to each worker with doParallel
        vertex <- ss[, v]

        # Loop through imputed datasets and run analyses
        out_stats <- lapply(model_template, refit_lmm, y = vertex, 
          cov_eff = cov_terms)

        # Pool results
        pooled_stats <- vw_pool(out_stats, m = m, n_terms = fe_n, 
          pvalue_method="t-as-z", cov_eff = cov_terms)

        # Log errors (if any)
        if (is.character(pooled_stats)) {
          cat(paste0(v, "\t", pooled_stats, "\n"), file = log_file, append = TRUE)
          # & skip to the next value of v
          next
        }

        # Log warnings (if any)
        if (pooled_stats$warning != "") {
          cat(paste0(v, "\t", pooled_stats$warning, "\n"), file = log_file, append = TRUE)
        }

        # Write results to their respective FBM
        c_vw[, v] <- pooled_stats$coef
        s_vw[, v] <- pooled_stats$se
        p_vw[, v] <- pooled_stats$p # -1 * log10(pooled_stats$p) # convert later
        f_vw[, v] <- pooled_stats$mfit
        r_vw[, v] <- pooled_stats$resid
        if (!is.null(cov_terms)) cov_vw[, v] <-  pooled_stats$cov
      }
    }
  })

  if (verbose) cli::cli_progress_done()
  
  # "coefficients", "standard_errors", "p_values", "fit_statistics", "residuals"
  out <- list(coef = c_vw, se = s_vw, p = p_vw, mfit = f_vw, resid = r_vw)
  if (!is.null(cov_terms)) out[['cov']] <- cov_vw

  # Post-processing ============================================================
  vw_message("Post-processing", type='step', verbose = verbose)

  pval_trans <- switch(mtc, 
    fdr = 'fdr',
    mcz = '-log10p')

  # Save model statistics into separate .mgh files
  convert_to_mgh(out, result_path,
                 fixed_terms = fixed_terms,
                 random_terms = names(n_random_groups),
                 stat_names = c(names(out), pval_trans), # do not save resid as mgh, leave as matrix
                 save_resid = save_residuals,
                 verbose = verbose)
  
  # clusters
  ct_vw <- NULL 
  
  if (mtc == 'mcz') {
    # Estimate full-width half maximum (using FreeSurfer)
    fwhm <- estimate_fwhm(result_path = result_path,
                          hemi = hemi,
                          mask = good_verts,
                          fs_template = fs_template)

    # Clamp fwhm to [1, 30]
    fwhm_clamped <- min(max(fwhm, 1), 30)

    if (fwhm != fwhm_clamped) {
      direction <- if (fwhm > 30) "high. Reduced to 30." else "low. Increased to 1."
      vw_message("! estimated smoothness is {.val {fwhm}}, which is really {direction}",
                verbose = verbose)
      fwhm <- fwhm_clamped
    } else {
      vw_message("Estimated smoothness = {.val {fwhm}}", type = 'note', verbose = verbose)
    }

    # Apply cluster-wise correction (using FreeSurfer) ===========================
    # vw_message("Clusterwise correction...", verbose = verbose)

    if (verbose) cli::cli_progress_step("Clusterwise correction", spinner=TRUE)

    for (stack_n in seq_along(fixed_terms)){
      stack_path <- paste0(result_path, ".stack", stack_n)
      fs_verbosity <- FALSE # if(stack_n == 1) verbose else FALSE

      ocn <- compute_clusters(stack_path = stack_path,
                      hemi = hemi,
                      fwhm = fwhm,
                      FS_HOME = FS_HOME,
                      mcz_thr = mcz_thr,
                      cwp_thr = cwp_thr,
                      fs_template = fs_template,
                      full_surfcluster_output = save_optional_cluster_info,
                      mask = paste0(result_path, ".finalMask.mgh"),
                      verbose = fs_verbosity)
      
      if (is.null(ocn)) break # did not compute clusters

      # else 
      if (is.null(ct_vw)) {
        # create cluster storage (once)
        ct_vw <- build_output_fbm(result_path, 'clust', nrow = fe_n, ncol = vw_n) 
      }

      ct_vw[stack_n, ] <- ocn
    }
    if (verbose) cli::cli_progress_done()
  }

  file.remove(paste(result_path, "residuals.mgh", sep = "."))
  if (save_residuals) {
    r_vw$save()
    vw_message("Residual matrix saved to {.file {r_vw$rds}}",
               verbose = verbose, type = 'note')
  }

  # Print summary 
  model_fit_summary <- vw_summarize_model_fit(fitstats = out$mfit, 
    random_terms = names(n_random_groups), verbose = verbose)
  
  if (!is.null(ct_vw)) {
    out[['clust']] <- ct_vw

    model_res_summary <- vw_summarize_model_clusters(coef = out$coef, clust = out$clust, 
      term_names = fixed_terms, result_path = result_path, verbose = verbose)
    
    cluster_correction_info <- list(
      'fwhm' = fwhm,
      'mcz_threshold' = mcz_thr,
      'cwp_threshold' = cwp_thr)
        
    files_to_remove <- c(
      paste0(result_path, ".clust.bk"),
      list.files(path = outp_dir,
        pattern = paste0("^", basename(result_path), ".*\\.cluster\\.summary$|^",
                              basename(result_path), ".*\\.-log10p\\.mgh$"), 
        recursive = TRUE, full.names = TRUE)
    )
    if (!save_optional_cluster_info) {
    files_to_remove <- c(files_to_remove, paste0(result_path, c('.finalMask.mgh','.inputMask.mgh'))) 
      # paste0(result_path, c('.fwhm.dat', '.finalMask.mgh'))) 
    }

    on.exit(file.remove(files_to_remove), add = TRUE)
    
  } else {
    model_res_summary <- vw_summarize_model_est(coef = out$coef, term_names = fixed_terms, verbose = verbose)

    cluster_correction_info <- ''
  }

  end.time <- Sys.time()

  yaml::write_yaml(
    list(
      'model'=model_desc,
      'n_datasets'=m,
      'n_observations'=n_obs,
      'n_observations_effective'=n_obs_effective,
      'n_groups'=as.list(n_random_groups),
      'n_vertices'=as.integer(count_vertices(fs_template)),
      'n_vertices_effective'=length(good_verts),
      'model_fit' = model_fit_summary,
      'results'= model_res_summary,
      'cluster_correction' = cluster_correction_info,
      'covariance_between_terms'=paste(cov_terms, collapse='and'),
      'random_seed' = as.integer(seed),
      'date'=as.character(Sys.Date()),
      'computation_time_min'=difftime(end.time, start.time, units = "mins"),
      'verywise_version'=as.character(utils::packageVersion('verywise'))),
    file=paste0(result_path, ".model.summary.yml"), column.major = FALSE)

  vw_message("Done! :)", type='step', verbose = verbose)

  return(out)
}
