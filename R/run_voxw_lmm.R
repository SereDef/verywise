#' @title
#' Run voxel-wise linear mixed model using \code{lme4::lmer()}
#'
#' @description
#' This is is the main function for conducting voxel-wise linear mixed model
#' analyses on brain morphology. It will first check use inputs, prepare
#' the phenotype data(list) and run a linear mixed model at each voxel using
#' the [single_lmm()] function.
#'
#' The function supports analysis of both single and multiple imputed datasets.
#' --TODO-- It also automatically handles roi masking, and provides cluster-wise
#' correction for multiple testing using cluster size permutation.
#'
#' @param formula A model formula object. This should specify a linear mixed
#'   model \code{lme4} syntax. The outcome variable should be one of the
#'   supported brain surface metrics (see Details). Example:
#'   \code{vw_value ~ age * sex + site + (1|participant_id)}.
#' @param pheno Either a \code{data.frame}/\code{tibble} containing the
#'   "phenotype" data (i.e., already loaded in the global environment), or a
#'   string specifying the file path to phenotype data. Supported file formats:
#'   .rds, .csv, .txt, .sav (SPSS).
#'   The data should be in \strong{long} format and it should contain all the
#'   variables specified in the left-hand side of the \code{formula} (i.e., after the `~`) 
#'   --TODO: no checking done so far-- plus the \code{obs_id} column.
#' @param ss_file A path to the super-subject matrix file (observations by voxels). The 
#'   rows are assumed in the same order as the phenotype. .csv formats are currently 
#'   supported.
#' @param outp_dir Character string specifying the output directory for results.
#'   If \code{NULL} (default), creates a "results" sub-directory in the
#'   current working directory (not recommended).
#' @param brain_template Character string specifying the brain template for
#'   voxel registration. TMP: number of voxels. Options: --TODO--
#' @param apply_mask Logical vector for ROI masking --TODO--
#' @param weights Optional string or numeric vector of weights for the linear mixed model.
#'   You can use this argument to specify inverse-probability weights. If this
#'   is a string, the function look for a column with that name in the phenotype
#'   data. Note that these are not normalized or standardized in any way.
#'   Default: \code{NULL} (no weights).
#' @param REML Logical specifying whether to optimize the REML criterion (as opposed to
#'   the log-likelihood). Default: TRUE. Use `REML = FALSE` if you intend to do model
#'   comparison (using AIC output).
#' @param lmm_control Optional list (of correct class, resulting from
#'   \code{lmerControl()} containing control parameters to be passed to
#'   \code{lme4::lmer()} (e.g. optimizer choice, convergence criteria,
#'   see the \code{?lmerControl} documentation for details.
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
#' @param save_residuals Logical indicating whether to save the residuals.csv
#'   file. Default: \code{FALSE}.
#' @param verbose Logical indicating whether to display progress messages.
#'   Default: \code{TRUE}.
#'
#' \strong{Statistical Approach:}
#' The function uses \code{lme4::lmer()} for mixed-effects modeling, enabling
#' analysis of longitudinal and hierarchical data. P-values are computed using
#' the t-as-z approximation, with cluster-wise correction applied using
#' --TODO--.
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
#' Parallel processing of multiple models should be handled by the user 
#' (e.g., using SLURM job arrays or similar, see vignette on parallelisation).
#' Within a single run_voxw_lmm call, voxels are divided into chunks of size
#' \code{chunk_size} and processed in parallel across \code{n_cores} workers
#' (when \code{n_cores > 1}). When multiple imputed datasets are present,
#' these are processed sequentially within each voxel.
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
#' Results are saved in --TODO-- format for visualization
#' with --TODO--
#'
#' @return A list of file-backed matrices (\code{bigstatsr::FBM} objects)
#' containing pooled coefficients, SEs, t- and p- values and residuals.
#' Results are also automatically saved to disk in --TODO-- format.
#'
#' @note
#' \itemize{
#'   \item Large datasets may require substantial memory. Consider adjusting
#'     \code{chunk_size} and \code{n_cores} based on your system specifications.
#'   \item For reproducibility, always specify a \code{seed}.
#' }
#'
#' @seealso
#' [single_lmm()] for single-voxel modeling
#'
#'
#' @author Serena Defina, 2026.
#'
#' @export
#'
run_voxw_lmm <- function(
  # Basic settings
  formula,
  pheno,
  ss_file,
  outp_dir = NULL,
  # Brain data processing
  brain_template = NULL,
  apply_mask = NULL,
  # Modeling settings
  weights = NULL,
  REML = TRUE,
  lmm_control = lme4::lmerControl(),
  # Reproducibility and parallel processing
  seed = 3108,
  n_cores = 1,
  chunk_size = 1000,
  # Cluster estimation TODO
  # Output control
  save_residuals = FALSE,
  verbose = TRUE) {

  vw_init_message('Linear mixed model', verbose = verbose)

  require_packages('bigreadr', call_fn = 'run_voxw_lmm')

  measure <- check_formula(formula, measure_control = FALSE)
  model_desc <- paste(as.character(formula)[c(1,3)], collapse = ' ') # Only lhs

  vw_message('* Outcome: {.val2 voxel values ({measure})}', verbose = verbose)
  vw_message('* Model:   {.val2 {model_desc}}', verbose = verbose)

  # Check user input ===========================================================
  vw_message("User input validation and set-up", type='step', verbose = verbose)

  if (verbose) cli::cli_progress_step('Check settings and prepare environment', spinner=TRUE)

  outp_dir <- check_path(outp_dir, create_if_not = TRUE)

  check_numeric_param(brain_template, integer = TRUE, lower = 1)

  check_numeric_param(seed, integer = TRUE, lower = 0)
  check_numeric_param(chunk_size, integer = TRUE, lower = 1,
                      upper = 5000) # for memory safety
  
  n_cores <- check_cores(n_cores)

  # Other set-up stuff =========================================================

  # Avoid bigstatsr wanrining about lost precision (float vs. double)
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

  # Extract first dataset
  data1 <- data_list[[1]]

  if (verbose) cli::cli_progress_done()

  vw_message(" * Phenotype: {.val {length(data_list)}} dataset{?s} of dimensions
             {.val2 { nrow(data1) }} x {.val2 { ncol(data1) }}.", verbose = verbose)

  check_weights(weights, data1)

  fixed_terms <- unpack_formula(formula, data1)

  # Read and clean vertex data =================================================

  vw_message("Brain data processing", type='step', verbose = verbose)
  vw_message("* reading super-subject file from {.path {ss_file}}", verbose = verbose)

  ss <- bigstatsr::big_read(ss_file, type = 'float', select = 1:brain_template,
                            backingfile = file.path(outp_dir, "ss"))
  
  on.exit(unlink(file.path(outp_dir, "ss.bk")))

  if (verbose) cli::cli_progress_step('Clean and chunk super-subject matrix', spinner=TRUE)

  vw_message(c("!" = "Note that the brain data values are checked in this pipeline."),
             verbose = verbose)
  good_voxels <- 1:brain_template # TMP

  # Ensure phenotype and ss row order matches ==================================
  vw_message(c("!" = "Note that the match between rows in the phenotype vs. brain data is not checked in this pipeline."),
             verbose = verbose)

  data1 <- data_list[[1]]
  
  # Prepare chunk sequence =====================================================
  chunk_seq <- make_chunk_sequence(good_voxels, chunk_size = chunk_size)

  # Number of vertices
  vw_n <- length(good_voxels)
  # Number of terms (excluding random terms)
  fe_n <- length(fixed_terms)
  # Number of participants*timepoint (long format)
  n_obs <- nrow(data1)
  # Number of (imputed) datasets
  m <- length(data_list)

  if (verbose) cli::cli_progress_done()

  vw_message(c(">" = "Ready to run {.val {length(good_voxels)}} models 
                     (split in {.val2 {length(chunk_seq)}} chunks)."),
             verbose = verbose)

  vw_message("Statistical model fitting", type='step', verbose = verbose)
  
  # Cache the model frame: `refit_lmm` uses an "update"-based workflow to minimize
  # repeated parsing and model construction overhead
  model_template <- precompile_model(formula = formula, data_list = data_list, 
    tmp_y = ss[, 1], measure = measure, weights = weights,
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

  result_path <- file.path(outp_dir, measure)

  # Temporary output matrices # note: default single precision (32 bits)
  # Coefficients, SE and P-values 
  c_vw <- build_output_fbm(result_path, "coef", nrow = fe_n, ncol = vw_n, verbose = verbose) 
  s_vw <- build_output_fbm(result_path, "se",   nrow = fe_n, ncol = vw_n) 
  p_vw <- build_output_fbm(result_path, "p",    nrow = fe_n, ncol = vw_n)
  # Residuals
  r_vw <- build_output_fbm(result_path, "resid", nrow = n_obs_effective, ncol = vw_n)
  # Fit statistics: singular_fits, aic, r2_conditional, r2_marginal, icc by group
  f_vw <- build_output_fbm(result_path, "mfit", nrow = (4L + length(n_random_groups)), ncol = vw_n)

  log_file <- paste0(result_path, ".issues.log") # Log model fitting issues

  # Parallel analyses ==========================================================

  progress_file <- paste0(result_path, ".progress.log")
  on.exit(if(file.exists(progress_file)) file.remove(progress_file), add = TRUE)

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
        .export = c("single_lmm", "vw_pool",
                    "init_progress_tracker", "update_progress_tracker")
    ) %dopar% { # Only parallel if n_cores > 1

      # Progress updates
      progress_tracker <- init_progress_tracker(chunk, chunk_seq, 
        progress_file = progress_file, verbose=verbose)

      for (v in chunk) {

        update_progress_tracker(v, progress_tracker, 
          progress_file = progress_file, verbose = verbose)

        # NOTE: ss does not need to be copied to each worker with doParallel
        voxel <- ss[, v]

        # Loop through imputed datasets and run analyses
        out_stats <- lapply(model_template, refit_lmm, y = voxel, 
          cov_eff = NULL)
        # out_stats <- lapply(data_list, single_lmm,
        #                     y = voxel,
        #                     y_name = paste0("vw_", measure),
        #                     model_formula = formula,
        #                     REML = REML, 
        #                     lmm_control = lmm_control,
        #                     weights = weights)

        # Pool results
        pooled_stats <- vw_pool(out_stats, m = m, n_terms = fe_n, 
          pvalue_method="t-as-z", cov_eff = NULL)

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
        p_vw[, v] <- pooled_stats$p
        f_vw[, v] <- pooled_stats$mfit
        r_vw[, v] <- pooled_stats$resid
      }
    }
  })

  if (verbose) cli::cli_progress_done()

  out <- list(coef = c_vw, se = s_vw, p = p_vw, mfit = f_vw, resid = r_vw)

  # Post-processing ==========================================================
  vw_message("Post-processing", type='step', verbose = verbose)

  # Save model statistics into separate .mgh files
  convert_to_mgh(out, result_path,
                 fixed_terms = fixed_terms,
                 random_terms = names(n_random_groups),
                 stat_names = c('coef','se','p','fdr','mfit'),
                 verbose = verbose)

  # TODO: other multiple testing correction? =======================================
  model_fit_summary <- vw_summarize_model_fit(fitstats = out$mfit, 
    random_terms = names(n_random_groups), verbose = verbose)
  
  model_res_summary <- vw_summarize_model_est(coef = out$coef, term_names = fixed_terms, verbose = verbose)

  end.time <- Sys.time()

  yaml::write_yaml(
    list(
      'model'=model_desc,
      'n_datasets'=m,
      'n_observations'=n_obs,
      'n_observations_effective'=n_obs_effective,
      'n_groups'=as.list(n_random_groups),
      'n_voxels'=brain_template,
      'n_voxels_effective'=length(good_voxels),
      'model_fit' = model_fit_summary,
      'results'= model_res_summary,
      'random_seed' = as.integer(seed),
      'date'=as.character(Sys.Date()),
      'computation_time_min'=difftime(end.time, start.time, units = "mins"),
      'verywise_version'=as.character(utils::packageVersion('verywise'))),
    file=paste0(result_path, ".model.summary.yml"), column.major = FALSE)

  vw_message("Done! :)", type='step', verbose = verbose)

  return(out)
}
