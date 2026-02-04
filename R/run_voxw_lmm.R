#' @title
#' Run voxel-wise linear mixed model using \code{lme4::lmer()}
#'
#' @description
#' This is is the main function for conducting voxel-wise linear mixed model
#' analyses on brain morphology. It will first check use inputs, prepare
#' the phenotype data(list) and run a linear mixed model at each voxel using
#' the \code{\link{single_lmm}} function.
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
#' \code{\link{single_lmm}} for single-voxel modeling
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
  lmm_control = lme4::lmerControl(),
  # Reproducibility and parallel processing
  seed = 3108,
  n_cores = 1,
  chunk_size = 1000,
  # Cluster estimation TODO
  # Output control
  save_residuals = FALSE,
  verbose = TRUE) {

  # vw_message(pretty_message(paste(hemi_name, "hemisphere")), verbose = verbose)

  # Check user input ===========================================================
  vw_message("Checking user inputs...", verbose = verbose)

  # measure <- check_formula(formula)

  # subj_dir <- check_path(subj_dir)
  # ss_exists <- check_ss_exists(subj_dir, ss_file)

  outp_dir <- check_path(outp_dir, create_if_not = TRUE)

  # check_freesurfer_setup(FS_HOME, verbose = verbose)

  n_cores <- check_cores(n_cores)

  check_numeric_param(seed, integer = TRUE, lower = 0)
  check_numeric_param(chunk_size, integer = TRUE, lower = 1,
                      upper = 5000) # for memory safety
  # check_numeric_param(fwhm, lower = 1, upper = 30)
  # check_numeric_param(mcz_thr, lower = 0)
  # check_numeric_param(cwp_thr, lower = 0)

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

  # Transform to list of dataframes (imputed and single datasets alike)
  data_list <- imp2list(pheno); rm(pheno)

  # Check that the data_list is not empty, it contains data.frames of the same
  # dims, and that "folder_id" and all variables specified in the formula are
  # present in the data
  # TODO: simplify check 
  # check_data_list(data_list, folder_id, formula)

  # Extract first dataset
  data1 <- data_list[[1]]

  vw_message(" * ", length(data_list), " datasets of dimention: ", nrow(data1),
             " x ", ncol(data1), verbose = verbose)

  check_weights(weights, data1)

  fixed_terms <- get_terms(formula, data1)

  # Check that the stacks are not overwritten by mistake and
  # Save the stack names (i.e. fixed terms) to a lookup file
  # check_stack_file(fixed_terms, outp_dir)

  # Read and clean vertex data =================================================

  # vw_message("Checking and preparing brain surface data...", verbose = verbose)


  vw_message(" * reading super-subject file: ", ss_file, verbose = verbose)

  ss <- bigstatsr::big_read(ss_file, type = 'float', select = 1:brain_template,
                            backingfile = file.path(outp_dir, "ss"))
  
  on.exit(unlink(file.path(outp_dir, "ss.bk")))

  vw_message(" * cleaning super-subject matrix...", verbose = verbose)

  # Cortical mask
  # is_cortex <- mask_cortex(hemi = hemi, fs_template = fs_template)

  # good_verts <- which(!is_cortex);
  good_voxels <- 1:brain_template # TMP

  # Ensure phenotype and ss row order matches ==================================
  # rownames_file <- gsub('.fsaverage\\d*\\.supersubject.bk',
  #                       '.ss.rownames.csv', ss$bk)
  # data_list <- check_row_match(rownames_file = rownames_file,
  #                              data_list = data_list,
  #                              folder_ids = folder_ids)

  data1 <- data_list[[1]]

  # Unpack model ===============================================================
  vw_message("Statistical model preparation...",
             "\n * call: ", deparse(formula), verbose = verbose)

  # Number of vertices
  vw_n <- length(good_voxels)
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
  model_template <- precompile_model(use_model_template = TRUE,
    formula = formula, tmp_data = data1, tmp_y = ss[,1],
    measure = 'value', lmm_control = lmm_control, verbose = verbose)

  # Prepare FBM output =========================================================

  result_path <- outp_dir # file.path(outp_dir, paste(hemi, measure, sep = "."))

  # Temporary output matrices
  res_bk_names <- c("coef", "se", "p", "resid") # "t",
  res_bk_paths <- build_output_bks(result_path, res_bk_names = res_bk_names,
                                   verbose = verbose)
  # These files will be removed by "on.exit" by convert_to_mgh

  fbm_precision <- "float" # single precision – 32 bits

  c_vw <- bigstatsr::FBM(fe_n, vw_n, init = 0, type = fbm_precision,
                         backingfile = res_bk_paths["coef"])  # Coefficients
  s_vw <- bigstatsr::FBM(fe_n, vw_n, init = 0, type = fbm_precision,
                         backingfile = res_bk_paths["se"])    # Standard errors
  # t_vw <- bigstatsr::FBM(fe_n, vw_n, init = 0,
  #                        backingfile = res_bk_paths["t"])     # t values
  p_vw <- bigstatsr::FBM(fe_n, vw_n, init = 1, type = fbm_precision,
                         backingfile = res_bk_paths["p"])     # P values
  r_vw <- bigstatsr::FBM(n_obs, vw_n, init = 0, type = fbm_precision,
                         backingfile = res_bk_paths["resid"]) # Residuals

  log_file <- paste0(result_path, "/issues.log") # Log model fitting issues

  # Prepare chunk sequence =====================================================
  vw_message(" * chunk dataset", verbose = verbose)
  chunk_seq <- make_chunk_sequence(good_voxels, chunk_size = chunk_size)

  # Parallel analyses ==========================================================
  vw_message("Running analyses...\n",
             " * dimentions: ", n_obs, " observations x ", length(good_voxels),
             " (of ", vw_n, " total) vertices.", verbose = verbose)

  progress_file <- paste0(result_path, "/progress.log")
  on.exit(file.remove(progress_file), add = TRUE)


  # Progress bar setup # note progressr only works with doFuture not doParallel
  vw_message(" * fitting linear mixed models...\n",
             "   this may take some time, check the ", basename(progress_file),
             " file for updates.", verbose = verbose)
  
  with_parallel(n_cores = n_cores, 
    progress_file = progress_file,
    seed = seed,
    verbose = verbose, 
    expr = {
      foreach::foreach(chunk = chunk_seq, 
        .packages = c("bigstatsr"), 
        .export = c("single_lmm", "vw_pool",
                    "init_progress_tracker", "update_progress_tracker")
    ) %dopar% { # Only parallel if n_cores > 1

      # Progress updates
      progress_tracker <- init_progress_tracker(chunk, chunk_seq, verbose)

      for (v in chunk) {

        update_progress_tracker(v, progress_tracker, verbose)

        # NOTE: ss does not need to be copied to each worker with doParallel
        voxel <- ss[, v]

        # Loop through imputed datasets and run analyses
        out_stats <- lapply(data_list, single_lmm,
                            y = voxel,
                            formula = formula,
                            model_template = model_template,
                            weights = weights,
                            lmm_control = lmm_control)

        # Pool results
        pooled_stats <- vw_pool(out_stats, m = m, pvalue_method="t-as-z")

        # Log errors (if any)
        if (is.character(pooled_stats)) {
          cat(paste0(v, "\t", pooled_stats, "\n"), file = log_file,
              append = TRUE)
          # & skip to the next value of v
          next
        }

        # Log warnings (if any)
        if (pooled_stats$warning != "") {
          cat(paste0(v, "\t", pooled_stats$warning, "\n"), file = log_file,
              append = TRUE)
        }

        # Write results to their respective FBM
        c_vw[, v] <- pooled_stats$coef
        s_vw[, v] <- pooled_stats$se
        p_vw[, v] <- pooled_stats$p # -1 * log10(pooled_stats$p) # convert later
        r_vw[, v] <- pooled_stats$resid
      }
    }
  })

  out <- list(c_vw, s_vw, p_vw, r_vw) # t_vw,
  # "coefficients", "standard_errors", "t_values", "p_values", "residuals"
  names(out) <- res_bk_names

  # Post-processing ============================================================

  # Save model statistics into separate .mgh files
  # convert_to_mgh(out,
  #                result_path,
  #                stacks = seq_along(fixed_terms),
  #                stat_names = c(res_bk_names,'-log10p'),
  #                verbose = verbose)

  # resid_mgh_path <- paste(result_path, "residuals.mgh", sep = ".")
  # if (!save_residuals) on.exit(file.remove(resid_mgh_path), add = TRUE)

  # # Estimate full-width half maximum (using FreeSurfer) ========================
  # vw_message("Estimating data smoothness for multiple testing correction...",
  #            verbose = verbose)

  # fwhm <- estimate_fwhm(result_path = result_path,
  #                       hemi = hemi,
  #                       mask = good_verts,
  #                       fs_template = fs_template)

  # # Clamp fwhm to [1, 30]
  # fwhm_clamped <- min(max(fwhm, 1), 30)

  # if (fwhm != fwhm_clamped) {
  #   direction <- if (fwhm > 30) "high. Reduced to 30." else "low. Increased to 1."
  #   vw_message(sprintf(" * estimated smoothness is %s, which is really %s",
  #                      fwhm, direction), verbose = verbose)
  #   fwhm <- fwhm_clamped
  # } else {
  #   vw_message(sprintf(" * estimated smoothness = %s", fwhm), verbose = verbose)
  # }

  # # Apply cluster-wise correction (using FreeSurfer) ===========================
  # vw_message("Clusterwise correction...", verbose = verbose)

  # for (stack_n in seq_along(fixed_terms)){
  #   stack_path <- paste0(result_path, ".stack", stack_n)
  #   fs_verbosity <- FALSE # if(stack_n == 1) verbose else FALSE
  #   compute_clusters(stack_path = stack_path,
  #                    hemi = hemi,
  #                    fwhm = fwhm,
  #                    FS_HOME = FS_HOME,
  #                    mcz_thr = mcz_thr,
  #                    cwp_thr = cwp_thr,
  #                    full_surfcluster_output = save_optional_cluster_info,
  #                    mask = paste0(result_path, ".finalMask.mgh"),
  #                    verbose = fs_verbosity)
  # }

  vw_message(pretty_message("All done! :)"))

  return(out)
}
