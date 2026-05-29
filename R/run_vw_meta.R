#' @title
#' Vertex-wise Random-Effects Meta-Analysis Across Studies
#'
#' @description
#' This function runs a vertex-wise random-effects meta-analysis across multiple
#' studies (sites or cohorts). Is expects individual study results as they are
#' outputted but `verywise::run_vw_lmm` or by `QDECR`.
#'
#' @param term Character. Name of the model term to visualize (matches entries in `stack_names.txt`).
#' @param hemi Character string specifying which hemisphere to analyze.
#'   Options: `"lh"` (left hemisphere: default), `"rh"` (right hemisphere).
#' @param measure Character. Surface measure, e.g. `'area'`, `'thickness'`, `'volume'`. Defaults to `'area'`.
#' @param study_names Character vector of study names (must match the length of \code{res_dirs}).
#' @param study_weights Numeric vector of study weights (e.g. sample size)
#' @param res_dirs Character vector. Path to the directories containing vertex-wise result files (`*.mgh`) of each study.
#' @param outp_dir Character string specifying the output directory for results.
#'   If \code{NULL} (default), creates a "verywise_results" sub-directory in the
#'   current working directory (not recommended).
#' @param mtc Character string: multiple testing correction strategy. 
#'   Options: `"fdr"` (False Discovery Rate: default), `"fs"` (FreeSurfer cluster correction).
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
#'  output form \code{mri_surfcluster} call. See \code{\link{compute_clusters}}
#'  for details. Default: \code{FALSE}.
#' @param verbose Logical. verbose execution (default: TRUE)
#'
#' @return A named list with three \code{FBM} objects:
#'   \describe{
#'     \item{coef}{Meta-analytic effect size estimates at each vertex.}
#'     \item{se}{Meta-analytic standard errors at each vertex.}
#'     \item{p}{Meta-analytic p-values at each vertex.}
#'   }
#'   Additionally, these are exported as MGH files in the working directory.
#'
#' @details
#' For each vertex, the function loads the effect size and standard error from each study, computes the variance,
#' and runs a random-effects meta-analysis using \code{\link[metafor]{rma}}. Results are saved as
#' Filebacked Big Matrices (FBM) and as MGH files for downstream neuroimaging analysis.
#'
#' @seealso \code{\link[metafor]{rma}}, \code{\link[bigstatsr]{FBM}}
#'
#' @examplesIf dir.exists("study1/results")
#' run_vw_meta(
#'   term = "age",
#'   hemi = "lh",
#'   measure = "area",
#'   study_names = c("Study1", "Study2"),
#'   res_dirs = c("study1/results", "study2/results")
#' )
#'
#' @export
run_vw_meta <- function(term, 
                        hemi = c('lh','rh'),
                        measure = "area",
                        study_names, 
                        study_weights = NULL,
                        res_dirs,
                        outp_dir = NULL,
                        mtc = 'fdr',
                        fs_template = 'fsaverage',
                        # Reproducibility and parallel processing
                        seed = 3108,
                        n_cores = 1,
                        chunk_size = 1000,
                        # Cluster estimation
                        FS_HOME = Sys.getenv("FREESURFER_HOME"),
                        mcz_thr = 30,
                        cwp_thr = 0.025,
                        save_optional_cluster_info = FALSE,
                        verbose = TRUE) {
  
  vw_init_message('Meta-analysis', verbose = verbose)
  
  # Validate input
  hemi <- match.arg(hemi)
  check_measure(measure)
  n_studies <- length(study_names)

  vw_message('* Outcome: {.val2 {outcome_name(hemi, measure)}}', verbose = verbose)
  vw_message('* Predictor: {.val2 {term}}', verbose = verbose)
  vw_message('* {.strong {n_studies}} studies: {study_names}', verbose = verbose)
  vw_message(' ', verbose = verbose)

  if (verbose) cli::cli_progress_step('User input validation and set-up', spinner=TRUE)
  
  outp_dir <- check_path(outp_dir, create_if_not = TRUE)
  
  if (length(res_dirs) != n_studies) {
    vw_error('The number of result directories ({length(res_dirs)}) does not match the number of studies ({n_studies}).')
  }

  if (!is.null(study_weights) & length(study_weights) != n_studies) {
    vw_error('The number of weights provided ({length(study_weights)}) does not match the number of studies ({n_studies}).')
  }

  if (length(term) == 1) term <- rep(term, n_studies) else if (length(term) != n_studies) {
    vw_error("The number of terms ({length(term)}) must be 1 or match the number of studies ({n_studies}).")
  }

  stack_paths <- vapply(seq_len(n_studies), function(s) {
  
    res_dir <- check_path(res_dirs[s])
    s_term <- term[s]
    
    # Find term stack index
    stack_file <- file.path(res_dir, "stack_names.txt")
    
    if (!file.exists(stack_file)) {
      vw_error("Cannot find `stack_names.txt` in {.file {res_dir}}")
    }
    
    stack_ids <- utils::read.table(stack_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    if (!s_term %in% stack_ids$stack_name) {
      vw_error("Term {.warn {s_term}} not found in {.file {res_dir}/stack_names.txt}")
    }

    s_stack <- paste0('stack', stack_ids[stack_ids$stack_name == s_term, 'stack_number'])
    
    # Return the constructed file path to vapply
    file.path(res_dir, paste(hemi, measure, s_stack, sep = '.'))
    
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)

  n_cores <- check_cores(n_cores)
  check_numeric_param(seed, integer = TRUE, lower = 0)
  check_numeric_param(chunk_size, integer = TRUE, lower = 1, upper = 5000) # for memory safety

  if (mtc == 'fs') {
    check_numeric_param(mcz_thr, set=c(13, 20, 23, 30, 33, 40))
    check_numeric_param(cwp_thr, set=c(0.025, 0.05))

    check_freesurfer_setup(FS_HOME, verbose = verbose)

    fwhms <- vapply(res_dirs, function(s) {
      as.numeric(readLines(file.path(s, paste(hemi, measure, 'fwhm.dat', sep='.'))))
    }, FUN.VALUE = numeric(1), USE.NAMES = FALSE)

    fwhm <- mean(fwhms)

    vw_message(c("i"="Average smoothness = {.val {fwhm}}", 
                 " "="input values {fwhms}"))
  }

  # Cortical mask
  is_cortex <- mask_cortex(hemi = hemi, fs_template = fs_template)

  # n_verts <- count_vertices(fs_template)
  n_verts <- length(is_cortex)

  if (verbose) cli::cli_progress_done()
  
  # Avoid bigstatsr warning about lost precision (float vs. double)
  old_opts <- options(bigstatsr.downcast.warning = FALSE)
  on.exit(options(old_opts), add = TRUE)

  # Esure reproducible seeds in parallel settings
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  if (verbose) cli::cli_progress_step('Load study estimates', spinner=TRUE)
  
  # Temporary output matrices
  result_path <- file.path(outp_dir, paste(hemi, measure, term[1], sep = "."))

  res_bk_names <- c(paste0("stud.", c("coef", "se")),
                    paste0("meta.", c("coef", "se", "p")))
  res_bk_paths <- build_output_bks(result_path, res_bk_names = res_bk_names,
                                   verbose = verbose)
  
  fbm_precision <- "float" # single precision – 32 bits
  
  # Effect sizes
  ef_vw <- bigstatsr::FBM(n_studies, n_verts, init = NA_real_, type = fbm_precision,
                          backingfile = res_bk_paths["stud.coef"])
  # Variance of effect sizes
  se_vw <- bigstatsr::FBM(n_studies, n_verts, init = NA_real_, type = fbm_precision,
                          backingfile = res_bk_paths["stud.se"])
  
  on.exit(file.remove(paste0(res_bk_paths[c("stud.coef","stud.se")],'.bk')), add = TRUE)

  for (s in seq_len(n_studies)) {

    res_path <- stack_paths[s]

    ef_vw[s, ] <- load.mgh(paste0(res_path, ".coef.mgh"))$x
    se_vw[s, ] <- load.mgh(paste0(res_path, ".se.mgh"))$x

  }

  # Pooled coefficients, standard errors and p-values 
  c_vw <- bigstatsr::FBM(1, n_verts, init = NA_real_, type = fbm_precision,
                         backingfile = res_bk_paths["meta.coef"])
  s_vw <- bigstatsr::FBM(1, n_verts, init = 0, type = fbm_precision,
                         backingfile = res_bk_paths["meta.se"])
  p_vw <- bigstatsr::FBM(1, n_verts, init = 1, type = fbm_precision,
                         backingfile = res_bk_paths["meta.p"])
  
  # Only process vertex where more than one study has estimate
  problem_verts <- fbm_col_has_na(ef_vw, n_cores = 1L, verbose = verbose)

  good_verts <- which(!problem_verts & is_cortex); rm(problem_verts)

  chunk_seq <- make_chunk_sequence(good_verts, chunk_size = chunk_size)

  vw_message(c(">" = "Ready to meta-analyze {.val {length(good_verts)}} vertices 
                     (split in {.val2 {length(chunk_seq)}} chunks)."),
             verbose = verbose)

  if (verbose) cli::cli_progress_step("Meta-analyzing results", spinner=TRUE)
  
  # Parallel analyses ==========================================================

  log_file <- paste0(result_path, ".issues.log") # Log model fitting issues

  progress_file <- paste0(outp_dir, ".progress.log")
  on.exit(if (file.exists(progress_file)) file.remove(progress_file), add = TRUE)

  if (verbose) {
    cli::cli_progress_step("Meta-analyzing results... this may take some time, check the {.file {basename(progress_file)}} file for updates.", 
    spinner=TRUE)
  }

  with_parallel(n_cores = n_cores, 
    seed = seed,
    verbose = verbose, 
    expr = {
      foreach::foreach(chunk = chunk_seq, 
        .packages = c("bigstatsr", "metafor"), 
        .export = c("init_progress_tracker")
    ) %dopar% { # Only parallel if n_cores > 1

      # Progress updates
      progress_tracker <- init_progress_tracker(chunk, chunk_seq, 
        progress_file = progress_file, verbose = verbose)

      for (v in chunk) {

        rma_args <- list(
          yi = ef_vw[, v], 
          sei = se_vw[, v], 
          method = "REML")

        # Only add weights to the argument list if they were provided
        if (!is.null(study_weights)) rma_args$weights <- study_weights

        meta_res <- tryCatch(
          do.call(metafor::rma, rma_args),
          error = function(e) {
            # log error
            cat(paste0(v, "\t", e, "\n"), file = log_file, append = TRUE)

            # return empty 
            list(beta = c(NA_real_), se = NA_real_, pval = NA_real_)
          }
        )
        
           c_vw[, v] <- meta_res$beta[1]
           s_vw[, v] <- meta_res$se
           p_vw[, v] <- meta_res$pval
        }
      }
  })

  if (verbose) cli::cli_progress_done()

  out <- list(coef = c_vw, se = s_vw, p = p_vw)

  # Post-processing ============================================================
  vw_message("Post-processing", type='step', verbose = verbose)

  # Save model statistics into separate .mgh files
  pval_trans <- switch(mtc, 
    fdr = 'fdr',
    fs = '-log10p')
  
  convert_to_mgh(out,
                 result_path,
                 stacks = NULL,
                 stat_names = c("coef", "se", "p", pval_trans),
                 verbose = verbose)

  if (mtc == 'fs') {
    if (verbose) cli::cli_progress_step("Clusterwise correction", spinner=TRUE)
    
    ocn <- compute_clusters(stack_path = result_path,
                      hemi = hemi,
                      fwhm = fwhm,
                      FS_HOME = FS_HOME,
                      mcz_thr = mcz_thr,
                      cwp_thr = cwp_thr,
                      fs_template = fs_template,
                      full_surfcluster_output = save_optional_cluster_info,
                      mask = NULL,
                      verbose = TRUE)
      
    if (!is.null(ocn)) {
      ct_bk_path <- build_output_bks(result_path, res_bk_names = c('meta.clust'), verbose = FALSE)
      ct_vw <- bigstatsr::FBM(1, n_verts, init = NA_real_, 
        type = fbm_precision, backingfile = ct_bk_path["meta.clust"])
      
      ct_vw[1,] <- ocn
    }

    if (verbose) cli::cli_progress_done()
    
    out[['clust']] <- ct_vw

    vw_summarize_model_clusters(coef = out$coef, clust = out$clust, 
      term_names = term[1], verbose = verbose)
  } else {

    # TODO better summary for FDR corrected 
    vw_summarize_model_est(coef = out$coef, term_names = term[1], verbose = verbose)
    
  }

  vw_message("Done! :)", type='step', verbose = verbose)

  return(out)
  }
