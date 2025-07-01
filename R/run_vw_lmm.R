#' @title
#' Run vertex-wise linear mixed model using \code{lme4::lmer()}
#'
#' @description
#' This is is the main function for conducting vertex-wise linear mixed model analyses
#' on brain surface metrics. It will first check use inputs, prepare the phenotype data(list) 
#' and run the linear mixed model in each of the specified hemispheres using the 
#' \code{\link{hemi_vw_lmm}} function. 
#' 
#' The function supports analysis of both single and multiple imputed datasets.
#' It also automatically handles cortical masking, and provides cluster-wise correction
#' for multiple testing using FreeSurfer's Monte Carlo simulation approach.
#'
#' @param formula A model formula object. This should specify a linear mixed model
#'   \code{lme4} syntax. The outcome variable should be one of the supported brain
#'   surface metrics (see Details). Example: 
#'   \code{vw_thickness ~ age * sex + site + (1|participant_id)}.
#' @param pheno Either a `data.frame`/`tibble` containing the "phenotype" data 
#'   (i.e., already loaded in the global environment), or a string specifying the
#'   file path to phenotype data. Supported file formats: .rds, 
#'   .csv, .txt, .sav (SPSS). 
#'   The data should be in **long** format and it should contain all the variables 
#'   specified in the \code{formula} plus the \code{folder_id} column.
#' @param subj_dir Character string specifying the path to FreeSurfer data directory.
#'   Must follow the verywise directory structure (see package vignette for details).
#' @param outp_dir Character string specifying the output directory for results.
#'   If \code{NULL} (default), creates a "verywise_results" subdirectory within 
#'   \code{subj_dir}.
#' @param hemi Character string specifying which hemisphere(s) to analyze.
#'   Options: \code{"both"} (default), \code{"lh"} (left only), \code{"rh"} (right only).
#' @param folder_id Character string specifying the column name in \code{pheno}
#'   that contains subject directory names of the input neuroimaging data 
#'   (e.g., "sub-001_ses-baseline" or "site1/sub-010_ses-F1"). These are expected to be nested 
#'   inside \code{subj_dir}.
#'   Default: \code{"folder_id"}.
#' @param seed Integer specifying the random seed for reproducibility. Default: 3108.
#' @param n_cores Integer specifying the number of CPU cores for parallel processing.
#'   Default: 1.
#' @param chunk_size Integer specifying the number of vertices processed per chunk
#'   in parallel operations. Larger values use more memory but may be faster.
#'   Default: 1000.
#' @param FS_HOME Character string specifying the FreeSurfer home directory.
#'   Defaults to \code{FREESURFER_HOME} environment variable.
#' @param fs_template Character string specifying the FreeSurfer template for
#'   vertex registration. Options:
#'  * \code{"fsaverage"} (default) = 163842 vertices (highest resolution),
#'  * \code{"fsaverage6"} = 40962 vertices,
#'  * \code{"fsaverage5"} = 10242 vertices,
#'  * \code{"fsaverage4"} = 2562 vertices,
#'  * \code{"fsaverage3"} = 642 vertices
#' Note that lower resolutions should be only used to downsample the brain map, for faster
#' model tuning. The final analyses should also run using `fs_template = "fsaverage"`
#' to avoid (small) imprecisions in vertex registration and smoothing.
#' @param apply_cortical_mask Logical indicating whether to exclude non-cortical
#'   vertices from analysis. Default: \code{TRUE} (recommended).
#' @param use_model_template Logical indicating whether to pre-compile the model
#'   template for faster estimation. Default: \code{TRUE} (recommended).
#' @param save_ss Logical indicating whether to save the super-subject matrix as
#'   an .rds file for faster future processing. Default: \code{TRUE} (recommended).
#' @param verbose Logical indicating whether to display progress messages.
#'   Default: \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link{hemi_vw_lmm}}, including e.g.,
#'   \code{fwhm} (smoothing kernel size) and \code{model} (statistical model type).
#'
#' @details
#' \strong{Supported Brain Surface Metrics:}
#' The *outcome* specified in `formula` should be a brain surface metric among:
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
#' the t-as-z approximation, with cluster-wise correction applied using FreeSurfer's
#' Monte Carlo simulation approach.
#' 
#' \strong{Multiple Imputation:}
#' The function automatically detects and handles multiple imputed datasets
#' (created with \code{mice} or similar packages), pooling results according
#' to Rubin's rules.
#' 
#' \strong{Parallel processing:}
#' The \code{verywise} package employs a carefully designed parallelization strategy 
#' to maximize computational efficiency while avoiding the performance penalties 
#' associated with nested parallelization.
#' Therefore, when `hemi = 'both'`, left and right cortical hemispheres are processed sequentially.
#' Within each hemisphere, vertices are divided into chunks of size \code{chunk_size} and processed 
#' in parallel across \code{n_cores} workers (when \code{n_cores > 1}).
#' When multiple imputed datasets are present, they are processed sequentially within each vertex.
#' 
#' Note that, on some systems, implicit parallelism in low-level matrix algebra libraries 
#' (BLAS/LAPACK) can interfere with explicit parallelization. If you feel like processing is 
#' taking too long, we recommend disabling these implicit threading libraries before starting R.
#' For example:
#' \preformatted{
#' export OPENBLAS_NUM_THREADS=1
#' export OMP_NUM_THREADS=1  
#' export MKL_NUM_THREADS=1
#' export VECLIB_MAXIMUM_THREADS=1
#' export NUMEXPR_NUM_THREADS=1
#' }
#' 
#' \strong{Output Files:}
#' Results are saved in FreeSurfer-compatible .mgh format for visualization
#' with [verywiseWIZard](https://github.com/SereDef/verywise-wizard), FreeView or 
#' other neuroimaging software.
#'  
#' @return A list of file-backed matrices containing pooled coefficients, SEs,
#' t- and p- values and residuals.
#'
#' @author Serena Defina, 2024.
#'
#' @export
#'
run_vw_lmm <- function(formula,
                       pheno,
                       subj_dir,
                       outp_dir = NULL,
                       hemi = c("both", "lh","rh"),
                       folder_id = "folder_id",
                       seed = 3108,
                       n_cores = 1,
                       chunk_size = 1000,
                       FS_HOME = Sys.getenv("FREESURFER_HOME"),
                       fs_template = "fsaverage",
                       apply_cortical_mask = TRUE,
                       use_model_template = TRUE,
                       save_ss = TRUE,
                       verbose = TRUE,
                       ...
) {
  # Check user input ===========================================================

  measure <- check_formula(formula)

  check_path(subj_dir)
  if (!is.null(outp_dir)) check_path(outp_dir)

  check_freesurfer_setup(FS_HOME, verbose=verbose)

  n_cores <- check_cores(n_cores)

  # Read phenotype data (if not already loaded) ================================
  vw_message("Checking and preparing phenotype dataset...", verbose=verbose)

  if (is.character(pheno)) pheno <- load_pheno_file(pheno)

  # Transform to list of dataframes (imputed and single datasets alike)
  data_list <- imp2list(pheno)

  # Check that the data list is not empty and that "folder_id" and all variables
  # specified in the formula are present in the data
  check_data_list(data_list, folder_id, formula)

  # Determine hemisphere(s) ====================================================
  hemi <- match.arg(hemi)
  if (hemi == "both") { hemis <- c("lh","rh") } else { hemis <- as.vector(hemi) }

  # Run analyses ===============================================================
  n_cores <- as.integer(n_cores) # make sure n_cores is an integer

  set.seed(seed)

  # Run analysis in each hemisphere (in sequence)
  out <- lapply(hemis, function(h) {

    if (h == "lh") { hemi_name <- "Left"} else { hemi_name <- "Right" }

    vw_message(pretty_message(paste(hemi_name, "hemisphere")), verbose=verbose)

    hemi_vw_lmm(hemi = h,
                formula = formula,
                subj_dir = subj_dir,
                outp_dir = outp_dir,
                data_list = data_list,
                folder_id = folder_id,
                seed = seed,
                n_cores = n_cores,
                chunk_size = chunk_size,
                FS_HOME = FS_HOME,
                fs_template = fs_template,
                apply_cortical_mask = apply_cortical_mask,
                use_model_template = use_model_template,
                save_ss = save_ss,
                verbose = verbose,
                ...)
    })

  names(out) <- hemis

  out
}

#' @title
#' Run vertex-wise linear mixed model in one hemisphere using \code{lme4::lmer()}
#'
#' @description
#' It runs the \code{\link{single_lmm}} function across vertices in a single hemisphere.
#'
#' @inheritParams run_vw_lmm
#' @param data_list : a list of dataframes containing phenotype information
#' (generated by \code{\link{imp2list}})
#' @param fwhm : (default = 10) full-width half maximum value
#' @param mcz_thr : (default = 0.001) numeric value for the Monte Carlo simulation threshold.
#' Any of the following are accepted (equivalent values separate by `/`):
#'  * 13 / 1.3 / 0.05,
#'  * 20 / 2.0 / 0.01,
#'  * 23 / 2.3 / 0.005,
#'  * 30 / 3.0 / 0.001, \* default
#'  * 33 / 3.3 / 0.0005,
#'  * 40 / 4.0 / 0.0001.
#' @param cwp_thr : (default = 0.025, when both hemispheres are ran, else 0.05)
#' the cluster-wise p-value threshold on top of all corrections.
#' @param model : (default = \code{"lme4::lmer"}) # "stats::lm"
#'
#' @return A list of file-backed matrices containing pooled coefficients, SEs,
#' t- and p- values and residuals.
#'
#' @author Serena Defina, 2024.
#'
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#'
#' @export
#'
hemi_vw_lmm <- function(formula,
                        data_list,
                        subj_dir,
                        outp_dir = NULL,
                        FS_HOME = "",
                        folder_id = "folder_id",
                        hemi,
                        fwhm = 10,
                        fs_template = "fsaverage",
                        mcz_thr = 30,
                        cwp_thr = 0.025,
                        seed = 3108,
                        apply_cortical_mask = TRUE,
                        save_ss = TRUE,
                        model = "lme4::lmer", # "stats::lm"
                        use_model_template = TRUE,
                        n_cores,
                        chunk_size = 1000,
                        verbose = TRUE
) {

  # TMP: Assume the brain measure is always the OUTCOME
  measure <- gsub("vw_", "", all.vars(formula)[1])

  if (is.null(outp_dir)) {
    outp_dir <- file.path(subj_dir, "verywise_results")
    dir.create(outp_dir, showWarnings = FALSE)
  }

  # Read and clean vertex data =================================================
  ss_file_name <- file.path(outp_dir, paste(hemi, measure, fs_template,
                            "supersubject.rds", sep="."))

  if (file.exists(ss_file_name)) {
    vw_message("Reading super-subject file from: ", ss_file_name, verbose=verbose)

    ss <- bigstatsr::big_attach(ss_file_name)

  } else {

    ss <- build_supersubject(subj_dir = subj_dir,
                             folder_ids = data_list[[1]][, folder_id],
                             outp_dir = outp_dir,
                             measure = measure,
                             hemi = hemi,
                             n_cores = n_cores,
                             fwhmc = paste0("fwhm", fwhm),
                             fs_template = fs_template,
                             save_rds = save_ss,
    )
  }

  vw_message("Cleaning super-subject matrix...", verbose=verbose)

  # Cortical mask
  is_cortex <- mask_cortex(hemi = hemi, fs_template = fs_template)

  # Additionally check that there are no vertices that contain any 0s. These may be
  # located at the edge of the cortical map and are potentially problematic
  problem_verts <- fbm_col_has_0(ss)
  if (sum(problem_verts) > 0) {
    vw_message("Ignoring ", sum(problem_verts)," vertices that contained 0 values.",
               verbose=verbose)
  }

  good_verts <- which(!problem_verts & is_cortex)

  # Unpack model ===============================================================
  vw_message("Statistical model preparation...", verbose=verbose)

  model_info <- get_terms(formula, data_list)

  fixed_terms <- model_info[[1]]

  # Number of vertices
  vw_n <- length(is_cortex)
  # Number of terms (excluding random terms)
  fe_n <- length(fixed_terms)
  # Number of participants*timepoint (long format)
  n_obs <- nrow(model_info[[2]])
  # Number of (imputed) datasets
  m <- length(data_list)

  if (use_model_template) {
    # "Pre-compile the model"
    # cache the model frame to avoid re-generating it them each time
    # single_lmm can leverage an "update"-based workflow to minimize repeated parsing
    # and model construction overhead
    vw_message(" * construct model template", verbose=verbose)
    tmp_data <- data_list[[1]]
    tmp_data[paste0('vw_', measure)] <- ss[,1]  # dummy outcome

    # Fit model once
    model_template <- lme4::lmer(formula = formula,
                                 data = tmp_data)
  } else {
    model_template <- NULL
  }

  # Prepare FBM output =========================================================
  vw_message(" * generate file-backed output containers", verbose=verbose)

  res_bk_names <- c("coef", "se", "t", "p", "resid")
  res_bk_paths <- file.path(
    outp_dir,
    paste(hemi, measure, res_bk_names, sep = ".")
  )

  # TMP: remove files if they exist
  for (bk in res_bk_paths) {
    # warning("Overwriting results")
    if (file.exists(paste0(bk,".bk"))) file.remove(paste0(bk,".bk"))
  }
  # Remove temporary matrices files after computations are done
  on.exit({
    for (bk in res_bk_paths) {
      if (file.exists(paste0(bk,".bk"))) file.remove(paste0(bk,".bk"))
    }
  })

  # Coefficients
  c_vw <- bigstatsr::FBM(fe_n, vw_n, init = 0, backingfile = res_bk_paths[1])
  # Standard errors
  s_vw <- bigstatsr::FBM(fe_n, vw_n, init = 0, backingfile = res_bk_paths[2])
  # t values
  t_vw <- bigstatsr::FBM(fe_n, vw_n, init = 0, backingfile = res_bk_paths[3])
  # P values
  p_vw <- bigstatsr::FBM(fe_n, vw_n, init = 1, backingfile = res_bk_paths[4])
  # Residuals
  r_vw <- bigstatsr::FBM(n_obs, vw_n, init = 0, backingfile = res_bk_paths[5])

  # Prepare chunk sequence =====================================================
  vw_message(" * chunk dataset", verbose=verbose)
  chunk_seq <- make_chunk_sequence(good_verts, chunk_size=chunk_size)

  # Parallel analyses ==========================================================
  vw_message("Running analyses...", verbose=verbose)
  vw_message("* dimentions: ", n_obs, " observations x ",
             length(good_verts), " (of ", vw_n, " total) vertices.", verbose=verbose)

  # if (requireNamespace("progressr", quietly = TRUE)) {
  #   # progressr::handlers(global = TRUE)
  #   p <- progressr::progressor(along = good_verts) #1:nrow(chunk_seq)
  # }
  # log_file <- file.path(outp_dir, paste0(hemi,".", measure,"model.log"))

  set.seed(seed)
  if (n_cores > 1 ) vw_message("* preparing cluster of ", n_cores, " workers...",
                               verbose=verbose)
  # Set up parallel processing
  cluster <- parallel::makeCluster(n_cores, type = "FORK")
  # TODO: would not work on Windows anyway till freesurfer dependency is needed
                                   # type = ifelse(.Platform$OS.type == "unix",
                                   #               "FORK", "PSOCK"))
  doParallel::registerDoParallel(cluster)
  on.exit({
    message(pretty_message("All done! :)"))
    parallel::stopCluster(cluster)
  }, add = TRUE)

  chunk = NULL
  foreach::foreach(chunk = chunk_seq,
                   .packages = c("bigstatsr"),
                   .export = c("single_lmm", "vw_pool", "vw_message")) %dopar% {

    worker_id <- Sys.getpid() # TMP for debugging

    for (v in chunk) {

      vertex <- ss[, v] # ss should not get copied n_cores times because it is a memory map

      # Loop through imputed datasets and run analyses
      out_stats <- lapply(data_list, single_lmm,
                          y = vertex,
                          formula = formula,
                          model_template = model_template,
                          pvalues = (m == 1))

      # Pool results
      pooled_stats <- vw_pool(out_stats, m = m)

      # Write results to their respective FBM
      c_vw[, v] <- pooled_stats$coef
      s_vw[, v] <- pooled_stats$se
      t_vw[, v] <- pooled_stats$t
      p_vw[, v] <- -1 * log10(pooled_stats$p)
      r_vw[, v] <- pooled_stats$resid # ("+", res) / length(res)

      # TMP for debugging
      vw_message(worker_id, ': vertex', v, 'done.')
    }
  }

  out <- list(c_vw, s_vw, t_vw, p_vw, r_vw)
  names(out) <- res_bk_names # c("coefficients", "standard_errors", "t_values", "p_values", "residuals")

  # Post-processing ============================================================

  # Save the stack names (i.e. fixed terms) to a lookup file
  stack_ids <- data.frame('stack_number' = 1:(length(fixed_terms)),
                          'stack_name' = fixed_terms)
  utils::write.table(stack_ids, file.path(outp_dir, "stack_names.txt"),
                     sep = "\t", row.names = FALSE)

  # Split all coefficients into separate files and save them in the final directory
  vw_message("Converting coefficients, SEs, t- and p-values to .mgh format...",
             verbose=verbose)

  results_grid <- expand.grid(stack_ids$stack_number, # how many terms
                    res_bk_names[-length(res_bk_names)]) # all statistics except the residuals

  stats_mgh_paths <- file.path(
    outp_dir,
    paste(hemi, paste0("stack",results_grid[,1]), results_grid[,2], "mgh", sep = ".")
  )

  lapply(seq_len(nrow(results_grid)), function(grid_row) {
    fbm2mgh(fbm = out[[ results_grid[grid_row, 2] ]],
            fname = stats_mgh_paths[grid_row],
            filter = as.integer(results_grid[grid_row, 1]))
    })

  vw_message("Converting residuals to .mgh format...", verbose=verbose)
  resid_mgh_path <- file.path(outp_dir, paste(hemi, "residuals.mgh", sep = "."))
  fbm2mgh(fbm = out[[ res_bk_names[length(res_bk_names)] ]],
          fname = resid_mgh_path)

  # Estimate full-width half maximum (using FreeSurfer)
  vw_message("Estimating data smoothness for multiple testing correction...",
             verbose=verbose)

  fwhm <- estimate_fwhm(outp_dir = outp_dir,
                        hemi = hemi,
                        resid_file = resid_mgh_path,
                        mask = good_verts,
                        fs_template = fs_template)
  if (fwhm > 30) {
    vw_message("Estimated smoothness is ", fwhm, ", which is really high. Reduced to 30.",
               verbose=verbose)
    fwhm <- 30
  } else if (fwhm < 1) {
    vw_message("Estimated smoothness is ", fwhm, ", which is really low. Increased to 1.",
               verbose=verbose)
    fwhm <- 1
  }

  vw_message("Clusterwise correction...", verbose=verbose)

  for ( stack_n in seq_along(fixed_terms) ){
    # message2("\n", verbose = verbose)
    compute_clusters(outp_dir = outp_dir,
                     hemi = hemi,
                     term_number = stack_n,
                     fwhm = fwhm,
                     FS_HOME = FS_HOME,
                     mcz_thr = mcz_thr,
                     cwp_thr = cwp_thr,
                     mask = file.path(outp_dir, "finalMask.mgh"))
    }

  return(out)
}

#' @title
#' Run a single \code{lme4::lmer} model and extract stats
#'
#' @param imp : The phenotype dataset (in verywise format).
#' @param y : A vector of outcome values (i.e. a single vertex from the
#' supersubject matrix).
#' @param formula : R formula describing the linear mixed model (using \code{lme4} notation)
#' @param model_template : (default = NULL) `single_lmm` can leverage an "update"-based
#' workflow to minimize repeated parsing and model construction overhead.
#' @param pvalues : (default : TRUE) whether to include p-values computed using the
#' *t-as-z* approach (see Details).
#'
#' @details
#' No additional parameters are currently passed to the \code{lme4::lmer} call
#' P-values are estimated using the *t-as-z* approach at the moment. This is known
#' to the anti-conservative for small sample sizes. However we preferred a
#' relatively lenient (and computationally inexpensive) solution at this stage.
#' We will be addressing Type I error mores strictly at the cluster forming stage.
#'
#' @return A list with two elements:
#' \enumerate{
#' \item \code{"stats"}: a dataframe with estimates, SEs and p-values for each
#' fixed effect term)
#' \item \code{"resid"}: a vector of residuals for the given model.
# #' \item \code{"df"}: residual degrees of freedom for the give model.
#' }
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
single_lmm <- function(imp, y, formula, model_template = NULL, pvalues = TRUE) {
  # Add (vertex) outcome to (single) dataset
  imp[all.vars(formula)[1]] <- y

  # Fit linear mixed model using `lme4::lmer`
  if (!is.null(model_template)) {
    fit <- stats::update(model_template, data = imp)
  } else {
    fit <- suppressMessages(lmer(formula = formula, data = imp)) }

  # coef(fit)$id # A matrix of effects by random variable
  # lme4::ranef(fit) # extract random effects (these should sum to 0)
  # lme4::VarCorr(fit) # estimated variances, SDs, and correlations between the random-effects terms

  # Extract estimates, standard errors and p-values
  stats <- data.frame(
    "qhat" = as.matrix(fixef(fit)), # Fixed effects estimates
    "se" = as.matrix(summary(fit)$coefficients[, "Std. Error"]) # corresponding standard errors
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
  if (pvalues){
    stats$tval <- stats$qhat/stats$se
    stats$pval <- 2 * (1 - pnorm(abs(stats$tval)))
  }

  # Save row.names (i.e. terms) as a column so these can be grouped later
  stats <- data.frame("term" = row.names(stats), stats, row.names=NULL)

  # Also extract model residuals for smoothness estimation
  resid <- residuals(fit)

  # Finally, extract residual degrees of freedom for Barnard-Rubin adjustment
  # Normally, this would be the number of independent observation minus the
  # number of fitted parameters, but not exactly what is done here.
  # Following the `broom.mixed` package approach, which `mice::pool` relies on
  # df <- df.residual(fit)

 return(list("stats" = stats, "resid" = resid))
}
