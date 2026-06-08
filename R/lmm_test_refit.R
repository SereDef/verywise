precompile_model2 <- function(formula,
                             data_list,
                             tmp_y,
                             measure,
                             weights,
                             REML = TRUE,
                             lmm_control = lme4::lmerControl(calc.derivs=FALSE),
                             verbose) {
  
  if (verbose) cli::cli_progress_step('Construct model template', spinner=TRUE)

  # Fit model template (once per imputed dataset)
  model_template <- lapply(data_list, function(imp) {

    imp[paste0("vw_", measure)] <- tmp_y  # dummy outcome

    if (!is.null(weights) && !is.numeric(weights)) {
      weight_vec <- imp[[weights]]
    } else {
      weight_vec <- weights
    }

    tmp_fit <- suppressMessages(
      lme4::lmer(formula = formula,
                 data = imp,
                 REML = REML,
                 control = lmm_control, 
                 weights = weight_vec)
    )
    # Embed the control settings into the call so update() can access them ?  
    # tmp_fit@call$control <- lmm_control
    # mp_fit@call$REML <- REML

    return(tmp_fit)
  })

  if (verbose) cli::cli_progress_done()
  model_template
}


single_lmm2 <- function(model_template_i, y) {

  error_msg <- NULL
  warning_msg <- character(0)

  # Refit linear mixed model using `lme4::refit`
  fit <- withCallingHandlers(
    tryCatch(
      lme4::refit(model_template_i, newresp = y),
      error = function(e) { 
        error_msg <<- conditionMessage(e)
        NULL
      }),
    warning = function(w) {
      warning_msg <<- c(warning_msg, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      warning_msg <<- c(warning_msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  if (!is.null(error_msg)) {
    result <- list("error" = error_msg)
    return(result)
  }

  coefs <- fixef(fit) # Fixed effects estimates
  ses   <- sqrt(diag(as.matrix(vcov(fit)))) # Their standard errors
  
  fixed_stats <- data.frame(
    term = names(coefs),
    qhat = as.numeric(coefs),
    se   = as.numeric(ses),
    row.names = NULL,
    check.names = FALSE)
  
  # Also extract model residuals for smoothness estimation
  resid <- residuals(fit)

  # Model performance ------------------------------------
  
  is_singular <- as.numeric(isSingular(fit))

  # Add singularity to performance grid, but do not print the warning
  if (is_singular) {
    warning_msg <- warning_msg[!grepl("boundary (singular) fit", warning_msg, 
                                      fixed = TRUE)]
  }

  # Add model performance metrics 
  aic <- AIC(fit)
  # icc <- icc(fit)[['ICC_adjusted']]
  # r2 <- t(r2_nakagawa(fit))
  

  # Use lme4 directly – no extra dependencies, no interactive prompts:
  vc <- VarCorr(fit)
  var_rand <- sum(sapply(vc, function(x) sum(diag(as.matrix(x)))))
  var_resid <- attr(vc, "sc")^2
  var_fix <- var(predict(fit, re.form = NA)) # var(as.vector(X %*% beta))

  icc <- safe_calc(var_rand / (var_rand + var_resid))
  r2_margin <- safe_calc(var_fix  / (var_fix + var_rand + var_resid))
  r2_condit <- safe_calc((var_fix + var_rand) / (var_fix + var_rand + var_resid))
  
  # varpart <- get_variance(fit, tolerance = 1e-12) # more lenient than default

  # icc <- safe_calc(varpart$var.random / (varpart$var.random + varpart$var.residual))

  # r2_margin <- safe_calc(varpart$var.fixed / 
  #   (varpart$var.fixed + varpart$var.random + varpart$var.residual))
  
  # r2_condit <- safe_calc((varpart$var.fixed + varpart$var.random) /
  #       (varpart$var.fixed + varpart$var.random + varpart$var.residual))
  
  # dropping names: singularity, aic, icc, r2_marginal, r2_conditional
  perf <- c(is_singular, aic, icc, r2_margin, r2_condit)
  
  # Extract residual degrees of freedom for Barnard-Rubin adjustment
  # Normally, this would be the number of independent observation minus the
  # number of fitted parameters, but not exactly what is done here.
  # Following the `broom.mixed` package approach, which `mice::pool` relies on
  # resid_df <- df.residual(fit)

  # TODO: implement "stack of interest"?

  # coef(fit)$id # A matrix of effects by random variable
  # lme4::ranef(fit) # extract random effects (these should sum to 0)
  # lme4::VarCorr(fit) # estimated variances, SDs, and correlations between the random-effects terms

  list("stats" = fixed_stats, "resid" = resid, "model_fit" = perf,
       "warning" = warning_msg)

}

run_vw_lmm2 <- function(
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
  FS_HOME = Sys.getenv("FREESURFER_HOME"),
  fwhm = 10,
  mcz_thr = 30,
  cwp_thr = 0.025,
  # Output control
  save_optional_cluster_info = FALSE,
  save_ss = FALSE,
  save_residuals = FALSE,
  verbose = TRUE) {
  
  vw_init_message('Linear mixed model', verbose = verbose)

  # Make output look nice in non-interactive sessions
  # old_cli_opts <- vw_setup_cli_output()
  # if (!is.null(old_cli_opts)) on.exit(options(old_cli_opts), add = TRUE)

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
  check_numeric_param(fwhm, lower = 1, upper = 30)
  check_numeric_param(mcz_thr, set=c(13, 20, 23, 30, 33, 40))
  check_numeric_param(cwp_thr, set=c(0.025, 0.05))

  n_cores <- check_cores(n_cores)

  if (verbose) cli::cli_progress_done()

  # FreeSurfer
  check_freesurfer_setup(FS_HOME, verbose = verbose)

  # Avoid bigstatsr warning about lost precision (float vs. double)
  old_opts <- options(bigstatsr.downcast.warning = FALSE)
  on.exit(options(old_opts), add = TRUE)

  # Esure reproducible seeds in parallel settings
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  
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

  data1 <- data_list[[1]]

  # Prepare chunk sequence =====================================================
  chunk_seq <- make_chunk_sequence(good_verts, chunk_size = chunk_size)

  # Number of vertices
  vw_n <- length(is_cortex); rm(is_cortex)
  # Number of terms (excluding random terms)
  fe_n <- length(fixed_terms)
  # Number of participants*timepoint (long format)
  n_obs <- nrow(data1)
  # Number of (imputed) datasets
  m <- length(data_list)

  if (verbose) cli::cli_progress_done()

  vw_message(c(">" = "Ready to run {.val {length(good_verts)}} models 
                     (split in {.val2 {length(chunk_seq)}} chunks)."),
             verbose = verbose)

  vw_message("Statistical model fitting", type='step', verbose = verbose)

  # Cache the model frame: `single_lmm` uses an "update"-based workflow to minimize
  # repeated parsing and model construction overhead
  model_template <- precompile_model2(formula = formula, data_list = data_list, 
    tmp_y = ss[, good_verts[1]], measure = measure, weights = weights,
    lmm_control = lmm_control, REML = REML, verbose = verbose)
  
  vw_message(c("i"= "model includes {.val2 {fe_n}} fixed parameters and {.val2 {summary(model_template[[1]])$ngrps}} groups"), 
    verbose = verbose)

  # Prepare FBM output =========================================================

  result_path <- file.path(outp_dir, paste(hemi, measure, sep = "."))

  # Temporary output matrices
  res_bk_names <- c("coef", "se", "p", "fitstats", "resid")
  res_bk_paths <- build_output_bks(result_path, res_bk_names = res_bk_names,
                                   verbose = verbose)
  # These files will be removed "on.exit" by convert_to_mgh

  fbm_precision <- "float" # single precision – 32 bits

  c_vw <- bigstatsr::FBM(fe_n, vw_n, init = NA_real_, type = fbm_precision,
                         backingfile = res_bk_paths["coef"])  # Coefficients
  s_vw <- bigstatsr::FBM(fe_n, vw_n, init = 0, type = fbm_precision,
                         backingfile = res_bk_paths["se"])    # Standard errors
  p_vw <- bigstatsr::FBM(fe_n, vw_n, init = 1, type = fbm_precision,
                         backingfile = res_bk_paths["p"])     # P values
  r_vw <- bigstatsr::FBM(n_obs, vw_n, init = 0, type = fbm_precision,
                         backingfile = res_bk_paths["resid"]) # Residuals
  # Fit statistics: singular_fits, aic, icc, r2_marginal, r2_conditional
  f_vw <- bigstatsr::FBM(5, vw_n, init = NA_real_, type = fbm_precision,
                         backingfile = res_bk_paths["fitstats"])   
  
  log_file <- paste0(result_path, ".issues.log") # Log model fitting issues

  # Parallel analyses ==========================================================

  progress_file <- paste0(result_path, ".progress.log")
  on.exit(if (file.exists(progress_file)) file.remove(progress_file), add = TRUE)

  # Progress bar setup # note progressr only works with doFuture not doParallel
  if (verbose) {
    cli::cli_progress_step("Fitting linear mixed models... this may take some time, check the {.file {basename(progress_file)}} file for updates.", 
    spinner=TRUE)
  }
  with_parallel(n_cores = n_cores, 
    seed = seed,
    verbose = verbose, 
    expr = {
      foreach::foreach(chunk = chunk_seq, 
        .packages = c("bigstatsr"), 
        .export = c("single_lmm2", "vw_pool",
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
        out_stats <- lapply(model_template, single_lmm2,
                            y = vertex)

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
        f_vw[, v] <- pooled_stats$fitstats
        r_vw[, v] <- pooled_stats$resid
      }
    }
  })

  if (verbose) cli::cli_progress_done()

  out <- list(c_vw, s_vw, p_vw, f_vw, r_vw)
  # "coefficients", "standard_errors", "p_values", "fit_statistics", "residuals"
  names(out) <- res_bk_names

  # Post-processing ============================================================
  vw_message("Post-processing", type='step', verbose = verbose)

  # Save model statistics into separate .mgh files
  
  convert_to_mgh(out,
                 result_path,
                 stacks = seq_along(fixed_terms),
                 stat_names = c(res_bk_names,'-log10p'),
                 verbose = verbose)

  resid_mgh_path <- paste(result_path, "residuals.mgh", sep = ".")
  if (!save_residuals) on.exit(file.remove(resid_mgh_path), add = TRUE)

  # Estimate full-width half maximum (using FreeSurfer) ========================

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
  
  # clusters
  ct_vw <- NULL 

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
      # create storage once
      ct_bk_path <- build_output_bks(result_path, res_bk_names = c('clust'), verbose = FALSE)
      ct_vw <- bigstatsr::FBM(fe_n, vw_n, init = NA_real_, 
        type = fbm_precision, backingfile = ct_bk_path["clust"])
    }

    ct_vw[stack_n, ] <- ocn

  }

  if (verbose) cli::cli_progress_done()
  
  # Print summary 
  vw_summarize_model_fit(fitstats = out$fitstats, verbose = verbose)
  if (!is.null(ct_vw)) {
    out[['clust']] <- ct_vw
    vw_summarize_model_clusters(coef = out$coef, clust = out$clust, 
      term_names = fixed_terms, verbose = verbose)
  } else {
    # TODO: mask 0 verts (done... NA now )
    vw_summarize_model_est(coef = out$coef, term_names = fixed_terms, verbose = verbose)
  }

  vw_message("Done! :)", type='step', verbose = verbose)

  return(out)
}
