#' @title
#' Run vertex-wise linear mixed model using \code{lme4::lmer()}
#'
#' @description
#' This is is the main function in v0 of verywise runs main analyses.
#'
#' @param formula : model formula object (this should specify a LME model)
#' @param pheno : the phenotype data or a path to the data file.
#' @param subj_dir : path to the FreeSurfer data, this expects a verywise structure.
#' @param hemi : (default = "lh") hemisphere.
#' @param measure : (default = \code{gsub("vw_", "", all.vars(formula)[1])}) vertex-wise measure
#' this should be the same as specified in formula
#' @param fwhm : (default = 10) full-width half maximum value
#' @param target : (default = "fsaverage") template on which to register vertex-wise data.
#' @param n_cores : (default = 1) number of cores for parallel processing.
#' @param model : (default = \code{"lme4::lmer"}) # "stats::lm"
#'
#' @import bigstatsr
#' @import foreach
#' @import lme4
#' @import stats
#' @importFrom doParallel registerDoParallel
#' @importFrom car Anova
#'
#' @return A list of file-backed matrices containig pooled coefficients, SEs,
#' t- and p- values and residuals.
#'
#' @author Serena Defina, 2024.
#'
run_vw_lmm <- function(formula, # model formula
                       pheno = NULL,
                       subj_dir,
                       hemi = "lh",
                       measure = gsub("vw_", "", all.vars(formula)[1]),
                       fwhm = 10,
                       target = "fsaverage",
                       n_cores = 1,
                       model = "lme4::lmer" # "stats::lm"
) {
  # Read phenotype data (if not already loaded) ================================
  if (is.null(pheno) & !exists("pheno", mode="list", envir=globalenv())) {
    pheno <- utils::read.csv(file.path(subj_dir, "phenotype.csv"))
  } else if (is.character(pheno) & file.exists(pheno)) {
    pheno <- utils::read.csv(pheno)
  }

  # Transform to list of dataframes (imputed and single datasets alike)
  data_list <- imp2list(pheno)

  # Read and clean vertex data =================================================
  ss_file_name <- file.path(subj_dir, paste0(
    hemi, ".", measure,
    "supersubject.rds"
  ))
  if (file.exists(ss_file_name)) {
    message("Reading super-subject file from: ", ss_file_name)

    ss <- bigstatsr::big_attach(ss_file_name)
  } else {
    message("Building super-subject...")

    ss <- build_supersubject(subj_dir,
      folder_id = pheno$folder_id,
      files_list = list.dirs.till(subj_dir, n = 2),
      measure = measure,
      hemi = hemi,
      fwhmc = paste0("fwhm", fwhm),
      target = target,
      n_cores = n_cores,
      mask = TRUE,
      save_rds = TRUE
    )
  }

  # Additionally check that there are no vertices that contain any 0s. These may be
  # located at the edge of the cortical map and are potentially problematic
  problem_verts <- fbm_col_has_0(ss)
  if (sum(problem_verts) > 0) {
    message(
      "Removing ", sum(problem_verts),
      " vertices that contained 0 values."
    )
  }
  good_verts <- which(!problem_verts)

  # Prepare chunk sequence =====================================================
  chunk_seq <- make_chunk_sequence(good_verts)

  # Unpack model ===============================================================
  model_info <- get_terms(formula, data_list)

  fixed_terms <- model_info[[1]]

  # Number of vertices
  vw_n <- length(good_verts)
  # Number of terms (excluding the random!?)
  fe_n <- length(fixed_terms)
  # Number of participants*timepoint (long format)
  n_obs <- nrow(model_info[[2]]) # nrow(data_list[[1]])
  # Number of (imputed) datasets
  m <- length(data_list)

  # Prepare FBM output =========================================================
  res_bk_names <- c("coef", "se", "t", "p", "resid")
  res_bk_paths <- file.path(
    subj_dir,
    paste(hemi, measure, res_bk_names, sep = ".")
  )

  # TMP: remove files if they exist
  for (bk in res_bk_paths) {
    if (file.exists(paste0(bk,".bk"))) file.remove(paste0(bk,".bk"))
  }

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

  # Parallel loop per chunk  ===================================================
  message("Running analyses...\n")

  future::plan("multisession", workers = n_cores) # Should let the user do it instead..?
  # if (n_cores >= parallelly::availableCores(omit = 1)) { warning("Uff, slow down kid")}
  # future::plan(
  #   list(
  #     future::tweak(
  #       future::multisession,
  #       workers = 2),
  #     future::tweak(
  #       future::multisession,
  #       workers = 4)
  #   )
  # ) # 8 cores in total

  # cl <- parallel::makeForkCluster(n_cores, outfile = "")
  # doParallel::registerDoParallel(cl)
  # on.exit(parallel::stopCluster(cl))

  # utils::capture.output(
  # pb <- utils::txtProgressBar(0, nrow(chunk_seq), style = 3)
  # , file = "/dev/null")
  p <- progressr::progressor(along = 1:nrow(chunk_seq))

  # lapply(1:nrow(chunk_seq), function(chunk) {
  # parallel::parLapply(cl, 1:nrow(chunk_seq), function(chunk) {
  # chunk <- NULL
  # foreach::foreach(chunk = 1:nrow(chunk_seq), .combine = "c") %dopar% {
  # future.apply::future_lapply(1:nrow(chunk_seq), function(chunk){
  furrr::future_walk(1:nrow(chunk_seq), function(chunk) {
    # Progress bar
    p()
    # utils::setTxtProgressBar(pb, chunk)

    id <- good_verts[chunk_seq[chunk, 1]:chunk_seq[chunk, 2]]

    Y <- ss[, id]

    lapply(1:ncol(Y), function(vertex) { # TODO: fix nested parallelisation here
    # parallel::parLapply(cl, 1:ncol(Y), function(vertex) {
      # Fetch brain data
      y <- Y[, vertex]

      # Empty list for storing results
      qhat <- se <- pval <- resid <- vector(mode = "list", length = m)

      # Loop through imputed datasets
      lapply(seq_along(data_list), function(imp) {
        dset <- data_list[[imp]]
        dset[all.vars(formula)[1]] <- y

        fit <- suppressMessages(lme4::lmer(formula = formula, data = dset))
        # coef(fit)$id is a matrix of effects by random variable?
        qhat[[imp]] <<- as.matrix(lme4::fixef(fit))
        # lme4::ranef(fit) to extract random effects (these should sum to 0)
        # lme4::VarCorr(fit)  # estimated variances, SDs, and correlations between the random-effects terms
        se[[imp]] <<- as.matrix(sqrt(diag(as.matrix(stats::vcov(fit)))))
        resid[[imp]] <<- stats::residuals(fit)
        # TODO: implement stack of interest?
        pval[[imp]] <<- as.matrix(car::Anova(fit, type = 3)[, "Pr(>Chisq)"]) # -log10()
      })

      out_stats <- vw_pool(qhat = qhat, se = se, m = m, fe_n = fe_n)
      # Average residuals across imputed datasets. TODO: pool this instead?
      out_stats$resid <- as.matrix(colMeans(do.call(rbind, resid)))

      # Write results to their respective FBM
      c_vw[, vertex] <<- out_stats$coef
      s_vw[, vertex] <<- out_stats$se
      t_vw[, vertex] <<- out_stats$t
      p_vw[, vertex] <<- -1 * log10(out_stats$p)
      r_vw[, vertex] <<- out_stats$resid # ("+", res) / length(res)

      NULL
    })
    # NULL
  })
  # on.exit(invisible(NULL))

  out <- list(c_vw, s_vw, t_vw, p_vw, r_vw)
  names(out) <- c("coefficients", "standard_errors", "t_values", "p_values", "residuals")

  message("All done!")
  return(out)
}
