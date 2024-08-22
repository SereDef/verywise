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
#' @param seed : (default = 3108) random seed.
#' @param model : (default = \code{"lme4::lmer"}) # "stats::lm"
#'
#'
#' @return A list of file-backed matrices containig pooled coefficients, SEs,
#' t- and p- values and residuals.
#'
#' @author Serena Defina, 2024.
#' @export
#'
run_vw_lmm <- function(formula, # model formula
                       subj_dir,
                       pheno = NULL,
                       hemi = "lh",
                       measure = gsub("vw_", "", all.vars(formula)[1]),
                       fwhm = 10,
                       target = "fsaverage",
                       n_cores = 1,
                       seed = 3108,
                       model = "lme4::lmer" # "stats::lm"
) {
  # Read phenotype data (if not already loaded) ================================
  if (is.null(pheno) & !exists("pheno", mode="list", envir=globalenv())) {
    pheno <- utils::read.csv(file.path(subj_dir, "phenotype.csv"))
  } else if (is.character(pheno) & file.exists(pheno)) {
    pheno <- utils::read.csv(pheno)
  } else  if (exists("pheno", mode="list", envir=globalenv())) {
    pheno <- get(pheno, envir = globalenv())
  } else {
    stop("Provide phenotype data pls.")
  }

  # Transform to list of dataframes (imputed and single datasets alike)
  data_list <- list(pheno, pheno) # imp2list(pheno)

  # Read and clean vertex data =================================================
  ss_file_name <- file.path(subj_dir, paste0(
    hemi, ".", measure,
    "supersubject.rds"
  ))
  if (file.exists(ss_file_name)) {
    message("Reading super-subject file from: ", ss_file_name)

    ss <- bigstatsr::big_attach(ss_file_name)
    # filter it

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

  # Unpack model ===============================================================
  model_info <- get_terms(formula, data_list)

  fixed_terms <- model_info[[1]]

  # Number of vertices
  vw_n <- length(good_verts)
  # Number of terms (excluding the random!?)
  fe_n <- length(fixed_terms)
  # Number of participants*timepoint (long format)
  n_obs <- nrow(model_info[[2]]) # nrow(data_list[[1]])
  if (n_obs != nrow(ss)) {
    stop("qualquadra non cosa")
  }
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

  # Prepare chunk sequence =====================================================
  # chunk_seq <- make_chunk_sequence(good_verts)

  # Parallel analyses ==========================================================
  message("Running analyses...\n")

  # Check number of cores
  if (n_cores > parallelly::availableCores()) {
    warning("You only have access to ", parallelly::availableCores(), " cores. ",
            parallelly::availableCores(omit = 1), " will be used.")
    n_cores <- parallelly::availableCores(omit = 1)
  } else if (n_cores >= parallelly::availableCores(omit = 1)) {
    warning("You are using ", n_cores, " cores, but you only have ",
            parallelly::availableCores(), " in total, other processes may get slower.")
  }

  future::plan("multisession", workers = n_cores) # Should let the user do it instead..?
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

  # utils::capture.output(
  # pb <- utils::txtProgressBar(0, nrow(chunk_seq), style = 3)
  # , file = "/dev/null")
  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::handlers(global = TRUE)
    p <- progressr::progressor(along = good_verts) #1:nrow(chunk_seq)
  }

  # Don't do chucking to avoid nested parallel processes that would require
  # more involved future::plan spec or force an even number
  # furrr::future_walk(1:nrow(chunk_seq), function(chunk) {
  #   # Progress bar
  #   if (requireNamespace("progressr", quietly = TRUE)) p()
  #   # utils::setTxtProgressBar(pb, chunk)
  #
  #   id <- good_verts[chunk_seq[chunk, 1]:chunk_seq[chunk, 2]]
  #
  #   Y <- ss[, id]
  set.seed(seed)

  furrr::future_walk(good_verts, function(vertex) {
      p()
      # Fetch brain data
      y <- ss[, vertex] # Y

      # Loop through imputed datasets and run analyses
      out_stats <- lapply(data_list, single_lmm, y = y, formula = formula)

      # Pool results
      to_pool <- do.call(rbind, lapply(out_stats, `[[`, 1))
      pooled_stats <- vw_pool(qhat = to_pool[,"qhat"], se = to_pool[,"se"],
                              m = m, fe_n = fe_n)
      # Average residuals across imputed datasets.
      # TODO: pool this instead?
      pooled_stats$resid <- as.matrix(colMeans(do.call(rbind,
                                                       lapply(out_stats, `[[`, 2))))

      # Write results to their respective FBM
      c_vw[, vertex] <<- pooled_stats$coef
      s_vw[, vertex] <<- pooled_stats$se
      t_vw[, vertex] <<- pooled_stats$t
      p_vw[, vertex] <<- -1 * log10(pooled_stats$p)
      message("now residuals")
      r_vw[, vertex] <<- pooled_stats$resid # ("+", res) / length(res)
  },
  .options = furrr::furrr_options(seed = TRUE),
  .progress = TRUE)

  out <- list(c_vw, s_vw, t_vw, p_vw, r_vw)
  names(out) <- c("coefficients", "standard_errors", "t_values", "p_values", "residuals")

  message("All done!")
  return(out)
}

#' Run a single \code{lme4::lmer} model and extract stats
#'
#' @param imp : dataset (containing the phenotype data)
#' @param y : vertex outcome (from supersubject)
#' @param formula : formula describing the model
#'
#' @details
#' No additional parameters are currently passed to the \code{lme4::lmer} call
#' P-values are estimated using \code{car::Anova(., type = 3)} at the moment
#'
#' @return A list with two elements: "stats" (a dataframe with estimate, SE and p-value)
#' and "resid" a vector of redisuals for the given model
#'
#' @importFrom lme4 lmer
#' @importFrom lme4 fixef
#' @importFrom stats vcov
#' @importFrom stats residuals
#' @importFrom car Anova
#'
#' @export
#'
single_lmm <- function(imp, y, formula) {
  # Add (vertex) outcome to (single) dataset
  imp[all.vars(formula)[1]] <- y

  # Fit linear mixed model using `lme4::lmer`
  fit <- suppressMessages(lmer(formula = formula, data = imp))

  # coef(fit)$id # A matrix of effects by random variable
  # lme4::ranef(fit) # extract random effects (these should sum to 0)
  # lme4::VarCorr(fit) # estimated variances, SDs, and correlations between the random-effects terms

  stats <- data.frame(
    "qhat" = as.matrix(fixef(fit)), # Fixed effects estimates
    "se" = as.matrix(sqrt(diag(as.matrix(vcov(fit))))), # corresponding standard errors
    "pval" = as.matrix(Anova(fit, type = 3)[, "Pr(>Chisq)"]) # Wald chi-square tests
    # TODO: check pbkrtest:: Kenward-Roger and Satterthwaite Based estimation
    # TODO: ask Bing: why type 3 and not 2
    # -log10() ?
    # TODO: implement stack of interest?
  )

  resid <- residuals(fit)

  return(list(stats, resid))
}
