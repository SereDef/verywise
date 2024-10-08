#' @title
#' Run vertex-wise linear mixed model using \code{lme4::lmer()}
#'
#' @description
#' This is is the main function in v0 of \code{verywise}.
#' It runs the \code{\link{single_lmm}} function across all vertices.
#'
#' @param formula : model formula object (this should specify a LME model)
#' @param pheno : the phenotype data or a path to the data file.
#' @param subj_dir : path to the FreeSurfer data, this expects a verywise structure.
#' @param hemi : (default = "both") which hemispheres to run.
#' @param seed : (default = 3108) random seed.
#' @inheritDotParams hemi_vw_lmm measure fwhm target model
#'
#' @return A list of file-backed matrices containing pooled coefficients, SEs,
#' t- and p- values and residuals.
#'
#' @author Serena Defina, 2024.
#'
#' @export
#'
run_vw_lmm <- function(formula, # model formula
                       subj_dir,
                       pheno = NULL,
                       hemi = c("both", "lh","rh"),
                       seed = 3108,
                       ...
) {
  # Read phenotype data (if not already loaded) ================================
  if (is.null(pheno)) {
    if (exists("pheno", mode="list", envir=globalenv())) {
      pheno <- get("pheno", envir = globalenv())
    } else if (file.exists(file.path(subj_dir, "phenotype.csv"))){
      message("Reading phenotype file from subject directory folder.")
      pheno <- utils::read.csv(file.path(subj_dir, "phenotype.csv"))
    } else {
      stop("Provide phenotype data pls.")
    }
  } else if (is.character(pheno) & file.exists(pheno)) {
    pheno <- utils::read.csv(pheno)
  }

  # Transform to list of dataframes (imputed and single datasets alike)
  data_list <- imp2list(pheno)

  hemi <- match.arg(hemi)

  set.seed(seed)

  if (hemi == "both") {
    out <- furrr::future_map(c("lh","rh"), function(h) {

      hemi_vw_lmm(formula = formula,
                  subj_dir = subj_dir,
                  data_list = data_list,
                  hemi = h,
                  seed = seed,
                  ...)},
      .options = furrr::furrr_options(seed = TRUE))

  } else {
    out <- hemi_vw_lmm(formula = formula,
                subj_dir = subj_dir,
                data_list = data_list,
                hemi = hemi,
                seed = seed,
                ...)
  }
  message("All done!")
  out
}

#' @title
#' Run vertex-wise linear mixed model in one hemisphere using \code{lme4::lmer()}
#'
#' @description
#' It runs the \code{\link{single_lmm}} function across vertices in a single hemisphere.
#'
#' @param formula : model formula object (this should specify a LME model)
#' @param data_list : the phenotype as a list of dataframes (generated by \code{\link{imp2list}})
#' @param subj_dir : path to the FreeSurfer data, this expects a verywise structure.
#' @param hemi : (default = "both") hemispheres to run.
#' @param measure : (default = \code{gsub("vw_", "", all.vars(formula)[1])}) vertex-wise measure
#' this should be the same as specified in formula
#' @param fwhm : (default = 10) full-width half maximum value
#' @param target : (default = "fsaverage") template on which to register vertex-wise data.
#' @param seed : (default = 3108) random seed.
#' @param apply_cortical_mask : (default = TRUE) remove vertices that are not on the cortex.
#' @param model : (default = \code{"lme4::lmer"}) # "stats::lm"
#'
#' @return A list of file-backed matrices containing pooled coefficients, SEs,
#' t- and p- values and residuals.
#'
#' @author Serena Defina, 2024.
#'
#' @export
#'
hemi_vw_lmm <- function(formula, # model formula
                        subj_dir,
                        data_list,
                        hemi,
                        measure = gsub("vw_", "", all.vars(formula)[1]),
                        fwhm = 10,
                        target = "fsaverage",
                        seed = 3108,
                        apply_cortical_mask = TRUE,
                        model = "lme4::lmer" # "stats::lm"
) {

  # Read and clean vertex data =================================================
  ss_file_name <- file.path(subj_dir, paste0(
    hemi, ".", measure,
    ".supersubject.rds"
  ))

  if (file.exists(ss_file_name)) {
    message("Reading super-subject file from: ", ss_file_name)

    ss <- bigstatsr::big_attach(ss_file_name)

  } else {

    ss <- build_supersubject(subj_dir,
      folder_id = data_list[[1]]$folder_id,
      files_list = list.dirs.till(subj_dir, n = 2),
      measure = measure,
      hemi = hemi,
      fwhmc = paste0("fwhm", fwhm),
      target = target,
      mask = apply_cortical_mask,
      save_rds = FALSE
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

  # if (requireNamespace("progressr", quietly = TRUE)) {
  #   # progressr::handlers(global = TRUE)
  #   p <- progressr::progressor(along = good_verts) #1:nrow(chunk_seq)
  # }

  # Don't do chucking to avoid nested parallel processes that would require
  # more involved future::plan spec

  set.seed(seed)

  furrr::future_walk(good_verts, function(vertex) {

    # if (requireNamespace("progressr", quietly = TRUE)) { p() }
    # Fetch brain data

    y <- ss[, vertex]

    # Loop through imputed datasets and run analyses
    out_stats <- lapply(data_list, single_lmm, y = y, formula = formula,
                        pvalues = (m == 1))

    # Pool results
    pooled_stats <- vw_pool(out_stats, m = m)

    # Write results to their respective FBM
    c_vw[, vertex] <<- pooled_stats$coef
    s_vw[, vertex] <<- pooled_stats$se
    t_vw[, vertex] <<- pooled_stats$t
    p_vw[, vertex] <<- -1 * log10(pooled_stats$p)
    r_vw[, vertex] <<- pooled_stats$resid # ("+", res) / length(res)

  },
  .options = furrr::furrr_options(seed = TRUE),
  .progress = TRUE)

  out <- list(c_vw, s_vw, t_vw, p_vw, r_vw)
  names(out) <- c("coefficients", "standard_errors", "t_values", "p_values", "residuals")

  return(out)
}

#' @title
#' Run a single \code{lme4::lmer} model and extract stats
#'
#' @param imp : The phenotype dataset (in verywise format).
#' @param y : A vector of outcome values (i.e. a single vertex from the
#' supersubject matrix).
#' @param formula : R formula describing the linear mixed model (using \code{lme4} notation)
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
single_lmm <- function(imp, y, formula, pvalues = TRUE) {
  # Add (vertex) outcome to (single) dataset
  imp[all.vars(formula)[1]] <- y

  # Fit linear mixed model using `lme4::lmer`
  fit <- suppressMessages(lmer(formula = formula, data = imp))

  # coef(fit)$id # A matrix of effects by random variable
  # lme4::ranef(fit) # extract random effects (these should sum to 0)
  # lme4::VarCorr(fit) # estimated variances, SDs, and correlations between the random-effects terms

  # Extract estimates, standard errors and p-values
  stats <- data.frame(
    "qhat" = as.matrix(fixef(fit)), # Fixed effects estimates
    "se" = as.matrix(sqrt(diag(as.matrix(vcov(fit))))) # corresponding standard errors
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

 return(list("stats" = stats, "resid" = resid)) # , df
}
