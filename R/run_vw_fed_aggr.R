#' @title Run federated aggregation for vertex-wise distributed LMM
#' 
#' @description
#' Performs the aggregation step of a privacy-preserving federated linear mixed
#' model (LMM) on vertex-wise neuroimaging data. This function never sees raw
#' participant data — it operates entirely on pre-computed sufficient statistics
#' (site payloads) uploaded by each site. It fits a random-intercept LMM with a
#' single shared variance component ratio \eqn{\lambda = \sigma_b^2 / \sigma_\varepsilon^2},
#' selected via a 1-D profile likelihood grid search (rather than iterative
#' site communication).
#'
#' @details
#' **Model structure.** The covariance matrix is assumed to follow
#' \deqn{V = \sigma_\varepsilon^2 (I + \lambda Z Z^\top)}
#' where \eqn{Z} is a block-diagonal indicator of site membership. Under this
#' structure, the Sherman–Morrison–Woodbury identity reduces the \eqn{N \times N}
#' inverse to site-level rank-1 updates, making the algorithm \eqn{O(K p^2)}
#' rather than \eqn{O(N^3)}.
#'
#' **Grid search.** For each candidate \eqn{\lambda}, the function reconstructs
#' the global Hessian \eqn{A(\lambda) = X^\top V^{-1} X} and its Cholesky factor,
#' then computes the profile (RE)ML log-likelihood per vertex in chunks. The
#' \eqn{\lambda} that maximises the log-likelihood is selected independently for
#' each vertex.
#'
#' **Chunked processing.** Vertices are processed in chunks of size
#' \code{chunk_size} to keep memory usage bounded regardless of surface
#' resolution. Results are written directly to memory-mapped \pkg{bigstatsr}
#' FBM files on disk.
#'
#' @param site_names Character vector of site identifiers. Each site must have a
#'   corresponding \code{<site>.tar.gz} archive in \code{inpt_dir}.
#' @param formula A two-sided \code{\link[stats]{formula}} specifying the fixed
#'   effects. The LHS names the surface measure (e.g. \code{thickness ~ age + sex}).
#' @param inpt_dir Character. Directory containing the site payload archives.
#' @param outp_dir Character or \code{NULL}. Directory for output FBM backing
#'   files. Defaults to \code{inpt_dir} if \code{NULL}.
#' @param hemi Character. Hemisphere: \code{"lh"} (default) or \code{"rh"}.
#' @param fs_template Character. FreeSurfer template surface used for the cortical
#'   mask (e.g. \code{"fsaverage"}, \code{"fsaverage5"}). Default: \code{"fsaverage"}.
#' @param fwhm Numeric. Full-width at half-maximum (mm) of the smoothing kernel
#'   applied at the site level. Used for informational messages only; smoothing
#'   must be applied before payload generation. Default: \code{10}.
#' @param seed Integer. Random seed passed to the parallel backend for
#'   reproducibility. Default: \code{3108}.
#' @param n_cores Integer. Number of parallel workers. \code{1} (default) runs
#'   sequentially; values \code{> 1} use \pkg{doParallel}. 
#' @param chunk_size Integer. Number of vertices processed per parallel job.
#'   Larger values are faster but use more memory. Default: \code{1000}.
#' @param REML Logical. If \code{TRUE} (default), use restricted maximum
#'   likelihood (REML); otherwise use ML.
#' @param lambda_grid Numeric vector of candidate \eqn{\lambda} values for the
#'   grid search, or \code{NULL} to use the default log-spaced grid spanning
#'   \eqn{[10^{-6}, 10^{2}]}.
#' @param chol_tol Numeric. Minimum diagonal pivot in the Cholesky factor below
#'   which a \eqn{\lambda} value is rejected as numerically singular.
#'   Default: \code{1e-12}.
#' @param ridge Numeric. Ridge penalty added to the diagonal of \eqn{A(\lambda)}
#'   before factorisation. Use a small positive value (e.g. \code{1e-6}) when
#'   predictors are near-collinear. Default: \code{0} (no regularisation).
#' @param verbose Logical. Whether to print progress messages. Default: \code{TRUE}.
#'
#' @return A named list of four \code{\link[bigstatsr]{FBM}} objects, each of
#'   dimension described below. Rows correspond to model terms; columns to
#'   vertices.
#' \describe{
#'   \item{\code{coef}}{\eqn{p \times V} matrix of fixed-effect estimates \eqn{\hat\beta}.}
#'   \item{\code{se}}{\eqn{p \times V} matrix of standard errors.}
#'   \item{\code{pval}}{\eqn{p \times V} matrix of two-sided p-values (t-distribution, \code{df} degrees of freedom).}
#'   \item{\code{variances}}{\eqn{2 \times V} matrix; row 1 = \eqn{\hat\sigma^2_\varepsilon}, row 2 = \eqn{\hat\tau^2 = \hat\lambda \hat\sigma^2_\varepsilon}.}
#' }
#'
#' @seealso \code{\link{chunk_dlmm}}, \code{\link{decompress_site_payloads}},
#'   \code{\link{vw_logLL}}
#' 
#' @author Serena Defina, 2026.
#' 
#' @export
#' 
run_vw_fed_aggr <- function(
  site_names,  
  formula,
  inpt_dir,
  outp_dir = NULL,
  # Brain data processing
  hemi = c("lh", "rh"), 
  fs_template = "fsaverage",
  fwhm = 10,
  # Reproducibility and parallel processing
  seed = 3108,
  n_cores = 1,
  chunk_size = 1000,
  # Modeling 
  REML = TRUE,
  lambda_grid = NULL,
  # refine_steps = 2L,
  chol_tol = 1e-12,
  ridge = 0,
  verbose = TRUE) {
  
  # Check user input ===========================================================

  hemi <- match.arg(hemi)
  measure <- check_formula(formula)

  # TODO enforce site_name; define term names to ensure harmonisarion

  vw_pretty_message(outcome_name(hemi, measure), verbose = verbose)
  vw_message('* Model: ', deparse(formula), verbose = verbose)

  site_payloads <- decompress_site_payloads(site_names, inpt_dir, hemi, measure)

  K <- length(site_names) # Number of sites
  N <- sum(site_payloads$n_obs) # Total N
  p <- length(site_payloads$terms) # Number of terms

  df <- if (REML) (N - p) else N

  # if (df <= 0) .vw_stop("Not enough df: N - p <= 0.")

  vw_message('* Combined N = ', N, ' (across ', K, ' sites).', verbose = verbose)

  # ==== Let's get started ======================================================
  
  # Auto grid: wide but not absurd; assumes lambda mostly in [1e-6, 1e2]
  lambda_grid <- check_lambda_grid(lambda_grid)

  # Profile lambda (static): w(λ), A(λ),
  lambda_prof <- Filter(Negate(is.null), lapply(lambda_grid, function(lambda) {

    # Site weight: how much "shrinkage" or correlation is attributed
    # to the random effect for each site (under this ratio)
    site_weights <- lambda / (1 + site_payloads$n_obs * lambda)  # dim: K

    # Reconstruct the global Hessian matrix A(λ) from site-specific X matrices
    Amat <- Reduce("+", Map(\(XtX, X1, wk) XtX - wk * tcrossprod(X1),
                   site_payloads$XtX, site_payloads$X1, as.list(site_weights))) # dim: p x p

    # (Optional) regularization
    # Add a ridge penalty to the diagonal for numerical stability or regularization, 
    # common in high-dimensional settings or when A is near-singular.
    if (ridge > 0) Amat <- Amat + diag(ridge, p)
    
    # Cholesky factorization: decompose A into RtR
    # This factorized form is used later to efficiently solve for β_hat(λ)
    # via forward/backward substitution.
    Rmat <- vw_chol(Amat, tol = chol_tol) # dim: p x p (upper triangle)

    # If this A(λ) is not positive-definite (or is nearly singular),
    # skip this λ value. I only consider variance ratios where the resulting model 
    # is identifiable and stable.
    if (is.null(Rmat)) return(NULL)
    
    # Compute log-determinant of A(λ), given its Cholesky factor R
    # Note: I use the sum of the natural log of each diagonal element insted of 
    # det(A) to avoid extreme numbers
    logdetA <- 2 * sum(log(diag(Rmat))) # dim: 1
    
    # Compute log-determinant of the block-diagonal covariance matrix V(λ)
    # (ignoring constant terms). In a LMM with random intercepts, the 
    # covariance matrix for site i is Vi = σϵ^2 (Ini + λJni) where 
    # Jni is ni × ni matrix of 1...
    logdetV <- sum(log1p(site_payloads$n_obs * lambda)) # dim: 1

    # logdetA and logdetV are used to compute the REML or ML profile likelihood 
    # value for the current λ (so one can pick the best λ).

    list(lambda = lambda, w = site_weights, R = Rmat, logdetA = logdetA, logdetV = logdetV)
  }))

  # TODO: if all grid points fail (e.g., due to highly collinear predictors or insufficient data), 
  # return an edge-case result?
  if (length(lambda_prof) == 0) stop('Aggregate deisgn matrix is not positive-definite under any lambda assumption.')
  
  # Cortical mask
  is_cortex <- mask_cortex(hemi = hemi, fs_template = fs_template)

  # Additionally check that there are no vertices that have mean 0 across all sites
  problem_verts <- Reduce("&", lapply(site_payloads$psumsY, fbm_col_has_0)) # or "|" for any site?

  good_verts <- which(!problem_verts & is_cortex); rm(problem_verts)

  # Number of vertices
  vw_n <- length(is_cortex); rm(is_cortex)
  
  # Prepare FBM output =========================================================

  result_path <- file.path(outp_dir, paste(hemi, measure, sep = "."))

  # Temporary output matrices
  res_bk_names <- c("coef", "se", "p", 'var') #, "fitstats", "resid") # "t",
  res_bk_paths <- build_output_bks(result_path, res_bk_names = res_bk_names,
                                   verbose = verbose)

  fbm_precision <- "float" # single precision – 32 bits

  c_vw <- bigstatsr::FBM(p, vw_n, init = 0, type = fbm_precision,
                         backingfile = res_bk_paths["coef"])  # Coefficients
  s_vw <- bigstatsr::FBM(p, vw_n, init = 0, type = fbm_precision,
                         backingfile = res_bk_paths["se"])    # Standard errors
  p_vw <- bigstatsr::FBM(p, vw_n, init = 1, type = fbm_precision,
                         backingfile = res_bk_paths["p"])     # P values
  v_vw <- bigstatsr::FBM(2, vw_n, init = 0, type = fbm_precision,
                         backingfile = res_bk_paths["var"]) # variances
  
  
  # log_file <- paste0(result_path, ".issues.log") # Log model fitting issues

  # Prepare chunk sequence =====================================================
  vw_message(" * chunk dataset", verbose = verbose)
  chunk_seq <- make_chunk_sequence(good_verts, chunk_size = chunk_size)

  # Parallel analyses ==========================================================
  vw_message("Running analyses...\n",
             " * dimentions: ", N, " observations x ", length(good_verts),
             " (of ", vw_n, " total) vertices.", verbose = verbose)

  progress_file <- paste0(result_path, ".progress.log")
  on.exit(if (file.exists(progress_file)) file.remove(progress_file), add = TRUE)

  with_parallel(n_cores = n_cores, 
    progress_file = progress_file,
    seed = seed,
    verbose = verbose, 
    expr = {
      foreach::foreach(chunk = chunk_seq, 
        .packages = c("bigstatsr"), 
        .export = c("chunk_dlmm")
                   # "init_progress_tracker", "update_progress_tracker")
    ) %dopar% { # Only parallel if n_cores > 1

      # Progress updates
      # progress_tracker <- init_progress_tracker(chunk, chunk_seq, verbose=TRUE)

      actual_chunk_size <- length(chunk)

      XtY_chunk <- lapply(site_payloads$XtY, function(mat) mat[, chunk]) # K: p x chunk
      sumY_chunk <- lapply(site_payloads$psumsY, function(mat) mat[1, chunk]) # 1tY
      YtY_chunk <- lapply(site_payloads$psumsY, function(mat) mat[2, chunk]) # sum sq Y

      chunk_stats <- chunk_dlmm(lambda_prof = lambda_prof, 
                              X1 = site_payloads$X1,
                              XtY = XtY_chunk,
                              sumY = sumY_chunk,
                              YtY = YtY_chunk,
                              n_terms = p, 
                              chunk_size = actual_chunk_size,
                              n_sites = K, 
                              REML = REML, df = df)
      
      # # Write results to their respective FBM
      c_vw[, chunk] <- chunk_stats$coef
      s_vw[, chunk] <- chunk_stats$se
      p_vw[, chunk] <- chunk_stats$pval
      v_vw[, chunk] <- chunk_stats$bw_var

    } 
  })
  
  list(coef = c_vw, se = s_vw, pval = p_vw, variances = v_vw)
}