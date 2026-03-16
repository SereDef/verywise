
# ══ LOCAL SITE ════════════════════════════════════════════════════════════════

#' @title Compute the sum and sum of squares of Y
#' 
#' @description
#' Calculates the cross-product of the design matrix X with the vertex chunk Y (XtY),
#' along with the sum (1'Y) and sum of squares (Y'Y) for each vertex in the chunk.
#'
#' @param X Numeric matrix. The design matrix of dimension \eqn{n \times p}.
#' @param ss_chunk Numeric matrix. A chunk of the ss matrix of dimension \eqn{n \times B},
#'   where B is the block size (number of vertices in the chunk).
#'
#' @return A list containing:
#' \describe{
#'   \item{XtY}{A numeric matrix of dimension \eqn{p \times B}, representing \eqn{X^T Y}.}
#'   \item{psumsY}{A numeric matrix of dimension \eqn{2 \times B}. The first row contains
#'     the column sums of Y (\eqn{1^T Y}), and the second row contains the column sums
#'     of squared Y (\eqn{Y^T Y}).}
#' }
#' 
#' @author Serena Defina, 2026.
#' 
chunk_Ymats <- function(X, ss_chunk) {

  stopifnot(nrow(X) == nrow(ss_chunk))

  # .vw_assert_no_na(Y_block, "Y_block")
  # sum of Y
  Ysum = colSums(ss_chunk) # dim: B (1tY)
  # sum of squared Y
  Ysum2 = colSums(ss_chunk * ss_chunk) # dim: B (YtY)

  list(XtY = crossprod(X, ss_chunk), # dim: p x B
       psumsY = rbind(Ysum, Ysum2))
}

#' @title Compress a site's results into a tar.gz archive
#'
#' @description
#' Creates a gzip-compressed tar archive containing a sufficient statistics
#' from the local site. The archive is written to the output directory and 
#' can be safely shared among insititutions
#'
#' @param outp_dir Character scalar. Directory containing the local output files
#'   to be packaged (e.g., FBM `.rds` + `.bk` pairs and other `.rds` objects).
#' @param site_name Character scalar. Site identifier used as the archive basename
#'   (output file will be `{site_name}.tar.gz`).
#' @param then_rm_files Logical. Remove input files after compressing them? 
#'   Default: TRUE.
#'
#' @return Invisibly returns `NULL`. Called for its side effect of creating the
#'   `{site_name}.tar.gz` file.
#'
#' @examples
#' \dontrun{
#' # Suppose outp_dir contains:
#' #   XtY.rds, XtY.bk, psumsY.rds, psums.bk, static.rds
#' compress_local(outp_dir = "path/to/site_outputs", site_name = "siteA")
#' # Creates "siteA.tar.gz" in the current working directory.
#' }
#'
#' @author Serena Defina, 2026.
#' 
#' @export
#' 
compress_local <- function(outp_dir, site_name, then_rm_files = TRUE) {

  files_to_tar <- list.files(
    path = outp_dir,
    pattern = "(XtY.rds|XtY.bk|psumsY.rds|psumsY.bk|static.rds)$",
    full.names = FALSE)
  
  # Define the output tar path (absolute path so it's valid after setwd)
  tar_path <- file.path(normalizePath(outp_dir), paste0(site_name, ".tar.gz"))
  
  # Temporarily change working directory so files are stored without path prefixes
  old_wd <- setwd(outp_dir)
  on.exit(setwd(old_wd))
  
  # Create tar
  utils::tar(
    tarfile = tar_path,
    files = files_to_tar,
    compression = "gzip",
    tar = "internal" # cross-platform compatibility
  )

  if (then_rm_files) file.remove(files_to_tar)
    
  invisible(NULL)
}

# ══ AGGREGATOR ════════════════════════════════════════════════════════════════

#' @title Load and harmonise site payload archives
#'
#' @description
#' Extracts each site's \code{.tar.gz} archive, loads the static sufficient
#' statistics (\code{.static.rds}) and attaches the memory-mapped FBM objects
#' (\eqn{X^\top Y} and partial sums of \eqn{Y}), then transposes the resulting
#' list-of-sites into a list-of-components ready for the aggregation step.
#'
#' @details
#' Components that must be **identical** across sites (e.g. \code{terms},
#' \code{fs_template}) are validated and reduced to a single value; an error is
#' thrown if any site disagrees. Scalar components (e.g. \code{n_obs}) are
#' unlisted into named vectors. Matrix and FBM components are kept as lists of
#' length \eqn{K}.
#'
#' @param site_names Character vector of site identifiers.
#' @param input_dir Character. Directory containing \code{<site>.tar.gz} archives.
#' @param hemi Character. Hemisphere (\code{"lh"} or \code{"rh"}).
#' @param measure Character. Surface measure (derived from the model formula LHS).
#'
#' @return A named list with one element per payload component. Key elements:
#' \describe{
#'   \item{\code{n_obs}}{Named integer vector of per-site sample sizes.}
#'   \item{\code{terms}}{Character vector of fixed-effect term names (shared).}
#'   \item{\code{fs_template}}{FreeSurfer template name (shared).}
#'   \item{\code{XtX}}{List of \eqn{K} matrices \eqn{X_k^\top X_k} (\eqn{p \times p}).}
#'   \item{\code{X1}}{List of \eqn{K} vectors \eqn{X_k^\top \mathbf{1}_{n_k}} (length \eqn{p}).}
#'   \item{\code{XtY}}{List of \eqn{K} FBMs \eqn{X_k^\top Y_k} (\eqn{p \times V}).}
#'   \item{\code{psumsY}}{List of \eqn{K} FBMs with row 1 = \eqn{\mathbf{1}^\top Y_k}, row 2 = \eqn{Y_k^\top Y_k} (\eqn{2 \times V}).}
#' }
#'
#' @keywords internal
#' 
decompress_site_payloads <- function(site_names, input_dir, hemi, measure) {

  # Ensure path is correctly expanded if necessary
  input_dir <- normalizePath(input_dir)

  # Step 1: Extract and load each site information
  sites_list <- sapply(site_names, function(site) {

    site_results <- file.path(input_dir, site, paste(site, hemi, measure, sep = '.'))

    # Extract archive
    utils::untar(tarfile = file.path(input_dir, paste0(site, ".tar.gz")), 
                 exdir = file.path(input_dir, site))
    
    # Load static info
    sinfo <- readRDS(paste0(site_results, ".static.rds"))
    
    sinfo[['XtY']] <- bigstatsr::big_attach(paste0(site_results, ".XtY.rds"))
    sinfo[['psumsY']] <- bigstatsr::big_attach(paste0(site_results, ".psumsY.rds"))

    sinfo
  }, USE.NAMES = TRUE, simplify = FALSE)

  # Step 2: Transpose 
  component_names <- names(sites_list[[1]])

  site_payloads <- sapply(component_names, function(component) {
    values <- lapply(sites_list, `[[`, component)
    if (component %in% c('n_obs', 'n_good_vx', 'date_created', 'verywise_version')) {
      unlist(values)  # Convert list of "scalars" -> named vector
    } else if (component  %in% c('terms', 'fs_template')) {
      # They should all be identical
      first <- values[[1]]
      is_identical <- vapply(values[-1], identical, logical(1), y = first)
      if (!all(is_identical)) stop(paste(component, "must be consistent across sites"))
      first
    } else {
      values # Keep as list (e.g., for matrices like XtX)
    }
  }, USE.NAMES = TRUE, simplify = FALSE)
  
  # TODO: warn if verywise version is not the same?
  # if not warn, the matrices have all the same dimensions 

  # Display some info:
  print(as.data.frame(site_payloads[
    c('n_obs', 'fs_template', 'n_good_vx', 'verywise_version', 'date_created')]))
  
  site_payloads

}


#' @title Vertex-wise distributed LMM: per-chunk estimation
#'
#' @description
#' Core computational workhorse called by \code{\link{run_vw_fed_aggr}} for
#' each chunk of vertices. Given the pre-computed static lambda profile
#' \code{lambda_prof} and the chunk-specific sufficient statistics
#' (\eqn{X^\top Y}, \eqn{\mathbf{1}^\top Y}, \eqn{Y^\top Y}),
#' this function selects the optimal \eqn{\lambda} per vertex and returns
#' fixed-effect estimates, standard errors, p-values, and variance components.
#'
#' @details
#' For each candidate \eqn{\lambda} the dynamic quantities are:
#' \deqn{B(\lambda) = X^\top V^{-1} Y = \sum_k \bigl(X_k^\top Y_k - w_k \mathbf{x}_{k1} \mathbf{1}_k^\top Y_k \bigr)}
#' \deqn{C(\lambda) = Y^\top V^{-1} Y = \sum_k \bigl(Y_k^\top Y_k - w_k (\mathbf{1}_k^\top Y_k)^2 \bigr)}
#' where \eqn{w_k = \lambda / (1 + n_k \lambda)} is the site weight.
#' The profile (RE)ML log-likelihood is evaluated at each candidate \eqn{\lambda}
#' and the maximiser selected per vertex via \code{\link[base]{max.col}}.
#'
#' Standard errors are derived from \eqn{\text{Var}(\hat\beta) = \hat\sigma^2 A^{-1}},
#' with \eqn{A^{-1}} obtained via \code{\link[base]{chol2inv}} from the
#' pre-computed Cholesky factor.
#'
#' @param lambda_prof List. Output of the static lambda profile loop in
#'   \code{\link{run_vw_fed_aggr}}. Each element is a named list with fields
#'   \code{lambda}, \code{w} (site weights), \code{R} (Cholesky factor of
#'   \eqn{A}), \code{logdetA}, and \code{logdetV}.
#' @param X1 List of length \eqn{K}. Each element is the column vector of
#'   row-sums of the design matrix for site \eqn{k} (\eqn{X_k^\top \mathbf{1}_{n_k}}).
#' @param XtY List of length \eqn{K}. Each element is a \eqn{p \times C} matrix
#'   \eqn{X_k^\top Y_k} for the current chunk of \eqn{C} vertices.
#' @param sumY List of length \eqn{K}. Each element is a length-\eqn{C} vector
#'   \eqn{\mathbf{1}^\top Y_k} (column sums of the outcome chunk).
#' @param YtY List of length \eqn{K}. Each element is a length-\eqn{C} vector
#'   of element-wise squared column sums \eqn{Y_k^\top Y_k}.
#' @param n_terms Integer. Number of fixed-effect terms \eqn{p}.
#' @param chunk_size Integer. Number of vertices in the current chunk \eqn{C}.
#' @param n_sites Integer. Number of sites \eqn{K}.
#' @param REML Logical. Use REML (\code{TRUE}) or ML (\code{FALSE}) likelihood.
#' @param df Integer. Residual degrees of freedom (\eqn{N - p} for REML, \eqn{N} for ML).
#'
#' @return A named list with four matrices, all with \eqn{C} columns:
#' \describe{
#'   \item{\code{coef}}{\eqn{p \times C} fixed-effect estimates.}
#'   \item{\code{se}}{\eqn{p \times C} standard errors.}
#'   \item{\code{pval}}{\eqn{p \times C} two-sided p-values.}
#'   \item{\code{bw_var}}{\eqn{2 \times C} variance components; row 1 = \eqn{\hat\sigma^2_\varepsilon}, row 2 = \eqn{\hat\tau^2}.}
#' }
#'
#' @keywords internal
#' 
chunk_dlmm <- function(lambda_prof, X1, XtY, sumY, YtY, n_terms, chunk_size, n_sites, REML, df) {

  # (to preallocate dims)
  # Y_size <- ncol(XtY[[1]]) # effective chunk size ~ 1000 

  # Profile lambda (dynamic): B, C, log-likelihood
  
  lambda_profY <- lapply(lambda_prof, function(lambda) {
    # static elements: lambda, w, R, logdetA, logdetV

    site_weights <- as.list(lambda$w)

    # for each variance ratio λ, compute: XtV-1Y (Bmat) and YtV-1Y (Cvec)
    Bmat <- Reduce("+", Map(\(XtY, X1, sumY, wk) XtY - wk * (X1 %o% sumY),
                        XtY, X1, sumY, site_weights)) # dim: p x chunk_size
    
    Cvec <- Reduce("+", Map(\(YtY, sumY, wk) YtY - wk * sumY^2,
                        YtY, sumY, site_weights)) # dim: chunk

    beta <- vw_chol_solve(lambda$R, Bmat)  # dim: p x chunk_size

    # Residual Sum of Squares (= sigma2 * df at optimum sigma2)
    Q <- Cvec - colSums(beta * Bmat)  # dim: chunk_size 
    # Numerical safety
    Q[Q <= 0 | !is.finite(Q)] <- NA_real_

    # Profile log-likelihood (constants dropped) --> maximize this
    logLL <- vw_logLL(REML = REML, df = df, Q = Q, 
      lambda$logdetV, lambda$logdetA)

    c(lambda, list(logLL = logLL, Bmat = Bmat, Cvec = Cvec, beta = beta, Q = Q))
    # logLL
  })

  # Pick best grid lambda per vertex
  lambda_logLL <- vapply(lambda_profY, `[[`, numeric(chunk_size), "logLL")
  # plot(lambda_logLL[1, ]) # see tuning
  maxlogLL_idx <- max.col(lambda_logLL)
  
  # lambda_hat  <- all_lambdas[maxlogLL_idx] # dim: chunk_size

  # Optional small refinement: local 1D search around lam_hat using a tiny bracket
  # (kept simple and robust; still block-wise)

  # Finally, compute at lamda_hat: vertex-wise but in grouped lambdas to stay efficient
  # Group by unique lambdas (rounded) reuse A factorization.
  # lambda_key <- round(lambda_hat, digits = 12)
  unique_best_lambdas <- sort(unique(maxlogLL_idx))

  # Storage
  coef <- matrix(NA_real_, nrow = n_terms, ncol = chunk_size)
  se <- matrix(NA_real_, nrow = n_terms, ncol = chunk_size)
  pval <- matrix(NA_real_, nrow = n_terms, ncol = chunk_size)
  bw_var <- matrix(NA_real_, nrow = 2, ncol = chunk_size) # sigma2 and tau2

  for (idx in unique_best_lambdas) {

    vw_filter <- maxlogLL_idx == idx

    stats <- with(lambda_profY[[idx]], {
      # lambda_hat <- lambda # lambda$lambda to be clear 

      # SE: Var(beta) = sigma2 * A^{-1}; diagonal of A^{-1} from cho
      sigma2 <- Q[vw_filter] / df
      tau2 <- lambda * sigma2
      Ainv_diag <- diag(chol2inv(R))

      se <- sqrt(outer(Ainv_diag, sigma2, `*`))

      coef <- beta[, vw_filter, drop = FALSE] 
      tval <- coef / se
      pval <- 2 * stats::pt(abs(tval), df = df, lower.tail = FALSE)

      list(coef = coef, se = se, pval = pval, sigma2 = sigma2, tau2 = tau2)
    })

    coef[, vw_filter] <- stats$coef
    se[, vw_filter] <- stats$se
    pval[, vw_filter] <- stats$pval
    bw_var[, vw_filter] <- rbind(stats$sigma2, stats$tau2)
  }
  
  list(coef = coef, se = se, pval = pval, bw_var = bw_var)
}


#' @title Robust Cholesky decomposition with positive-definiteness check
#'
#' @description
#' A safe wrapper around \code{\link[base]{chol}} that returns \code{NULL}
#' instead of throwing an error when the input matrix \code{A} is not
#' numerically positive-definite (PD). An additional near-singularity check
#' rejects decompositions whose smallest diagonal pivot falls below
#' \code{tol}, preventing numerically unstable downstream solves.
#'
#' @param A Numeric square matrix. Should be symmetric PD (e.g. a Hessian
#'   \eqn{X^\top V^{-1} X}).
#' @param tol Numeric. Minimum acceptable diagonal pivot. Default: \code{1e-12}.
#'
#' @return The upper-triangular Cholesky factor \eqn{R} such that
#'   \eqn{R^\top R = A}, or \code{NULL} if \code{A} is not sufficiently PD.
#'
#' @keywords internal
#' 
vw_chol <- function(A, tol = 1e-12) {

  R <- tryCatch(chol(A), error = function(e) NULL)
  
  # Failure: A is not PD ("Error: leading minor of order ... is not positive definite")
  if (is.null(R)) return(NULL) 
  
  # Extract the diagonal elements of the computed factor R
  d <- diag(R)

  # Check for:
  #  - empty matrices (edge case)
  #  - Inf or NaN values, which indicate numerical instability
  #  - Near-singularity. Even if chol() succeeds, very small diagonal entries imply the
  #    matrix is nearly singular (condition number is high). If the smallest pivot is 
  #    below tol (default 1e-12), treat the decomposition as failed to avoid numerical
  #    errors in downstream inverse calculations
  if (length(d) == 0L || any(!is.finite(d)) || min(d) < tol) return(NULL)
  
  # Happy path
  R
}

#' @title Solve a linear system given a Cholesky factor
#'
#' @description
#' Solves \eqn{A X = B} for \eqn{X} given the upper-triangular Cholesky factor
#' \eqn{R} of \eqn{A} (i.e. \eqn{R^\top R = A}), using forward and backward
#' substitution. This is numerically preferable to forming \eqn{A^{-1}} explicitly.
#'
#' @param R Upper-triangular matrix. Cholesky factor of \eqn{A} as returned by
#'   \code{\link{vw_chol}} or \code{\link[base]{chol}}.
#' @param B Numeric matrix (or vector). Right-hand side; must have the same
#'   number of rows as \code{R}.
#'
#' @return Numeric matrix \eqn{X = A^{-1} B} of the same dimensions as \code{B}.
#'
#' @keywords internal
#' 
vw_chol_solve <- function(R, B) {
  # Solve A X = B given chol(A)=R'R; base R: chol2inv(R) %*% B is ok for small p.
  # Use forwardsolve/backsolve for better numeric behavior.
  y <- forwardsolve(t(R), B, upper.tri = FALSE, transpose = FALSE)
  backsolve(R, y, upper.tri = TRUE, transpose = FALSE)
}

#' @title Profile (RE)ML log-likelihood for a single variance ratio
#'
#' @description
#' Computes the profile (RE)ML log-likelihood (up to constants) for a given
#' candidate \eqn{\lambda} and a vector of residual sum-of-squares values
#' \eqn{Q} across vertices. This is evaluated inside the \eqn{\lambda} grid
#' search to select the optimal variance ratio per vertex.
#'
#' @details
#' **REML** profile log-likelihood (constants dropped):
#' \deqn{\ell_{\text{REML}} = -\frac{1}{2} \left[ (N-p) \log\!\frac{Q}{N-p} + \log|V| + \log|A| \right]}
#'
#' **ML** profile log-likelihood (constants dropped):
#' \deqn{\ell_{\text{ML}} = -\frac{1}{2} \left[ N \log\!\frac{Q}{N} + \log|V| \right]}
#'
#' where \eqn{\log|V|} and \eqn{\log|A|} are the pre-computed log-determinants
#' passed as \code{logdetV} and \code{logdetA}.
#'
#' @param REML Logical. If \code{TRUE} (default), evaluate the REML objective;
#'   otherwise the ML objective.
#' @param df Integer. Degrees of freedom: \eqn{N - p} (REML) or \eqn{N} (ML).
#' @param Q Numeric vector of length \eqn{C} (chunk size). Residual sum of
#'   squares \eqn{Q = Y^\top V^{-1} Y - \hat\beta^\top X^\top V^{-1} Y} per vertex.
#'   Vertices with \code{Q <= 0} or non-finite values are pre-set to \code{NA}
#'   by the caller.
#' @param logdetV Numeric scalar. \eqn{\log |V(\lambda)|}, computed as
#'   \eqn{\sum_k \log(1 + n_k \lambda)}.
#' @param logdetA Numeric scalar. \eqn{\log |A(\lambda)|}, computed as
#'   \eqn{2 \sum_j \log R_{jj}} from the Cholesky factor.
#'
#' @return Numeric vector of length \eqn{C}; the profile log-likelihood value
#'   at this \eqn{\lambda} for each vertex in the chunk. Vertices with
#'   \code{NA} in \code{Q} propagate \code{NA}.
#'
#' @keywords internal
#' 
vw_logLL <- function(REML = TRUE, df, Q, logdetV, logdetA) {
    
    if (REML) {
      # REML: -0.5 * [ (N-p)*log(Q/(N-p)) + sum log(1+n_i*lam) + log|A| ]

      logLL <- -0.5 * (df * log(Q / df) + logdetV + logdetA)

    } else {
      # ML: -0.5 * [ N*log(Q/N) + sum log(1+n_i*lam) ]

      logLL <- -0.5 * (df * log(Q / df) + logdetV) 
    }
  
  logLL
}
