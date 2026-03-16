#' @title
#' Simulate a distributed (multi-site) cross-sectional brain surface dataset
#'
#' @description
#' Generates a synthetic distributed dataset where each site's data lives in its
#' own subdirectory, mirroring a real multi-site setup where no site can access
#' another's raw data. For each site \eqn{k} and each vertex \eqn{v}, the generative
#'  model is:
#'
#' \deqn{
#'   y_{ikv} = \mu_v
#'     + \sum_{j} \beta_j \cdot x_{ijk}
#'     + u_{kv}
#'     + \varepsilon_{ikv}
#' }
#'
#' where \eqn{\mu_v} is a vertex-specific intercept drawn once from
#' \eqn{\mathcal{N}(\code{vw\_mean},\, \code{vw\_sd}^2)},
#' \eqn{\beta_j} are fixed effects shared across all vertices and sites
#' (supplied via \code{betas}),
#' \eqn{u_{kv} \sim \mathcal{N}(0, \tau^2_v)} is a site-level random intercept,
#' and \eqn{\varepsilon_{ikv} \sim \mathcal{N}(0, \sigma^2_v)} is residual noise.
#'
#' The output folder layout is:
#' \preformatted{
#' <path>/
#'   <site1>/
#'     phenotype.csv
#'     sub-001/surf/<hemi>.<measure>.<fwhmc>.<fs_template>.mgh
#'     sub-002/surf/...
#'   <site2>/
#'     phenotype.csv
#'     ...
#' }
#'
#' Each site folder can therefore be passed directly to
#' \code{\link{run_vw_fed_local}} as the \code{subj_dir} and the matching
#' \file{phenotype.csv} as \code{pheno}.
#'
#' @param path Character string. Root directory for all site sub-folders.
#'   Created recursively if it does not exist.
#' @param site_sizes Named integer vector. Names become site/folder names;
#'   values give the number of subjects at each site.
#'   Default: \code{c(site1 = 80, site2 = 120, site3 = 60)}.
#' @param betas Named numeric vector of fixed-effect coefficients. Each name
#'   must correspond to a covariate that will be simulated in the phenotype
#'   data. Currently supported covariate names:
#'   \describe{
#'     \item{\code{age}}{Continuous, drawn from
#'       \eqn{\mathcal{N}(30, 5^2)} and then z-scored within each site.}
#'     \item{\code{sex}}{Binary (0/1), drawn from \eqn{\text{Bernoulli}(0.5)}.}
#'   }
#'   Example: \code{c(age = 0.3, sex = -0.2)}.
#'   Set any coefficient to \code{0} to include the covariate in the design
#'   matrix without injecting a signal (useful for null-effect benchmarks).
#' @param tau2 Numeric scalar or length-\eqn{V_{\text{roi}}} vector.
#'   **Between-site variance** of the random intercept \eqn{u_{kv}}.
#'   A value of \code{0} collapses the model to a fixed-effects OLS.
#'   Default: \code{0.5}.
#' @param sigma2 Numeric scalar or length-\eqn{V_{\text{roi}}} vector.
#'   **Within-site** residual variance \eqn{\varepsilon_{ikv}}.
#'   Default: \code{1}.
#'   \code{tau2} and \code{sigma2} together define the intra-class correlation
#'   \eqn{\text{ICC} = \tau^2 / (\tau^2 + \sigma^2)}.
#' @param fs_template Character string. FreeSurfer template; determines the
#'   total number of vertices in each \code{.mgh} file. Options:
#'   \itemize{
#'     \item \code{"fsaverage"}  = 163842 vertices
#'     \item \code{"fsaverage6"} = 40962 vertices
#'     \item \code{"fsaverage5"} = 10242 vertices
#'     \item \code{"fsaverage4"} = 2562 vertices
#'     \item \code{"fsaverage3"} = 642 vertices
#'   }
#'   Default: \code{"fsaverage"}.
#' @param roi_subset Character vector of ROI names used to restrict the active
#'   vertices (the rest are set to \code{0} and excluded from analysis).
#'   Vertex locations are read from the FreeSurfer annotation files stored in
#'   \file{R/sysdata.rda}.
#'   Default: \code{c("temporalpole", "frontalpole", "entorhinal")}.
#' @param location_association Character vector (optional). If supplied, the
#'   signal encoded in \code{betas} is injected \emph{only} within these ROIs;
#'   the remaining active vertices (\code{roi_subset}) are generated under a
#'   null model (\eqn{\beta = 0}). Useful for testing spatial localisation of
#'   discovered effects.
#' @param measure Character string. Surface measure; used for file naming only.
#'   Default: \code{"thickness"}.
#' @param hemi Character string. Hemisphere: \code{"lh"} or \code{"rh"}.
#'   Default: \code{"lh"}.
#' @param fwhmc Character string. Smoothing label; used for file naming only.
#'   Default: \code{"fwhm10"}.
#' @param vw_mean Numeric scalar. Mean of the vertex-specific intercept
#'   \eqn{\mu_v \sim \mathcal{N}(\code{vw\_mean},\, \code{vw\_sd}^2)}.
#'   For cortical thickness a realistic value is around \code{2.5} mm.
#'   Default: \code{2.5}.
#' @param vw_sd Numeric scalar. Standard deviation of the vertex-specific
#'   intercept distribution. Default: \code{0.5}.
#' @param overwrite Logical. If \code{FALSE} and \file{phenotype.csv} already
#'   exists in a site folder, the existing file is re-used and no new brain
#'   surface files are written for that site. Default: \code{TRUE}.
#' @param seed Integer. Random seed for reproducibility. Default: \code{3108}.
#' @param verbose Logical. Print progress messages. Default: \code{TRUE}.
#'
#' @return
#' Invisibly returns a list with the ground-truth parameters used during data
#' generation:
#' \describe{
#'   \item{\code{beta0}}{Length-\eqn{V_{\text{roi}}} vector of simulated
#'     vertex-specific intercepts \eqn{\mu_v}.}
#'   \item{\code{betas}}{The \code{betas} argument as supplied.}
#'   \item{\code{tau2}}{Length-\eqn{V_{\text{roi}}} vector of between-site
#'     variances (after recycling).}
#'   \item{\code{sigma2}}{Length-\eqn{V_{\text{roi}}} vector of residual
#'     variances (after recycling).}
#'   \item{\code{u}}{Numeric matrix \eqn{K \times V_{\text{roi}}} of realised
#'     site random intercepts.}
#'   \item{\code{icc}}{Length-\eqn{V_{\text{roi}}} vector of theoretical
#'     ICCs: \eqn{\tau^2_v / (\tau^2_v + \sigma^2_v)}.}
#'   \item{\code{signal_vertices}}{Logical vector of length \eqn{n_{\text{verts}}}
#'     indicating which global vertex indices received a non-zero fixed effect
#'     (i.e. the intersection of \code{roi_subset} and
#'     \code{location_association}, if supplied).}
#' }
#' Data and phenotype files are written to \code{path} as a side-effect.
#'
#' @seealso
#' \code{\link{run_vw_fed_local}} for the analysis function this feeds into,
#' \code{\link{simulate_longit_dataset}} for the longitudinal (multi-session) variant.
#'
#' @examples
#' \dontrun{
#' truth <- simulate_fed_dataset(
#'   path = tempfile("fed_sim_"),
#'   site_sizes = c(site1 = 80, site2 = 120, site3 = 60),
#'   betas = c(age = 0.3, sex = -0.2),
#'   tau2 = 0.5,
#'   sigma2 = 1,
#'   fs_template = "fsaverage5"  # fast; use "fsaverage" for final analyses
#' )
#'
#' # Ground-truth ICC summary across active vertices
#' summary(truth$icc)
#'
#' # Recovered site random intercepts (K x V_roi matrix)
#' dim(truth$u)
#' }
#'
#' @author Serena Defina, 2026.
#'
#' @export
simulate_distrib_dataset <- function(path,
    site_sizes = c(site1 = 80L, site2 = 120L, site3 = 60L),
    betas = c(sex = -0.2, age = 0.3),
    tau2 = 0.01,
    sigma2 = 0.01,
    fs_template = "fsaverage",
    roi_subset = c("temporalpole", "frontalpole", "entorhinal"),
    location_association = NULL,
    measure = "thickness",
    hemi = "lh",
    fwhmc = "fwhm10",
    vw_mean = 5,
    vw_sd = 0.5,
    overwrite = TRUE,
    seed = 3108,
    verbose = TRUE) {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    is.numeric(site_sizes), all(site_sizes > 0), !is.null(names(site_sizes)),
    is.numeric(betas),      !is.null(names(betas)),
    all(names(betas) %in% c("age", "sex")),
    is.numeric(tau2),       all(tau2   >= 0),
    is.numeric(sigma2),     all(sigma2 >  0)
  )

  K  <- length(site_sizes)
  site_names <- names(site_sizes)
  mgh_fname  <- paste(hemi, measure, fwhmc, fs_template, "mgh", sep = ".")
  n_verts  <- count_vertices(fs_template)

  check_path(path, create_if_not = TRUE)
  set.seed(seed)

  # ── Locate active vertices from annotation files ──────────────────────────
 
  roi_locs <- locate_roi(
    rois    = roi_subset,
    n_verts = n_verts,
    hemi    = hemi,
    verbose = verbose
  )
  V_roi <- sum(roi_locs)

  # Vertices that will carry the fixed-effect signal (may be a strict subset)
  if (!is.null(location_association)) {
    signal_locs <- locate_roi(
      rois    = location_association,
      n_verts = n_verts,
      hemi    = hemi,
      verbose = verbose
    )
    # Signal vertices must lie within roi_subset
    signal_locs <- signal_locs & roi_locs
  } else {
    signal_locs <- roi_locs
  }

  # ── Ground-truth parameters (drawn once, shared across sites) ────────────
  tau2   <- rep_len(tau2,   V_roi)
  sigma2 <- rep_len(sigma2, V_roi)
  icc    <- tau2 / (tau2 + sigma2)

  # Vertex-specific intercepts μ_v ~ N(vw_mean, vw_sd^2), length V_roi
  beta0 <- stats::rnorm(V_roi, mean = vw_mean, sd = vw_sd)

  # Site random intercepts u[k, v] ~ N(0, tau2[v]),  dim: K x V_roi
  u <- matrix(
    stats::rnorm(K * V_roi, mean = 0, sd = rep(sqrt(tau2), each = K)),
    nrow = K, ncol = V_roi
  )

  # ── Per-site data generation ──────────────────────────────────────────────
  for (k in seq_len(K)) {

    site   <- site_names[k]
    n_k    <- site_sizes[k]
    site_dir <- file.path(path, site)
    check_path(site_dir, create_if_not = TRUE)

    pheno_path <- file.path(site_dir, "phenotype.csv")

    # ── Phenotype ────────────────────────────────────────────────────────────
    if (!overwrite && file.exists(pheno_path)) {
      vw_message(
        " * [", site, "] re-using phenotype file from: ",
        format(file.info(pheno_path)[, "mtime"]),
        verbose = verbose
      )
      pheno <- utils::read.csv(pheno_path)

    } else {
      pheno <- simulate_pheno_data(
        site_name = site,
        n         = n_k,
        has_age   = "age" %in% names(betas),
        has_sex   = "sex" %in% names(betas)
      )
      utils::write.csv(pheno, pheno_path, row.names = FALSE)
    }

    vw_message(
      " * [", site, "] writing ", n_k, " surface files...",
      verbose = verbose
    )

    # Pre-compute fixed-effect contribution for each subject: n_k x V_roi
    # Only non-zero over signal_locs (columns indexed within roi_locs)
    signal_within_roi <- signal_locs[roi_locs]  # logical length V_roi

    fe_all <- matrix(0, nrow = n_k, ncol = V_roi)
    for (nm in names(betas)) {
      fe_all[, signal_within_roi] <-
        fe_all[, signal_within_roi] +
        betas[[nm]] * pheno[[nm]]  # outer product: n_k x sum(signal_within_roi)
    }

    # ── Write one .mgh per subject ───────────────────────────────────────────
    for (i in seq_len(n_k)) {

      sub_dir <- file.path(site_dir, pheno$folder_id[i])
      dir.create(file.path(sub_dir, "surf"), recursive = TRUE,
                 showWarnings = FALSE)

      vw_data <- numeric(n_verts)         # full surface, mostly 0
      eps <- stats::rnorm(V_roi, 0, sqrt(sigma2))

      vw_data[roi_locs] <- beta0 +                 # μ_v
                           fe_all[i, ] +           # Σ β_j x_ij (0 outside signal)
                           u[k, ] +                # u_kv
                           eps                     # ε_ikv

      vw_data[vw_data < 0] <- 0.01                 # non-negative surface metric

      save.mgh(as.mgh(vw_data), file.path(sub_dir, "surf", mgh_fname))
    }
  }

  vw_message("Done.", verbose = verbose)

  invisible(list(beta0 = beta0, betas = betas,
                 tau2 = tau2, sigma2 = sigma2, u = u, icc = icc,
                 signal_vertices = signal_locs))
}


simulate_pheno_data <- function(site_name, n, has_age, has_sex) {

  pheno <- data.frame(
    folder_id = sprintf("sub-%03d", seq_len(n))
  )

  if (has_age) {
    raw_age    <- stats::rnorm(n, mean = 30, sd = 5)
    pheno$age  <- (raw_age - mean(raw_age)) / stats::sd(raw_age) # z-scored within site
  }

  if (has_sex) {
    pheno$sex <- stats::rbinom(n, 1L, 0.5)
  }

  pheno
}
