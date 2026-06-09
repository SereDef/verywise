# ── Supported covariate names in `betas` ──────────────────────────────────────
#
#  "age"     Continuous. Baseline age ~ N(mean_age, sd_age^2), z-scored within
#            site. At each follow-up the value is incremented by interval_years
#            plus a small N(0, 0.5^2) jitter.
#  "sex"     Binary 0 / 1 ~ Bernoulli(0.5). Constant across sessions.
#  "time"    0-based integer session index (0, 1, 2, …). Encodes the
#            within-person linear time trend.
#  "wisdom"  Continuous nuisance covariate ~ N(10, 3^2), z-scored within site.
#            Constant across sessions.
#
# Any name not in the list above is ignored during phenotype generation but
# you may still set its beta to 0 for a null-effect benchmark.
# ─────────────────────────────────────────────────────────────────────────────


# ══════════════════════════════════════════════════════════════════════════════
#  simulate_longit_dataset2()
# ══════════════════════════════════════════════════════════════════════════════

#' @title
#' Simulate a longitudinal brain surface dataset with associated phenotype data
#'
#' @description
#' Generates a synthetic longitudinal dataset for multiple sites/cohorts, each
#' with multiple timepoints/sessions per subject.  The generative model for
#' vertex \eqn{v}, subject \eqn{i} nested in site \eqn{k}, at session \eqn{t}
#' is:
#'
#' \deqn{
#'   y_{itkv} = \mu_v
#'     + \sum_j \beta_j \cdot x_{ijtk}
#'     + b_{iv}
#'     + b^{\mathrm{time}}_{iv} \cdot t
#'     + u_{kv}
#'     + \varepsilon_{itkv}
#' }
#'
#' where
#' \eqn{\mu_v \sim \mathcal{N}(\code{vw\_mean},\, \code{vw\_sd}^2)} is a
#' vertex-specific intercept drawn once and shared across all subjects and sites,
#' \eqn{\beta_j} are fixed effects supplied via \code{betas},
#' \eqn{b_{iv} \sim \mathcal{N}(0,\, \code{subj\_intercept\_sd}^2)} is a
#' subject-level random intercept,
#' \eqn{b^{\mathrm{time}}_{iv} \sim \mathcal{N}(0,\, \code{subj\_slope\_sd}^2)}
#' is a subject-level random slope for time (set \code{subj\_slope\_sd = 0} to
#' suppress it),
#' \eqn{u_{kv} \sim \mathcal{N}(0,\, \code{site\_sd}^2)} is a site-level random
#' intercept,
#' and \eqn{\varepsilon_{itkv} \sim \mathcal{N}(0,\, \code{sigma2})} is residual
#' noise.
#'
#' The function writes:
#' \itemize{
#'   \item Brain surface data in FreeSurfer \code{.mgh} format organised in a
#'         \pkg{verywise} folder structure (see vignettes for details).
#'   \item A phenotype \code{data.frame} saved as \file{phenotype.csv} in
#'         \code{path}.
#' }
#'
#' @param path Character string. Root directory for the dataset.  Created
#'   recursively if it does not exist.
#' @param data_structure Named list defining cohorts/sites.  Each element is a
#'   list with:
#'   \describe{
#'     \item{\code{"sessions"}}{Character vector of session labels.}
#'     \item{\code{"n_subjects"}}{Integer number of subjects.}
#'   }
#' @param betas Named numeric vector of fixed-effect coefficients.  Names
#'   identify the covariates; currently supported names are \code{"age"},
#'   \code{"sex"}, \code{"time"}, and \code{"wisdom"}.  Set a coefficient to
#'   \code{0} to include the covariate without injecting a signal (useful for
#'   null-effect benchmarks).
#'   Example: \code{c(age = -0.1, sex = 0.05, time = -0.2)}.
#' @param sigma2 Numeric scalar (> 0). Within-subject residual variance
#'   \eqn{\varepsilon_{itkv}}.  Default: \code{0.25}.
#' @param subj_intercept_sd Numeric scalar (\eqn{\geq 0}).  Standard deviation
#'   of the subject-level random intercept \eqn{b_{iv}}.  Default: \code{0.2}.
#' @param subj_slope_sd Numeric scalar (\eqn{\geq 0}).  Standard deviation of
#'   the subject-level random slope for \code{time}.  Set to \code{0} to remove
#'   random slopes entirely.  Default: \code{0.1}.
#' @param site_sd Numeric scalar (\eqn{\geq 0}).  Standard deviation of the
#'   site-level random intercept \eqn{u_{kv}}.  Default: \code{0.1}.
#' @param vw_mean Numeric scalar.  Mean of the vertex-specific intercept
#'   \eqn{\mu_v}.  For cortical thickness a realistic value is around
#'   \code{2.5} mm.  Default: \code{2.5}.
#' @param vw_sd Numeric scalar.  Standard deviation of \eqn{\mu_v}.
#'   Default: \code{0.5}.
#' @param mean_age Numeric scalar.  Population mean of baseline age.
#'   Default: \code{20}.
#' @param sd_age Numeric scalar.  Population standard deviation of baseline
#'   age.  Default: \code{2}.
#' @param interval_years Numeric scalar.  Mean number of years between
#'   consecutive sessions.  Default: \code{2}.
#' @param fs_template Character string.  FreeSurfer template; determines the
#'   number of vertices per surface file.  Options:
#'   \itemize{
#'     \item \code{"fsaverage"}  — 163 842 vertices (default)
#'     \item \code{"fsaverage6"} — 40 962 vertices
#'     \item \code{"fsaverage5"} — 10 242 vertices
#'     \item \code{"fsaverage4"} —  2 562 vertices
#'     \item \code{"fsaverage3"} —    642 vertices
#'   }
#' @param roi_subset Character vector of ROI names that define the active
#'   vertices.  The remaining vertices are set to \code{0} and excluded from
#'   analysis.  Default: \code{c("temporalpole", "frontalpole", "entorhinal")}.
#' @param location_association Optional character vector.  If supplied, the
#'   fixed-effect signal from \code{betas} is injected \emph{only} within
#'   these ROIs; the remaining active vertices in \code{roi_subset} are
#'   generated under a null model (\eqn{\beta = 0}).  Useful for testing
#'   spatial localisation of discovered effects.
#' @param measure Character string.  Surface measure label; used for file
#'   naming only.  Default: \code{"thickness"}.
#' @param hemi Character string.  Hemisphere: \code{"lh"} or \code{"rh"}.
#'   Default: \code{"lh"}.
#' @param fwhmc Character string.  Smoothing label; used for file naming only.
#'   Default: \code{"fwhm10"}.
#' @param overwrite Logical.  If \code{FALSE} and \file{phenotype.csv} already
#'   exists in \code{path}, it is reused and no new brain files are written.
#'   Default: \code{TRUE}.
#' @param seed Integer.  Random seed for reproducibility.  Default: \code{3108}.
#' @param verbose Logical.  Print progress messages.  Default: \code{TRUE}.
#'
#' @return
#' Invisibly returns a named list of ground-truth parameters:
#' \describe{
#'   \item{\code{beta0}}{Length-\eqn{V_{\mathrm{roi}}} vector of simulated
#'     vertex-specific intercepts \eqn{\mu_v}.}
#'   \item{\code{betas}}{The \code{betas} argument as supplied.}
#'   \item{\code{sigma2}}{Scalar residual variance as supplied.}
#'   \item{\code{subj_intercept_sd}, \code{subj_slope_sd},
#'     \code{site_sd}}{Random-effect standard deviations as supplied.}
#'   \item{\code{u_site}}{Length-\eqn{K} vector of realised site random
#'     intercepts.}
#'   \item{\code{signal_vertices}}{Logical vector of length
#'     \eqn{n_{\mathrm{verts}}} indicating which global vertex indices carry a
#'     non-zero fixed effect.}
#' }
#' Data and phenotype files are written to \code{path} as a side-effect.
#'
#' @seealso
#' \code{\link{simulate_long_pheno_data2}},
#' \code{\link{simulate_freesurfer_data2}},
#' \code{\link{simulate_distrib_dataset}}
#'
#' @examples
#' truth <- simulate_longit_dataset2(
#'   path              = tempfile("longit_sim"),
#'   betas             = c(age = -0.1, sex = 0.05, time = -0.2),
#'   sigma2            = 0.25,
#'   subj_intercept_sd = 0.3,
#'   subj_slope_sd     = 0.05,
#'   site_sd           = 0.1,
#'   fs_template       = "fsaverage3"
#' )
#' truth$betas
#' sum(truth$signal_vertices)
#'
#' @author Serena Defina, 2024.
#'
#' @export
simulate_longit_dataset2 <- function(
    path,
    data_structure = list(
      "cohort1" = list("sessions" = c("01", "02", "03"), "n_subjects" = 10),
      "cohort2" = list("sessions" = c("01", "02"),       "n_subjects" = 20)
    ),
    betas             = c(age = -0.1, sex = 0.05, time = -0.2),
    sigma2            = 0.25,
    subj_intercept_sd = 0.2,
    subj_slope_sd     = 0.1,
    site_sd           = 0.1,
    vw_mean           = 2.5,
    vw_sd             = 0.5,
    mean_age          = 20,
    sd_age            = 2,
    interval_years    = 2,
    fs_template       = "fsaverage",
    roi_subset        = c("temporalpole", "frontalpole", "entorhinal"),
    location_association = NULL,
    measure           = "thickness",
    hemi              = "lh",
    fwhmc             = "fwhm10",
    overwrite         = TRUE,
    seed              = 3108,
    verbose           = TRUE) {

  vw_init_message("Simulating longitudinal dataset", verbose = verbose)

  # ── Input validation ─────────────────────────────────────────────────────
  stopifnot(
    is.list(data_structure),       !is.null(names(data_structure)),
    is.numeric(betas),             !is.null(names(betas)),
    all(names(betas) %in% c("age", "sex", "time", "wisdom")),
    is.numeric(sigma2),            sigma2 > 0,
    is.numeric(subj_intercept_sd), subj_intercept_sd >= 0,
    is.numeric(subj_slope_sd),     subj_slope_sd     >= 0,
    is.numeric(site_sd),           site_sd           >= 0
  )

  check_path(path, create_if_not = TRUE)
  set.seed(seed)

  pheno_path <- file.path(path, "phenotype.csv")

  # ── Phenotype ────────────────────────────────────────────────────────────
  if (!overwrite && file.exists(pheno_path)) {
    vw_message(
      "i" = "re-using existing phenotype file from: {format(file.info(pheno_path)[, 'mtime'])}",
      verbose = verbose
    )
    pheno <- utils::read.csv(pheno_path)

  } else {
    pheno <- simulate_long_pheno_data2(
      data_structure = data_structure,
      betas          = betas,
      mean_age       = mean_age,
      sd_age         = sd_age,
      interval_years = interval_years,
      seed           = seed,
      verbose        = verbose
    )
    utils::write.csv(pheno, pheno_path, row.names = FALSE)
  }

  # ── Brain surface data ───────────────────────────────────────────────────
  truth <- simulate_freesurfer_data2(
    path                 = path,
    data_structure       = data_structure,
    pheno                = pheno,
    betas                = betas,
    sigma2               = sigma2,
    subj_intercept_sd    = subj_intercept_sd,
    subj_slope_sd        = subj_slope_sd,
    site_sd              = site_sd,
    vw_mean              = vw_mean,
    vw_sd                = vw_sd,
    fs_template          = fs_template,
    roi_subset           = roi_subset,
    location_association = location_association,
    measure              = measure,
    hemi                 = hemi,
    fwhmc                = fwhmc,
    seed                 = seed,
    verbose              = verbose
  )

  vw_message("Done! :)", type = "step", verbose = verbose)
  invisible(truth)
}


# ══════════════════════════════════════════════════════════════════════════════
#  simulate_long_pheno_data2()
# ══════════════════════════════════════════════════════════════════════════════

#' @title
#' Simulate longitudinal phenotype data
#'
#' @description
#' Generates synthetic phenotype data in long format for multiple cohorts and
#' multiple sessions per subject.  The set of columns produced is driven by
#' \code{betas}: only the covariates whose names appear in \code{betas} are
#' generated, plus mandatory bookkeeping columns.
#'
#' Covariate generation rules:
#' \itemize{
#'   \item \code{age} — Baseline age
#'     \eqn{\sim \mathcal{N}(\code{mean\_age},\, \code{sd\_age}^2)},
#'     z-scored within site.  At each follow-up session the value is
#'     incremented by \code{interval_years} plus
#'     \eqn{\mathcal{N}(0,\, 0.5^2)} jitter.
#'   \item \code{sex} — Binary 0 / 1 \eqn{\sim \mathrm{Bernoulli}(0.5)};
#'     constant across sessions.
#'   \item \code{time} — 0-based integer session index
#'     (\code{0, 1, 2, \ldots}); enters the model as a within-person linear
#'     time trend.
#'   \item \code{wisdom} — Continuous nuisance covariate
#'     \eqn{\sim \mathcal{N}(10,\, 3^2)}, z-scored within site; constant
#'     across sessions.
#' }
#'
#' @inheritParams simulate_longit_dataset2
#'
#' @return
#' A \code{data.frame} in long format with mandatory columns \code{site},
#' \code{id}, \code{session} (original character label), \code{time}
#' (0-based integer index), \code{folder_id}, followed by any covariate
#' columns implied by \code{betas}.
#'
#' @seealso
#' \code{\link{simulate_longit_dataset2}}, \code{\link{simulate_freesurfer_data2}}
#'
#' @author Serena Defina, 2024.
#'
#' @export
simulate_long_pheno_data2 <- function(
    data_structure = list(
      "cohort1" = list("sessions" = c("01", "02"), "n_subjects" = 100),
      "cohort2" = list("sessions" = c("01", "02"), "n_subjects" = 150)
    ),
    betas          = c(age = -0.1, sex = 0.05, time = -0.2),
    mean_age       = 20,
    sd_age         = 2,
    interval_years = 2,
    seed           = 3108,
    verbose        = TRUE) {

  vw_message("* creating phenotype file...", verbose = verbose)
  set.seed(seed)

  # Which covariates are needed?
  has_age    <- "age"    %in% names(betas)
  has_sex    <- "sex"    %in% names(betas)
  has_wisdom <- "wisdom" %in% names(betas)
  # "time" is handled as a pure index column — always produced if data has
  # multiple sessions; its beta only affects brain data generation.

  site_frames <- lapply(names(data_structure), function(site) {

    n   <- data_structure[[site]]$n_subjects
    ss  <- data_structure[[site]]$sessions
    T_  <- length(ss)

    # ── Time-invariant baseline (one row per subject) ──────────────────────
    baseline <- data.frame(id = seq_len(n))

    if (has_sex) {
      baseline$sex <- stats::rbinom(n, 1L, 0.5)
    }
    if (has_age) {
      raw_age         <- stats::rnorm(n, mean = mean_age, sd = sd_age)
      baseline$age_bl <- (raw_age - mean(raw_age)) / stats::sd(raw_age)  # z-score within site
    }
    if (has_wisdom) {
      raw_wis         <- stats::rnorm(n, mean = 10, sd = 3)
      baseline$wisdom <- (raw_wis - mean(raw_wis)) / stats::sd(raw_wis)
    }

    # ── Expand to long format (one row per subject × session) ──────────────
    long <- do.call(rbind, lapply(seq_along(ss), function(t_idx) {
      row          <- baseline
      row$session  <- ss[t_idx]
      row$time     <- t_idx - 1L        # 0-based integer, matches brain-data loop

      if (has_age) {
        # Realistic longitudinal age: baseline + elapsed time + small noise
        jitter  <- stats::rnorm(n, mean = 0, sd = 0.5)
        row$age <- row$age_bl + (t_idx - 1L) * interval_years + jitter
        row$age_bl <- NULL              # drop helper; only 'age' is needed downstream
      }
      row
    }))

    long$site      <- site
    long$folder_id <- file.path(
      site,
      paste0("sub-", long$id, "_ses-", long$session)
    )
    long
  })

  long_df <- dplyr::bind_rows(site_frames)

  # Enforce a clean, readable column order
  id_cols    <- c("site", "id", "session", "time", "folder_id")
  covar_cols <- intersect(c("sex", "age", "wisdom"), names(long_df))
  long_df[, c(id_cols, covar_cols)]
}


# ══════════════════════════════════════════════════════════════════════════════
#  simulate_freesurfer_data2()
# ══════════════════════════════════════════════════════════════════════════════

#' @title
#' Simulate longitudinal FreeSurfer vertex-wise data
#'
#' @description
#' Writes one \code{.mgh} file per (subject, session) pair following the
#' generative model documented in \code{\link{simulate_longit_dataset2}}.
#' This function is exported so advanced users can generate brain surface data
#' from a pre-existing phenotype \code{data.frame} (e.g., after multiple
#' imputation).
#'
#' Random effects drawn per call (all seeded via \code{seed}):
#' \itemize{
#'   \item Vertex-specific intercepts \eqn{\mu_v} — drawn once, shared across
#'     sites and subjects.
#'   \item Site random intercepts \eqn{u_k} — one per site.
#'   \item Subject random intercepts \eqn{b_i} and slopes
#'     \eqn{b^{\mathrm{time}}_i} — one pair per subject per site.
#' }
#'
#' @inheritParams simulate_longit_dataset2
#' @param pheno \code{data.frame} as returned by
#'   \code{\link{simulate_long_pheno_data2}}, containing at minimum the columns
#'   \code{site}, \code{id}, \code{time} (0-based integer), \code{folder_id},
#'   plus any covariate columns named in \code{betas}.
#'
#' @return
#' Invisibly returns a named list of ground-truth simulation parameters:
#' \describe{
#'   \item{\code{beta0}}{Length-\eqn{V_{\mathrm{roi}}} vector of
#'     vertex-specific intercepts \eqn{\mu_v}.}
#'   \item{\code{betas}}{The \code{betas} argument as supplied.}
#'   \item{\code{sigma2}, \code{subj_intercept_sd}, \code{subj_slope_sd},
#'     \code{site_sd}}{Variance parameters as supplied.}
#'   \item{\code{u_site}}{Length-\eqn{K} vector of realised site random
#'     intercepts.}
#'   \item{\code{signal_vertices}}{Logical vector (length =
#'     total vertices) indicating which vertices carry a non-zero fixed effect.}
#' }
#'
#' @seealso \code{\link{simulate_longit_dataset2}}
#'
#' @author Serena Defina, 2024.
#'
#' @export
simulate_freesurfer_data2 <- function(
    path,
    data_structure    = list(
      "cohort1" = list("sessions" = c("01", "02"), "n_subjects" = 100),
      "cohort2" = list("sessions" = c("01", "02"), "n_subjects" = 150)
    ),
    pheno             = NULL,
    betas             = c(age = -0.1, sex = 0.05, time = -0.2),
    sigma2            = 0.25,
    subj_intercept_sd = 0.2,
    subj_slope_sd     = 0.1,
    site_sd           = 0.1,
    vw_mean           = 2.5,
    vw_sd             = 0.5,
    fs_template       = "fsaverage",
    roi_subset        = c("temporalpole", "frontalpole", "entorhinal"),
    location_association = NULL,
    measure           = "thickness",
    hemi              = "lh",
    fwhmc             = "fwhm10",
    seed              = 3108,
    verbose           = TRUE) {

  hemi_name <- if (hemi == "lh") "left" else "right"
  vw_message(
    "* creating FreeSurfer dataset (", hemi_name, " hemisphere)...",
    verbose = verbose
  )
  set.seed(seed)

  check_path(path, create_if_not = TRUE)

  mgh_fname <- paste(hemi, measure, fwhmc, fs_template, "mgh", sep = ".")
  n_verts   <- count_vertices(fs_template)

  # ── Locate active vertices from annotation files ──────────────────────────
  roi_locs <- locate_roi(
    rois    = roi_subset,
    n_verts = n_verts,
    hemi    = hemi,
    verbose = verbose
  )
  V_roi <- sum(roi_locs)

  # Vertices that carry the fixed-effect signal (may be a strict subset)
  signal_locs <- if (!is.null(location_association)) {
    sig <- locate_roi(
      rois    = location_association,
      n_verts = n_verts,
      hemi    = hemi,
      verbose = FALSE
    )
    sig & roi_locs                      # signal must lie within roi_subset
  } else {
    roi_locs
  }
  # Boolean mask of length V_roi: TRUE where betas are injected
  signal_within_roi <- signal_locs[roi_locs]

  # ── Ground-truth parameters drawn once (shared across all subjects/sites) ──
  beta0 <- stats::rnorm(V_roi, mean = vw_mean, sd = vw_sd)   # μ_v

  site_names <- names(data_structure)
  K          <- length(site_names)
  u_site     <- if (K > 1L) stats::rnorm(K, 0, site_sd) else 0  # u_k

  # ── Per-site loop ─────────────────────────────────────────────────────────
  for (k in seq_len(K)) {

    site <- site_names[k]
    n    <- data_structure[[site]]$n_subjects
    sess <- data_structure[[site]]$sessions

    vw_message(
      "* ", site, ": ", n, " (subjects) x ", length(sess),
      " (sessions) cortical ", measure, " files",
      verbose = verbose
    )

    site_dir   <- file.path(path, site)
    dir.create(site_dir, showWarnings = FALSE)

    pheno_site <- if (!is.null(pheno)) pheno[pheno$site == site, ] else NULL

    # Subject random effects — drawn once per subject, applied across sessions
    ri_subj <- stats::rnorm(n, 0, subj_intercept_sd)           # b_i
    rs_subj <- if (subj_slope_sd > 0) {
      stats::rnorm(n, 0, subj_slope_sd)                         # b_i^time
    } else {
      rep(0, n)
    }

    for (i in seq_len(n)) {
      for (t_idx in seq_along(sess)) {

        t          <- t_idx - 1L        # 0-based time index, matches pheno$time
        sess_label <- sess[t_idx]

        sub_dir <- file.path(site_dir, paste0("sub-", i, "_ses-", sess_label))
        dir.create(file.path(sub_dir, "surf"), recursive = TRUE,
                   showWarnings = FALSE)

        # ── Fixed-effect contribution: Σ β_j * x_ijt (zero outside signal) ──
        fe <- numeric(V_roi)

        if (!is.null(pheno_site)) {
          row <- pheno_site[pheno_site$id == i & pheno_site$time == t, ]

          for (nm in names(betas)) {
            if (nm == "time") {
              # time index may not be a pheno column; use t directly
              fe[signal_within_roi] <- fe[signal_within_roi] + betas[["time"]] * t
            } else if (nm %in% names(row)) {
              fe[signal_within_roi] <- fe[signal_within_roi] +
                betas[[nm]] * row[[nm]]
            }
          }
        }

        # ── Assemble vertex-wise surface vector ───────────────────────────
        eps <- stats::rnorm(V_roi, 0, sqrt(sigma2))

        vw_vec           <- numeric(n_verts)           # full surface, zeros outside ROI
        vw_vec[roi_locs] <- beta0 +                    # μ_v  (vertex intercept)
                            fe +                       # Σ β_j x_ijt (fixed effects)
                            ri_subj[i] + rs_subj[i] * t +   # b_i + b_i^time * t
                            u_site[k] +                # u_k  (site intercept)
                            eps                        # ε_itkv

        vw_vec[vw_vec < 0] <- 0.01      # surface metrics must be non-negative

        save.mgh(as.mgh(vw_vec), file.path(sub_dir, "surf", mgh_fname))
      }
    }
  }

  invisible(list(
    beta0             = beta0,
    betas             = betas,
    sigma2            = sigma2,
    subj_intercept_sd = subj_intercept_sd,
    subj_slope_sd     = subj_slope_sd,
    site_sd           = site_sd,
    u_site            = u_site,
    signal_vertices   = signal_locs
  ))
}
