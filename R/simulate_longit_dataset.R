# =============================================================================
# Longitudinal simulation helpers
# =============================================================================

#' @title Build a minimal phenotype data frame from a data structure specification
#' @description
#' Creates a long-format data frame with one row per subject × session
#' combination, containing bookkeeping columns used by the simulation
#' pipeline. No phenotypic variables are added.
#'
#' @param data_structure A named list where each element corresponds to one
#'   site/cohort. Every element must contain:
#'   \describe{
#'     \item{\code{sessions}}{Character vector of unique session labels
#'       (e.g. \code{c("01", "02")}).}
#'     \item{\code{n_subjects}}{Positive integer giving the number of
#'       subjects at this site.}
#'   }
#'
#' @return A \code{data.frame} in long format with columns:
#'   \describe{
#'     \item{\code{id}}{Integer subject identifier (within-site).}
#'     \item{\code{site}}{Site/cohort label (character).}
#'     \item{\code{time}}{Session label (character).}
#'     \item{\code{folder_id}}{FreeSurfer-style path fragment,
#'       e.g. \code{"cohort1/sub-1_ses-01"}.}
#'   }
#'
#' @keywords internal
build_minimal_pheno <- function(data_structure) {
  stopifnot(is.list(data_structure), !is.null(names(data_structure)))

  site_dfs <- lapply(names(data_structure), function(site) {
    spec <- data_structure[[site]]

    if (!is.list(spec) || is.null(spec$sessions) || is.null(spec$n_subjects)) {
      vw_error("Each element of data_structure must contain 'sessions' and 'n_subjects'.")
    }

    sess <- as.character(spec$sessions)
    n <- as.integer(spec$n_subjects)

    if (length(sess) < 1L || anyNA(sess) || anyDuplicated(sess)) {
      vw_error("Sessions must be a non-empty character vector with unique labels.")
    }
    if (n <= 0L) {
      vw_error("n_subjects must be a positive integer for site '{site}'")
    }

    min_pheno <- expand.grid(id = seq_len(n), site = site, time = sess, stringsAsFactors = FALSE)

    min_pheno$folder_id <- file.path(site, paste0("sub-", min_pheno$id, "_ses-", min_pheno$time))

    min_pheno
  })

  do.call(rbind, site_dfs)
}

#' @title Validate a user-supplied phenotype data frame
#' @description
#' Checks that `pheno` is a `data.frame` containing all columns required 
#' by `roi_associations` + a `folder_id` column. If `pheno = NULL` and 
#' `roi_associations` is empty, a minimal bookkeeping-only phenotype is built
#' from `data_structure`.
#'
#' @param pheno A `data.frame` in long format, or `NULL`.
#' @param data_structure Named list passed to [build_minimal_pheno()]
#'   when `pheno = NULL`.
#' @param roi_associations Named list of numeric vectors as accepted by
#'   [simulate_freesurfer_data()].
#'
#' @return The validated (or freshly built) `data.frame`.
#'
#' @keywords internal
validate_pheno <- function(pheno, data_structure, roi_associations) {
  
  if (is.null(pheno)) {
    if (length(roi_associations) > 0) vw_error('If associations are specified a matching phenotype must be provided')
    return(build_minimal_pheno(data_structure))
  }

  if (!is.data.frame(pheno)) vw_error("'pheno' must be a data.frame in long format.")

  required_vars <- c(unique(unlist(lapply(roi_associations, names))), 'folder_id')
  missing_cols  <- setdiff(required_vars, names(pheno))
  if (length(missing_cols) > 0L) vw_error("'pheno' is missing required column(s): {missing_cols}")
  
  pheno
}

#' @title Validate ROI association specifications
#' @description
#' Checks that `roi_associations` is a named list of named numeric vectors whose names
#' correspond to known Desikan-Killiany ROI labels.
#'
#' @param roi_associations Named list of named numeric vectors specifying fixed-effect 
#'   beta coefficients per ROI.  May be `NULL` or an empty `list()`.
#' @param simulate_other_rois Logical. If `FALSE`, an empty `roi_associations` is an 
#'   error because no data would be simulated.
#'
#' @return The validated `roi_associations` list (invisibly).
#'
#' @keywords internal
validate_roi_associations <- function(roi_associations, simulate_other_rois) {

  if (is.null(roi_associations) || length(roi_associations) == 0L) {
    if (!simulate_other_rois) vw_error("simulate_other_rois is FALSE but roi_associations is empty! No data is simulated.") 
    return(list())
  }
  if (!is.list(roi_associations) || is.null(names(roi_associations))) vw_error("'roi_associations' must be a named list of named numeric vectors.")
  
  if (anyDuplicated(names(roi_associations))) {
    vw_error("'roi_associations' contains duplicated ROI names. That is not allowed.")
  }

  dk_rois <- unique(locate_roi()$roi_label)
  dk_rois <- dk_rois[!is.na(dk_rois)]

  unknown_rois <- setdiff(names(roi_associations), dk_rois)
  if (length(unknown_rois) > 0L) {
    vw_error("Unknown ROI(s) in 'roi_associations': {unknown_rois}. Please use one of: {dk_rois}")
  }

  lapply(roi_associations, function(x) {
    stopifnot(is.numeric(x) && !is.null(names(x)))
  })
  
  roi_associations
}

# =============================================================================
#  simulate_long_pheno_data()
# =============================================================================

#' @title Simulate longitudinal phenotype data
#' @description
#' Generates a long-format phenotype `data.frame` for one or more cohorts/sites 
#' across multiple sessions / time points.
#'
#' @section Supported covariates:
#' Covariates are declared through the `baseline` argument. The type of covariate
#' is inferred from the names of each list element:
#' \describe{
#'   \item{Continuous (e.g. `age`, `wisdom`)}{Specify `c(mean = m, sd = s)`.
#'      Baseline values are drawn from \eqn{N(\var{m}, \var{s}^2)}.}
#'   \item{Categorical (e.g. `sex`)}{Specify `c(levels = c("Male", "Female"))`.
#'     Each subject is assigned a level with equal probability and this value 
#'     is held constant across sessions.}
#' }
#'
#' @section Change model:
#' Follow-up values for continuous variables listed in `change` are computed as:
#' \deqn{y_{i,s} = y_{i,1} + (s-1)\bar{\delta} + \varepsilon, \quad \varepsilon \sim N(0, \sigma_\delta^2)}
#' where \eqn{s} is the 1-based session index.  The deviation from baseline is
#' drawn independently at each session.
#' The variance \eqn{\sigma_\delta^2} is constant across sessions; only the mean 
#' shift accumulates linearly.
#'
#' @param data_structure Named list specifying site/cohort structure.  Each
#'   element must be a list with:
#'   \describe{
#'     \item{`sessions`}{Character vector of unique session labels.}
#'     \item{`n_subjects`}{Positive integer, number of subjects.}
#'   }
#' @param baseline Named list defining baseline distributions.  Supports continuous 
#'   covariates (`c(mean, sd)`) and categorical covariates (`c(levels = ...)`).
#'   See section **Supported covariates**.
#' @param change Named list of `c(mean, sd)` specifying the per-wave mean shift and 
#'   noise SD for longitudinal covariates.  Only continuous variables should appear 
#'   here.
#' @param seed Integer random seed for reproducibility.
#' @param verbose Logical
#'
#' @return A `data.frame` in long format with one row per subject × session.
#' Always contains columns `site`, `id`, `time`, `folder_id`, + all covariates 
#' declared in `baseline`.
#'
#' @seealso [simulate_freesurfer_data()], [simulate_longit_dataset()]
#'
#' @examples
#' pheno <- simulate_long_pheno_data(
#'   data_structure = list(
#'     GENR = list(sessions = c("01", "02"), n_subjects = 50)
#'   ),
#'   baseline = list(age = c(mean = 10, sd = 1), sex = c(levels = c("Male", "Female"))),
#'   change   = list(age = c(mean = 4, sd = 0.5))
#' )
#' head(pheno)
#'
#' @export
#' 
simulate_long_pheno_data <- function(
    data_structure = list(
      cohort1 = list(sessions = c("01", "02"), n_subjects = 100),
      cohort2 = list(sessions = c("01", "02"), n_subjects = 150)
    ),
    baseline = list(age = c(mean=10, sd=0.5), 
                    sex = c(levels=c('Male','Female')), 
                    wisdom = c(mean = 0, sd = 1)),
    change = list(age = c(mean=4, sd=0.5), 
                  wisdom = c(mean = 1, sd = 0.5)),
    seed           = 3108,
    verbose        = TRUE) {

  if (verbose) cli::cli_progress_step('Generate phenotype data', spinner=TRUE)
  set.seed(seed)

  site_frames <- lapply(names(data_structure), function(site) {
    spec <- data_structure[[site]]
    n <- as.integer(spec$n_subjects)
    sess <- spec$sessions

    baseline_df <- data.frame(id = seq_len(n), site = site, time = sess[1])

    for (var in names(baseline)) {
      var_def <- baseline[[var]]
      if (identical(sort(names(var_def)), c("mean", "sd"))) {
        baseline_df[[var]] <- stats::rnorm(n, mean = var_def['mean'], sd = var_def['sd'])
      } else {
        raw <- sample(unname(var_def), n, replace = TRUE)
        baseline_df[[var]] <- factor(raw, levels = unname(var_def))
      }
    }

    long_df <- if (length(sess) == 1) { baseline_df 
    } else {
      follow_dfs <- lapply(seq(2, length(sess)), function(s) {
        follow_df <- baseline_df
        follow_df$time <- sess[s]

        for (var in names(change)) {
          var_def <- change[[var]]
          follow_df[[var]] <- baseline_df[[var]] + stats::rnorm(n, mean = (s-1L)*var_def['mean'], sd = var_def['sd'])
        }

        follow_df
      })

      do.call(rbind, c(list(baseline_df), follow_dfs))
    }
  })
  
  pheno <- do.call(rbind, site_frames)

  pheno$folder_id = file.path(pheno$site, paste0("sub-", pheno$id, "_ses-", pheno$time))

  vw_message(c("i" = "{.val {nrow(pheno)}} total observations (from {.val2 {length(data_structure)}} site{?s} 
                       and {.val2 {length(unique(pheno$time))}} wave{?s})", 
               " " = "Variables: {.strong {names(pheno)}}"), verbose = verbose)

  pheno
}

# =============================================================================
#  simulate_freesurfer_data()
# =============================================================================

#' @title Simulate longitudinal FreeSurfer vertex-wise surface data
#' @description
#' Writes one `.mgh` file per observation (subject × session) to `path`, mimicking
#' the FreeSurfer recon-all output layout. Vertex values are generated from a 
#' mixed-effects model with:
#' \itemize{
#'   \item ROI-specific fixed effects specified via `roi_associations`;
#'   \item a scalar subject random intercept (shared across all vertices);
#'   \item an optional scalar site random intercept.
#' }
#'
#' @section ROI targeting:
#' \describe{
#'   \item{`simulate_other_rois = FALSE`}{Only ROIs in `roi_associations` are filled;
#'      all other vertices remain 0.}
#'   \item{`simulate_other_rois = TRUE`}{All Desikan-Killiany cortical ROIs are filled.
#'     ROIs absent from `roi_associations` are simulated under the null (no fixed effects,
#'     noise only).}
#' }
#'
#' @section Data-generating process:
#' For observation \eqn{i} and vertex \eqn{v} in ROI \eqn{r}: 
#' \deqn{y_{iv} = \mathbf{x}_i^\top \boldsymbol{\beta}_r + b_i + b_{\text{site}(i)} + \varepsilon_{iv}}
#' where \eqn{b_i \sim N(0, \sigma_\text{subj}^2)},
#' \eqn{b_s \sim N(0, \sigma_\text{site}^2)}, and
#' \eqn{\varepsilon_{iv} \sim N(\mu_\text{vw}, \sigma_\text{vw}^2)}.
#' Random effects are scalar (not vertex-specific).
#'
#' @param path Character string; root directory in which to write `.mgh` files. 
#'   Created if it does not exist.
#' @param pheno `data.frame` in long format from [simulate_long_pheno_data()], 
#'   or `NULL` to auto-generate a minimal bookkeeping frame from `data_structure`.
#' @param data_structure Named list; see [simulate_long_pheno_data()].
#'   Ignored when `pheno` is supplied.
#' @param roi_associations Named list of named numeric vectors. Names of the list
#'   are Desikan-Killiany ROI labels; names of each vector are column names in 
#'   `pheno`; values are beta coefficients. 
#'   Example: `list(temporalpole = c(age = 0.3, wisdom = 0.5))`.
#' @param simulate_other_rois Logical; if `TRUE` all non-association ROIs are 
#'   simulated under the null.
#' @param hemi One of `"lh"` or `"rh"`.
#' @param measure Character; FreeSurfer surface measure, e.g. `"thickness"` or
#'   `"area"`.
#' @param vw_mean Numeric; mean of the vertex-level residual noise (i.e. 
#'   grand-mean cortical thickness).
#' @param vw_sd Numeric \eqn{\geq 0}; SD of vertex-level residual noise.
#' @param subj_sd Numeric \eqn{\geq 0}; SD of the subject random intercept.
#' @param site_sd Numeric \eqn{\geq 0}; SD of the site random intercept.
#'   Ignored when there is only one site (intercept drawn as a single value).
#' @param fs_template Character; FreeSurfer template space, e.g. `"fsaverage"`.
#' @param fwhmc Character; smoothing kernel label, e.g. `"fwhm10"`.
#' @param seed Integer random seed.
#' @param verbose Logical; progress messages.
#'
#' @return `NULL` invisibly. Side effect: `.mgh` files written to
#'   `path/<folder_id>/surf/<hemi>.<measure>.<fwhmc>.<fs_template>.mgh`.
#' @seealso [simulate_long_pheno_data()], [simulate_longit_dataset()]
#'
#' @examplesIf dir.exists('path/to/simulated_brains')
#' simulate_freesurfer_data(
#'   path = 'path/to/simulated_brains',
#'   pheno = pheno,
#'   roi_associations = list(temporalpole = c(age = 0.3)),
#'   hemi = "both"
#' )
#'
#' @importFrom stats rnorm setNames
#' @export
#' 
simulate_freesurfer_data <- function(
    path,
    pheno = NULL,
    data_structure = list(
      cohort1 = list(sessions = c("01", "02"), n_subjects = 100),
      cohort2 = list(sessions = c("01", "02"), n_subjects = 150)
    ),
    roi_associations = list(),
    simulate_other_rois = FALSE,
    hemi = "lh",
    measure = "thickness",
    vw_mean = 2.5,
    vw_sd = 0.5,
    subj_sd = 0.2,
    site_sd = 0.1,
    fs_template = "fsaverage",
    fwhmc = "fwhm10",
    seed = 3108,
    verbose = TRUE) {

  hemi_name <- if (hemi == "lh") "left" else "right"

  if (verbose) cli::cli_progress_step('Generate FreeSurfer dataset ({hemi_name} hemisphere)', spinner=TRUE)
 
  stopifnot(
    is.numeric(vw_sd), vw_sd >= 0,
    is.numeric(subj_sd), subj_sd >= 0,
    is.numeric(site_sd), site_sd >= 0
  )

  set.seed(seed)

  check_path(path, create_if_not = TRUE)

  roi_associations <- validate_roi_associations(roi_associations = roi_associations, simulate_other_rois = simulate_other_rois)
  pheno <- validate_pheno(pheno = pheno, data_structure = data_structure, roi_associations = roi_associations)
  
  n_verts <- count_vertices(fs_template)
  n_obs <- nrow(pheno)

  ss <- matrix(0, nrow = n_obs, ncol = n_verts, dimnames = list(pheno$folder_id, NULL))

  # Build random terms 
    
  subj_key <- paste(pheno$site, pheno$id, sep = ":")  # e.g. "GENR:1"
  unique_subjs <- unique(subj_key)
  subj_re_map <- setNames(rnorm(length(unique_subjs), 0, subj_sd), unique_subjs)
  ri_subj <- subj_re_map[subj_key] # intercept: same subject across timepoints gets same RE
  
  unique_sites <- unique(pheno$site)
  if (length(unique_sites) == 1) { ri_site = 0L } else {
    site_re_map <- setNames(rnorm(length(unique_sites), 0, site_sd), unique_sites)
    ri_site <- site_re_map[pheno$site] 
  }

  roi_map <- aparc.annot[[hemi]]$label_names[1:n_verts]

  if (simulate_other_rois) {
    roi_locs <- which(!roi_map %in% names(roi_associations))
    ss[, roi_locs] <- ri_site + ri_subj + rnorm(n_obs * length(roi_locs), mean = vw_mean, sd = vw_sd)
    ss[, roi_locs] <- pmax(ss[, roi_locs], 0.001)
  }

  if (length(roi_associations) > 0) {
    for (roi in names(roi_associations)) {
      roi_spec <- roi_associations[[roi]]
      fixed_part <- as.vector(as.matrix(pheno[names(roi_spec)]) %*% roi_spec)
      roi_locs <- which(roi_map == roi)
      ss[, roi_locs] <- fixed_part + ri_site + ri_subj + rnorm(n_obs * length(roi_locs), mean = vw_mean, sd = vw_sd)
      ss[, roi_locs] <- pmax(ss[, roi_locs], 0.001)
    }
  }

  mgh_fname <- paste(hemi, measure, fwhmc, fs_template, "mgh", sep = ".")

  for (obs in pheno$folder_id) {
    obs_dir <- file.path(path, obs, "surf")
    dir.create(obs_dir, recursive = TRUE, showWarnings = FALSE)
    save.mgh(as.mgh(ss[obs, ]), file.path(obs_dir, mgh_fname))
  }

  invisible(NULL)
}

# =============================================================================
#  simulate_longit_dataset()
# =============================================================================

#' @title Simulate a complete longitudinal brain surface dataset
#' @description
#' High-level wrapper that:
#' \enumerate{
#'   \item Generates a longitudinal phenotype `data.frame` via
#'         [simulate_long_pheno_data()];
#'   \item Writes it to `path/phenotype.csv`;
#'   \item Calls [simulate_freesurfer_data()] for each requested hemisphere.
#' }
#'
#' @inheritParams simulate_long_pheno_data
#' @inheritParams simulate_freesurfer_data
#' @param hemi One of `"lh"`, `"rh"`, or `"both"` (default: `"both"`).
#'
#' @return `NULL` invisibly. Side effects: `phenotype.csv` and vertex-wise
#'   `.mgh` files written under `path`.
#'
#' @seealso [simulate_long_pheno_data()], [simulate_freesurfer_data()]
#'
#' @examplesIf dir.exists('path/to/simulated/dataset')
#' simulate_longit_dataset(
#'   path = 'path/to/simulated/dataset/',
#'   data_structure = list(
#'     GENR = list(sessions = c("01", "02", "03"), n_subjects = 50)
#'   ),
#'   roi_associations = list(temporalpole = c(age = 1.3)),
#'   hemi = "lh"
#' )
#'
#' @export
#' 
simulate_longit_dataset <- function(
    path,
    data_structure = list(
      cohort1 = list(sessions = c("01", "02", "03"), n_subjects = 100),
      cohort2 = list(sessions = c("01", "02"), n_subjects = 150)
    ),
    baseline = list(age = c(mean=10, sd=0.5), 
                    sex = c(levels=c('Male','Female')), 
                    wisdom = c(mean = 0, sd = 1)),
    change = list(age = c(mean=4, sd=0.5), 
                  wisdom = c(mean = 1, sd = 0.5)),
    roi_associations = list(temporalpole = c(age = 1.3, sex = 0.5), 
                            entorhinal = c(age = 0.9), 
                            frontalpole = c(wisdom = 0.7)),
    simulate_other_rois = FALSE,
    hemi = "both",
    measure = "thickness",
    vw_mean = 2.5,
    vw_sd = 0.5,
    subj_sd = 0.2,
    site_sd = 0.1,
    fs_template = "fsaverage",
    fwhmc = "fwhm10",
    seed = 3108,
    verbose = TRUE) {

  vw_init_message("Simulating longitudinal dataset", verbose = verbose)

  check_path(path, create_if_not = TRUE)

  pheno <- simulate_long_pheno_data(
    data_structure = data_structure,
    baseline = baseline,
    change = change,
    seed = seed,
    verbose = verbose)
  
  utils::write.csv(pheno, file.path(path,'phenotype.csv'), row.names = FALSE)
  
  hemis <- if (hemi == 'both') c('lh', 'rh') else c(hemi) 
  
  lapply(hemis, function(h) {
    simulate_freesurfer_data(path = path,
                              pheno = pheno,
                              roi_associations = roi_associations,
                              simulate_other_rois = simulate_other_rois,
                              hemi = h,
                              measure = measure,
                              vw_mean = vw_mean,
                              vw_sd = vw_sd,
                              subj_sd = subj_sd,
                              site_sd = site_sd,
                              fs_template = fs_template,
                              fwhmc = fwhmc,
                              seed = seed,
                              verbose = verbose)
  })
  
  vw_message("Done! :)", type = "step", verbose = verbose)
  invisible(NULL)
}