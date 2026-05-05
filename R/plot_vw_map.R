#' @title Plot vertex-wise coefficient maps on a 3D cortical surface
#'
#' @description
#' Locates the vertex-wise coefficient MGH files for a given model term and
#' surface measure, optionally applies cluster-wise significance (CWS) masking,
#' and renders the result on a standard fsaverage surface via
#' \code{\link{plot_vw_surf}}.
#'
#' Surface maps are loaded directly in R using \code{load.mgh()} without any
#' Python dependency at the file-loading stage. The Python/nilearn rendering
#' backend is invoked only through \code{\link{plot_vw_surf}}.
#'
#' @param res_dir Character. Path to the directory containing FreeSurfer-style
#'   vertex-wise result files (\code{*.mgh}) and \code{stack_names.txt}.
#' @param term Character. Name of the model term to visualize (matched against
#'   entries in \code{stack_names.txt}).
#' @param measure Character. Surface measure to load, e.g. \code{"area"},
#'   \code{"thickness"}, \code{"volume"}. Default \code{"area"}.
#' @param hemi Character. Which hemisphere(s) to plot: \code{"both"} (default),
#'   \code{"lh"}, or \code{"rh"}.
#' @param surface Character. Surface mesh: \code{"pial"} (default) or
#'   \code{"inflated"}.
#' @param threshold Controls vertex-level masking before plotting:
#'   \describe{
#'     \item{\code{"cws"} (default)}{Cluster-wise significance masking. Loads
#'       the matching \code{*.cache.*.sig.ocn.mgh} file and sets all vertices
#'       not belonging to a significant cluster (OCN label \code{== 0}) to
#'       \code{NA}. If no OCN file is found, the unmasked coefficients are
#'       plotted with a warning.}
#'     \item{Numeric}{Passed directly to \code{\link{plot_vw_surf}} as an
#'       absolute-value threshold (vertices with \code{|value| < threshold}
#'       are hidden).}
#'     \item{\code{NULL}}{No masking; all vertices are rendered.}
#'   }
#' @param ... Additional arguments forwarded to \code{\link{plot_vw_surf}},
#'   e.g. \code{views}, \code{cmap}, \code{vmin}, \code{vmax},
#'   \code{colorbar}, \code{colorbar_label}, \code{title}, \code{to_file},
#'   \code{dpi}, \code{fs_home}, \code{fs_template}.
#'
#' @return Invisibly: the output of \code{\link{plot_vw_surf}} — the temp HTML
#'   file path (interactive mode) or \code{to_file} path (static PNG mode).
#'   Called primarily for its side-effect of opening or saving the figure.
#'
#' @seealso \code{\link{plot_vw_surf}}, \code{\link{plot_vw_diff}}
#'
#' @examplesIf rlang::is_installed("reticulate") && dir.exists("~/results/fs_results")
#' # Both hemispheres, interctive plot, cluster-wise masking (default)
#' plot_vw_map(
#'   res_dir = "~/results/fs_results",
#'   term    = "age",
#'   measure = "thickness"
#' )
#'
#' # Left hemisphere only, numeric threshold, save to PNG
#' plot_vw_map(
#'   res_dir   = "~/results/fs_results",
#'   term      = "age",
#'   hemi      = "lh",
#'   threshold = 0.05,
#'   to_file   = "figures/age_lh.png",
#'   dpi       = 300L
#' )
#'
#' # No masking, custom colour map and title
#' plot_vw_map(
#'   res_dir   = "~/results/fs_results",
#'   term      = "sex",
#'   threshold = NULL,
#'   cmap      = "RdBu_r",
#'   title     = "Sex difference in cortical area"
#' )
#'
#' @export
#' 
plot_vw_map <- function(res_dir, term,  measure = 'area', 
                        hemi = c("both", "lh", "rh"),
                        surface = c('pial','inflated'),
                        threshold = 'cws',
                        ...
                        # outline_rois = NULL
                      ) {

  # Match mesh arguments
  hemi <- match.arg(hemi)
  surface <- match.arg(surface)

  # Validate paths
  if (!dir.exists(res_dir))
    vw_error("Results directory does not exist: {.file {res_dir}}")

  stack_file <- file.path(res_dir, "stack_names.txt")
  if (!file.exists(stack_file))
    vw_error(c(
      "Cannot find {.file stack_names.txt} in {.file {res_dir}}.",
      " " = "Results folder may be incorrect or corrupted."
    ))

  # Read existing file
  stack_ids <- utils::read.table(stack_file, header = TRUE, sep = "\t",
                                 stringsAsFactors = FALSE)
    if (!term %in% stack_ids$stack_name)
      vw_error(c(
        "Term {.val {term}} not found in {.file stack_names.txt}.",
        "i" = "Available terms: {.or {.val {stack_ids$stack_name}}}"
      ))
  
  # Extract stack
  stack <- paste0('stack', stack_ids[ stack_ids$stack_name == term, 'stack_number'])

  
  hemis_to_load <- if (hemi == "both") c("lh", "rh") else hemi

  .load_hemi_data <- function(h) {

    coef_file <- file.path(res_dir, paste(h, measure, stack, "coef.mgh", sep = "."))

    if (!file.exists(coef_file)) {
      vw_message(c("!" = "Coefficient file not found for {h}, skipping: {.file {coef_file}}"))
      return(NULL)
    }

    # load coef via nibabel through reticulate
    coef <- load.mgh(coef_file)$x

    if (identical(threshold, "cws")) {
      # ocn_file <- file.path(res_dir,
      #               paste(h, measure, stack, "cache.th30.abs.sig.ocn.mgh", sep = "."))
      ocn_file <- list.files(res_dir,
        pattern = paste0("^", h, "\\.", measure, "\\.", stack, "\\.cache\\..*\\.sig\\.ocn\\.mgh$"),
        full.names = TRUE)
      
      if (length(ocn_file) == 0) {
        vw_message(c("!" = "No OCN file found for {h}.", "i" = "Plotting unmasked coefficients."))
      } else {
        if (length(ocn_file) > 1) {
          vw_message("!" = "Multiple OCN files found for {h}, using: {.file {basename(ocn_file[1])}}")
          ocn_file <- ocn_file[1]
        }

        ocn <- load.mgh(ocn_file)$x
        # keep only vertices belonging to a significant cluster (ocn is a
        # positive integer label; non-significant vertices are 0)
        coef[is.na(ocn)] <- NA
      }

    }
    coef
  }

  lh_data <- if ("lh" %in% hemis_to_load) .load_hemi_data("lh") else NULL
  rh_data <- if ("rh" %in% hemis_to_load) .load_hemi_data("rh") else NULL

  if (is.null(lh_data) && is.null(rh_data)) {
    vw_error(c(
      "No coefficient MGH files found for term {.val {term}} / measure {.val {measure}}.",
      "i" = "Expected e.g. {.file lh.{measure}.{stack}.coef.mgh} in {.file {res_dir}}"))
  }
      
  surf_threshold <- if (identical(threshold, "cws")) NULL else threshold

  dots  <- list(...)
  title <- dots$title %||% paste("Effect of", term, "on", measure)
  # remove title from dots to avoid duplicate argument
  dots$title <- NULL

  do.call(plot_vw_surf, c(
    list(lh = lh_data, rh = rh_data, surface = surface, threshold = surf_threshold, title = title),
         dots))
}