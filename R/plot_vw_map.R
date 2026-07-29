#' @title Plot vertex-wise coefficient maps on a 3D cortical surface
#' 
#' @description
#' Locates the vertex-wise coefficient MGH files for a given model term and
#' surface measure, optionally applies cluster-wise significance (CWS) masking,
#' and renders the result on a standard fsaverage surface via
#' [plot_vw_surf()].
#'
#' Surface maps are loaded directly in R using [load.mgh()] without any
#' Python dependency at the file-loading stage. The `Python/nilearn` rendering
#' backend is invoked only through [plot_vw_surf()].
#'
#' @param res_dir Character. Path to the directory containing FreeSurfer-style
#'   vertex-wise result files (`*.mgh`) and `stack_names.txt`.
#' @param term Character. Name of the model term to visualize (matched against
#'   entries in `stack_names.txt`).
#' @param measure Character. Surface measure to load, e.g. `"area"`, `"thickness"`,
#'   `"volume"`. Default `"area"`.
#' @param hemi Character. Which hemisphere(s) to plot: `"both"` (default), `"lh"` or `"rh"`.
#' @param surface Character. Surface mesh: `"pial"` (default) or `"inflated"`.
#' @param threshold Controls vertex-level masking before plotting:
#'   \describe{
#'     \item{`"cws"` (default)}{Cluster-wise significant. Loads the matching 
#'       `*.cache.*.sig.ocn.mgh` file and masks all vertices not belonging to a significant 
#'       cluster. If no OCN file is found, the unmasked coefficients are plotted with a warning.}
#'    \item{`"fdr <= 0.05"` or `"p < 0.01"` etc}{FDR or P values below a threshold. Loads the matching 
#'       `*.fdr.mgh` or `*.p.mgh` file and masks all vertices not meeting the rule. 
#'       If no file is found, the unmasked coefficients are plotted with a warning.}
#'     \item{Numeric}{Absolute-value threshold (vertices with `|value| < threshold` are hidden).}
#'     \item{`NULL`}{No masking; all vertices are rendered.}
#'   }
#' @param ... Additional arguments forwarded to [plot_vw_surf()],
#'   e.g. `roi_outline`, `views`, `cmap`, `vmin`, `vmax`, `colorbar`, `colorbar_label`, `title`, `to_file`,
#'   `dpi`, `fs_home`, `fs_template.`
#'
#' @return Invisibly: the output of [plot_vw_surf()] — the temp HTML file path (interactive 
#'   mode) or `to_file` path (static PNG mode).
#'   Called primarily for its side-effect of opening or saving the figure.
#'
#' @seealso [plot_vw_surf()], [plot_vw_diff()]
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
  
  require_packages('reticulate', call_fn = 'plot_vw_map')

  # Match mesh arguments
  hemi <- match.arg(hemi)
  surface <- match.arg(surface)

  # Validate paths
  if (!dir.exists(res_dir))
    vw_error("Results directory does not exist: {.file {res_dir}}")

  stack_file <- file.path(res_dir, "stack_names.txt")
  if (!file.exists(stack_file)) {
    stack = term # assume this is a meta-analysis? 
    # vw_error(c(
    #   "Cannot find {.file stack_names.txt} in {.file {res_dir}}.",
    #   " " = "Results folder may be incorrect or corrupted."
    # ))
    vw_message(c(
      "i" = "Cannot find {.file stack_names.txt} in {.file {res_dir}}.",
      " " = "Assuming this is a meta-analysis."
    ))

  } else {
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
  }
  
  hemis_to_load <- if (hemi == "both") c("lh", "rh") else hemi

  lh_data <- if ("lh" %in% hemis_to_load) load_and_mask_coef("lh", measure, stack, res_dir, threshold)
  rh_data <- if ("rh" %in% hemis_to_load) load_and_mask_coef("rh", measure, stack, res_dir, threshold)

  if (is.null(lh_data[['coef']]) && is.null(rh_data[['coef']])) {
    vw_error(c(
      "No coefficient MGH files found for term {.val {term}} / measure {.val {measure}}.",
      "i" = "Expected e.g. {.file lh.{measure}.{stack}.coef.mgh} in {.file {res_dir}}"))
  }
  
  dots  <- list(...)
  title <- dots$title %||% paste("Effect of", term, "on", measure)
  # remove title from dots to avoid duplicate argument
  dots$title <- NULL

  do.call(plot_vw_surf, c(
    list(lh = lh_data[['coef']], rh = rh_data[['coef']], lh_mask = lh_data[['mask']], rh_mask = rh_data[['mask']], 
         surface = surface, title = title),
         dots))
}