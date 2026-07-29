#' @title Plot the vertex-wise difference between two scalar maps
#'
#' @description
#' Computes the vertex-wise difference (`a - b`) between two surface maps
#' (numeric vectors or MGH/GII file paths) for each hemisphere, then renders
#' the result on a standard fsaverage surface via [plot_vw_surf()].
#'
#' Useful for contrasting two conditions, time-points, groups, or model terms
#' without pre-computing the difference outside R.
#'
#' @param lh_a,lh_b Left-hemisphere map: numeric vector or path to an MGH/GII
#'   file. Both must have the same length / vertex count. Pass `NULL` to
#'   omit the left hemisphere entirely.
#' @param rh_a,rh_b Right-hemisphere map: same as above for the right
#'   hemisphere. Pass `NULL` to omit the right hemisphere.
#' @param label_a,label_b Short character labels used in the default figure
#'   title (e.g. `"group A"`, `"group B"`). Ignored when `title` is supplied via 
#'   `...`.
#' @param ... Additional arguments forwarded to [plot_vw_surf()]
#'   (e.g.`surface`,`views`,`cmap`,`vmin`,`vmax`,`threshold`,
#'   `colorbar`,`colorbar_label`,`title`, `to_file`,`dpi`,`fs_home`,`fs_template`).
#' 
#' @param ... Additional arguments forwarded to [plot_vw_surf()],
#'   e.g. `roi_outline`, `views`, `cmap`, `vmin`, `vmax`, `colorbar`, `colorbar_label`, `title`, `to_file`,
#'   `dpi`, `fs_home`, `fs_template.`
#'
#' @return Invisibly: the output of [plot_vw_surf()] — the temp HTML file path (interactive 
#'   mode) or `to_file` path (static PNG mode).
#'   Called primarily for its side-effect of opening or saving the figure.
#'
#' @seealso [plot_vw_surf()], [plot_vw_map()] [vw_diff()]
#'
#'
#' @examplesIf rlang::is_installed("reticulate") && dir.exists("path/to/model1/")
#' # Two numeric vectors
#' plot_vw_diff(
#'   lh_a = lh_term1, lh_b = lh_term2,
#'   rh_a = rh_term1, rh_b = rh_term2,
#'   label_a = "Term 1", label_b = "Term 2",
#'   cmap = "RdBu_r"
#' )
#'
#' # Two MGH files, left hemisphere only
#' plot_vw_diff(
#'   lh_a = "path/to/model1/lh.area.aic.mgh",
#'   lh_b = "path/to/model2/lh.area.aic.mgh",
#'   label_a = "model1", label_b = "model2"
#' )
#'
#' @export
#' 
plot_vw_diff <- function(lh_a = NULL, lh_b = NULL,
                         rh_a = NULL, rh_b = NULL,
                         label_a = "a", label_b = "b",
                         ...) {
  
  require_packages('reticulate', call_fn = 'plot_vw_diff')

  # helper: load from file path or return numeric vector as-is
  diff_map <- vw_diff(lh_a = lh_a, lh_b = lh_b, rh_a = rh_a, rh_b = rh_b,
                      label_a = label_a, label_b = label_b)

  dots <- list(...)

  dots$title  <- dots$title %||% paste0("Difference: ", label_a, " \u2212 ", label_b)

  do.call(plot_vw_surf, c(
    list(lh = diff_map[['lh']], rh = diff_map[['rh']]),
    dots
  ))

}