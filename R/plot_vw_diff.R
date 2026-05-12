#' @title Plot the vertex-wise difference between two scalar maps
#'
#' @description
#' Computes the vertex-wise difference (\code{a - b}) between two surface maps
#' (numeric vectors or MGH/GII file paths) for each hemisphere, then renders
#' the result on a standard fsaverage surface via \code{\link{plot_vw_surf}}.
#'
#' Useful for contrasting two conditions, time-points, groups, or model terms
#' without pre-computing the difference outside R.
#'
#' @param lh_a,lh_b Left-hemisphere map: numeric vector or path to an MGH/GII
#'   file. Both must have the same length / vertex count. Pass \code{NULL} to
#'   omit the left hemisphere entirely.
#' @param rh_a,rh_b Right-hemisphere map: same as above for the right
#'   hemisphere. Pass \code{NULL} to omit the right hemisphere.
#' @param label_a,label_b Short character labels used in the default figure
#'   title (e.g. \code{"group A"}, \code{"group B"}). Ignored when \code{title}
#'   is supplied via \code{...}.
#' @param ... Additional arguments forwarded to \code{\link{plot_vw_surf}}
#'   (e.g. \code{surface}, \code{views}, \code{cmap}, \code{vmin}, \code{vmax},
#'   \code{threshold}, \code{colorbar}, \code{colorbar_label}, \code{title},
#'   \code{to_file}, \code{dpi}, \code{fs_home}, \code{fs_template}).
#'
#' @return Invisibly: the output of \code{\link{plot_vw_surf}} (temp HTML path
#'   or \code{to_file} path).
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

  # helper: load from file path or return numeric vector as-is
  .resolve <- function(x, arg_name) {
    if (is.null(x))       return(NULL)
    if (is.character(x))  return(load.mgh(x)$x)
    if (is.numeric(x))    return(as.numeric(x))
    vw_error("{.arg {arg_name}} must be a numeric vector or a file path, not {.cls {class(x)}}.")
  }

  # helper: subtract b from a with sanity checks
  .diff_hemi <- function(a, b, hemi) {
    if (is.null(a) && is.null(b)) return(NULL)

    if (is.null(a))
      vw_error(c(
        "{hemi}_a is NULL but {hemi}_b is provided.",
        "i" = "Supply both {hemi}_a and {hemi}_b, or neither."
      ))
    if (is.null(b))
      vw_error(c(
        "{hemi}_b is NULL but {hemi}_a is provided.",
        "i" = "Supply both {hemi}_a and {hemi}_b, or neither."
      ))
    if (length(a) != length(b))
      vw_error(c(
        "{hemi}_a and {hemi}_b have different lengths ({length(a)} vs {length(b)}).",
        "i" = "Both maps must be on the same surface template."
      ))

    a - b
  }

  lh_diff <- .diff_hemi(.resolve(lh_a, "lh_a"), .resolve(lh_b, "lh_b"), "lh")
  rh_diff <- .diff_hemi(.resolve(rh_a, "rh_a"), .resolve(rh_b, "rh_b"), "rh")

  all_diff <- c(lh_diff, rh_diff)
  qs <- stats::quantile(all_diff, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

  min_diff <- qs[[1]]
  max_diff <- qs[[5]]
  
  vw_message(c(
    "i" = "Difference map summary ({sum(is.finite(all_diff))} vertices)",
    " " = "Min:    {.val2 {round(min_diff, 3)}}",
    " " = "Q1:     {.val2 {round(qs[[2]], 3)}}",
    " " = "Median: {.val2 {round(qs[[3]], 3)}}",
    " " = "Mean:   {.val2 {round(mean(all_diff, na.rm = TRUE), 3)}}",
    " " = "Q3:     {.val2 {round(qs[[4]], 3)}}",
    " " = "Max:    {.val2 {round(max_diff, 3)}}"
  ))

  dots <- list(...)

  dots$title  <- dots$title %||% paste0("Difference: ", label_a, " \u2212 ", label_b)

  dots$vmin <- dots$vmin %||% min_diff
  dots$vmax <- dots$vmax %||% max_diff

  do.call(plot_vw_surf, c(
    list(lh = lh_diff, rh = rh_diff),
    dots
  ))

  return(all_diff)
}