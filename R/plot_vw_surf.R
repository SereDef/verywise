#' @title Plot vertex-wise maps on brain surfaces
#'
#' @description
#' Renders left and/or right hemisphere vertex-wise scalar maps on a standard
#' fsaverage surface using Python/nilearn. Fully headless (no XQuartz, rgl,
#' or display server required).
#'
#' Surface meshes are resolved in order:
#' \enumerate{
#'   \item \env{FREESURFER_HOME}/subjects/\code{fs_template}/surf/ - used
#'         directly if the directory exists (no download needed).
#'   \item nilearn automatic download, cached permanently in
#'         \file{~/nilearn_data/} (download only happens once per template).
#' }
#'
#' \strong{Interactive} (\code{to_file = NULL}): generates a
#' self-contained WebGL HTML file opened in the IDE Viewer or default
#' browser.\cr
#' \strong{Static} (\code{to_file} supplied): saves a tiled PNG via
#' matplotlib.
#'
#' @param lh,rh Numeric vector, `.mgh`/`.gii`` file path, or
#'   `NULL`. At least one of `lh`/`rh`` must be supplied.
#' @param lh_mask,rh_mask Masking boolean maps, or `NULL` (no masking).
#' @param vmin,vmax Numeric colormap limits, or `NULL` for automatic
#'   scaling (minimum and maximum of the thresholded data).
#' @param views Character vector of camera angles - any subset of
#'   `"lateral"`, `"medial"`, `"dorsal"`, `"ventral"`,
#'   `"anterior"`, `"posterior"`. Default: `"all"`.
#' @param surface Surface mesh: `"pial"` (default) or `"inflated"`.
#' @param bg_map Background shading: `"sulc"` (default),
#'   `"curv"`, or `"none"`.
#' @param fs_template Resolution (fsaverage template). Must match the length 
#'   of `lh`/`rh`. Default: `"fsaverage"`.
#' @param fs_home (optional) location of FreeSurfer home for templates.
#' @param colorbar Logical. Draw a shared colour bar? Default `TRUE`.
#' @param colorbar_label Character label for the colour bar axis, or `NULL`.
#' @param colorbar_width Fraction of one brain panel's width reserved for the
#'   colorbar/density strip, or `NULL` for the mode default
#'   (static: 0.40, interactive: 0.20).
#' @param cmap Matplotlib colormap name. 
#' @param roi_outline Character vector of DK/aparc region names to outline
#'   with contour lines (e.g. `c("superiorfrontal", "precentral")`),
#'   or `NULL` for no ROI overlay.
#' @param title Character figure title, or `NULL`.
#' @param to_file Path ending in `.png` for static export, or
#'   `NULL` for interactive HTML mode.
#' @param dpi Integer output resolution for static PNG. Default `150`.
#' @param cell_px Brain panel size in pixels: a single number, a length-2 vector
#'   \code{c(width, height)}, or `NULL` for the renderer's default
#'   (static: 400x440, interactive: 500x350).
#' 
#' @return Invisibly: `to_file` path (static) or the temp HTML
#'   file path (interactive). Called primarily for the side-effect.
#'
#' @examples
#' \dontrun{
#' # interactive (opens in Viewer / browser) with file 
#' plot_vw_surf(lh = 'path/to/lh.coef.mgh', fs_template = "fsaverage5")
#'
#' # static PNG (4 views, 300 dpi...) with vector objects
#' plot_vw_surf(
#'   lh = lh_coef,
#'   rh = rh_coef,
#'   threshold = 0.05,
#'   views = c("lateral", "medial", "dorsal", "ventral"),
#'   cmap = "RdBu_r",
#'   title = "Effect of age on cortical thickness",
#'   to_file = "figures/age_thickness.png",
#'   dpi = 300
#' )
#' }
#'
#' @export
plot_vw_surf <- function(
    lh = NULL,
    rh = NULL,
    lh_mask = NULL,
    rh_mask = NULL,
    vmin = NULL,
    vmax = NULL,
    views = 'all',
    surface = c("pial", "inflated"),
    bg_map = c("sulc", "curv", "none"),
    fs_template = "fsaverage",
    fs_home = NULL,
    colorbar = TRUE,
    colorbar_label = NULL,
    colorbar_width = NULL,
    cmap = NULL,
    roi_outline = NULL,
    title = NULL,
    to_file = NULL,
    dpi = 150L,
    cell_px = NULL) {

  # --- initialise Python renderer (once per session) -----------------------

  .vw_surf_init_py()

  # --- input validation ----------------------------------------------------
  if (is.null(lh) && is.null(rh))
    vw_error("At least one of `lh` or `rh` must be supplied.")
  
  if (!is.null(to_file)) {
    if (!grepl("\\.png$", to_file, ignore.case = TRUE)) {
      vw_error("Output file name must end in '.png'.")
    }
  }
      
  surface <- match.arg(surface)
  bg_map  <- match.arg(bg_map)

  valid_views <- c("lateral", "dorsal", "anterior", 
                   "medial", "ventral", "posterior")
  if (identical(views, 'all')) {
    views <- if (!is.null(to_file)) valid_views else 'lateral'
  } else {
    bad_views <- setdiff(views, valid_views)
    if (length(bad_views))
        vw_error(c("Invalid view{?s}: {bad_views}",
        "i" = "Please choose from {.or {.strong {valid_views}}}"))
  }

  if (!is.null(roi_outline)){
    valid_rois <- locate_roi()$roi_label  
    bad_rois <- setdiff(roi_outline, valid_rois)
    if (length(bad_rois))
        vw_error(c("Invalid ROI{?s}: {bad_rois}",
        "i" = "Please choose from {.or {.strong {valid_rois[!is.na(valid_rois)]}}}"))
    if (is.null(fs_home)) {
      vw_error("We currently only support ROI outlines from FreeSurfer aparc, please provide `fs_home`")
    }
    roi_outline <- as.list(roi_outline)
  }

  lh <- check_hemi(lh, fs_template)
  rh <- check_hemi(rh, fs_template)

  lh_mask <- check_mask(lh_mask, fs_template)
  rh_mask <- check_mask(rh_mask, fs_template)

  where_is_my_mesh <- resolve_mesh(fs_template, fs_home)

  args <- list(lh = lh,
               rh = rh,
          lh_mask = lh_mask,
          rh_mask = rh_mask,
             vmin = vmin,
             vmax = vmax,
            views = as.list(views),
          surface = surface,
      bg_map_type = bg_map,
      fs_template = fs_template,
          fs_home = where_is_my_mesh,  # NULL: Python None (nilearn download)
         colorbar = colorbar,
   colorbar_label = colorbar_label,
   colorbar_width = colorbar_width,
             cmap = cmap,
        roi_names = roi_outline,
            title = title,
          cell_px = cell_px
  )

  if (!is.null(to_file)) {

    out_dir <- dirname(normalizePath(to_file, mustWork = FALSE))

    if (!dir.exists(out_dir))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    do.call(reticulate::py$vw_surf_static_plotly,
            c(args, list(output_file = to_file, dpi = as.integer(dpi))))
    
    vw_message("\u2714 Brain map saved to: {.file {to_file}}")

    return(invisible(to_file))

  } else {

    tmp_html <- tempfile(fileext = ".html")

    do.call(reticulate::py$vw_surf_interactive,
            c(args, list(output_html = tmp_html)))

    viewer <- getOption("viewer", utils::browseURL)
    viewer(tmp_html)
    vw_message("\u2714 Interactive brain map opened")

    return(invisible(tmp_html))
  }
}

# Session-level init flag
# Avoids querying reticulate::py (which is NULL before first Python call).
.vw_surf_env <- new.env(parent = emptyenv())
.vw_surf_env$ready <- FALSE

.vw_surf_init_py <- function() {
  if (isTRUE(.vw_surf_env$ready)) return(invisible(NULL))

  # Ensure reticulate is available
  require_packages('reticulate', call_fn = 'plot_vw_surf')

  req_pkgs <- c("nilearn", "numpy", "matplotlib", "plotly", "kaleido", 
                "choreographer", "logistro")
  
  # Declare requirements
  reticulate::py_require(req_pkgs)

  missing_pkgs <- Filter(
    function(x) !reticulate::py_module_available(x), req_pkgs)
  
  if (length(missing_pkgs) > 0) {
    vw_error(c(
      "Surface plotting requires Python packages that are not available:",
      "i" = "Missing: {.pkg {missing_pkgs}}",
      "i" = "{.pkg verywise} tries to auto-configure these, but this fails if you are using a custom Python setup (e.g., on an HPC cluster) or have {.envvar RETICULATE_PYTHON} set.",
      " " = "To fix this, choose one of the following options:",
      "*" = "Unset {.envvar RETICULATE_PYTHON} in your {.file .Renviron} and restart R to let {.pkg verywise} manage the environment.",
      "*" = "Install the missing packages in your active Python environment via the terminal:",
      " " = " {.code pip install {paste(missing_pkgs, collapse = ' ')}}",
      " " = " Then ensure your {.file .Renviron} points directly to that environment's Python:",
      " " = " {.code RETICULATE_PYTHON=/path/to/your/venv/bin/python}"
    ))
  }

  patch_file <- system.file("python", "patch_kaleido.py", package = "verywise")
  main_file <- system.file("python", "plot_vw_surf.py", package = "verywise")

  if (!nzchar(main_file))
    stop("verywise: could not find inst/python/plot_vw_surf.py - ",
         "try reinstalling the package.", call. = FALSE)
  
  # Apply the kaleido/orjson patch before loading the renderer
  reticulate::py_run_string(
    paste(readLines(patch_file), collapse = "\n")
  )
  
  py_code <- paste(readLines(main_file), collapse = "\n")

  tryCatch(
    reticulate::py_run_string(py_code),
    # reticulate::source_python(py_file),
    error = function(e) stop(
      "verywise: failed to load Python brain-plot renderer.\n",
      "Ensure nilearn, matplotlib, and numpy are installed.\n",
      "Original error: ", conditionMessage(e),
      call. = FALSE
    )
  )
  
  .vw_surf_env$ready <- TRUE
  invisible(NULL)
}