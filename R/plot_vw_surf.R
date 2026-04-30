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
#' @param lh Numeric vector, \code{.mgh}/\code{.gii} file path, or
#'   \code{NULL}.
#' @param rh Numeric vector, \code{.mgh}/\code{.gii} file path, or
#'   \code{NULL}. At least one of \code{lh}/\code{rh} must be supplied.
#' @param fs_template fsaverage template: \code{"fsaverage3"} through
#'   \code{"fsaverage"}. Must match the length of lh/rh. Default \code{"fsaverage"}.
#' @param fs_home (optional) location of FreeSurfer home for templates.
#' @param surface Surface mesh: \code{"inflated"} (default), \code{"pial"},
#'   or \code{"white"}.
#' @param cmap Matplotlib colormap name.
#' @param bg_map Background shading: \code{"sulc"} (default),
#'   \code{"curv"}, or \code{"none"}.
#' @param vmin,vmax Numeric colour limits, or \code{NULL} for automatic
#'   symmetric scaling.
#' @param threshold Absolute-value masking threshold, or \code{NULL}
#'   (no masking).
#' @param views Character vector of camera angles - any subset of
#'   \code{"lateral"}, \code{"medial"}, \code{"dorsal"}, \code{"ventral"},
#'   \code{"anterior"}, \code{"posterior"}. Default: \code{"all"}.
#' @param colorbar Logical. Draw a shared colour bar? Default \code{TRUE}.
#' @param colorbar_label Character label for the colour bar axis, or
#'   \code{NULL}.
#' @param title Character figure title, or \code{NULL}.
#' @param to_file Path ending in \code{.png} for static export, or
#'   \code{NULL} for interactive HTML mode.
#' @param dpi Integer output resolution for static PNG. Default \code{150L}.
#'
#' @return Invisibly: \code{to_file} path (static) or the temp HTML
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
#'   cmap = "RdBu_r",
#'   threshold = 0.05,
#'   views = c("lateral", "medial", "dorsal", "ventral"),
#'   title = "Effect of age on cortical thickness",
#'   to_file = "figures/age_thickness.png",
#'   dpi = 300L
#' )
#' }
#'
#' @export
plot_vw_surf <- function(
    lh = NULL,
    rh = NULL,
    fs_template = "fsaverage",
    fs_home = NULL,
    surface = c("inflated", "pial", "white"),
    cmap = NULL,
    bg_map = c("sulc", "curv", "none"),
    vmin = NULL,
    vmax = NULL,
    threshold = NULL,
    views = 'all',
    colorbar = TRUE,
    colorbar_label = NULL,
    title = NULL,
    to_file = NULL,
    dpi = 150L) {

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
  if (views == 'all') {
    views <- if (!is.null(to_file)) valid_views else 'lateral'
  } else {
    bad_views <- setdiff(views, valid_views)
    if (length(bad_views))
        vw_error(c("Invalid view{?s}: {bad_views}",
        "i" = "Please choose from {.or {.strong {valid_views}}}"))
  }

  lh <- check_hemi(lh, fs_template)
  rh <- check_hemi(rh, fs_template)

  where_is_my_mesh <- resolve_mesh(fs_template, fs_home)

  # --- initialise Python renderer (once per session) -----------------------
  reticulate::py_require(c("nilearn", "matplotlib", "numpy", "plotly", "kaleido"))
  .vw_surf_init_py()

  common <- list(
              lh = lh,
              rh = rh,
         surface = surface,
     bg_map_type = bg_map,
            cmap = cmap,
            vmin = vmin,
            vmax = vmax,
       threshold = threshold,
          views = as.list(views),
        colorbar = colorbar,
  colorbar_label = colorbar_label,
           title = title,
     fs_template = fs_template,
         fs_home = where_is_my_mesh  # NULL: Python None (nilearn download)
  )

  if (!is.null(to_file)) {

    out_dir <- dirname(normalizePath(to_file, mustWork = FALSE))

    if (!dir.exists(out_dir))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    do.call(reticulate::py$vw_surf_static_plotly,
            c(common, list(output_file = to_file, dpi = as.integer(dpi))))
    
    vw_message("\u2714 Brain map saved to: {.file {to_file}}")

    return(invisible(to_file))

  } else {

    reticulate::py_require("plotly")
    tmp_html <- tempfile(fileext = ".html")

    do.call(reticulate::py$vw_surf_interactive,
            c(common, list(output_html = tmp_html)))

    viewer <- getOption("viewer", utils::browseURL)
    viewer(tmp_html)
    vw_message("\u2714 Interactive brain map opened")

    return(invisible(tmp_html))
  }
}

check_hemi <- function(hemi, fs_template) {
  
  if (is.null(hemi)) return(hemi)  # empty or file path: skip check
  
  hemi_name <- deparse(substitute(hemi))
  
  if (is.character(hemi)) {
    if (!file.exists(hemi)) vw_error("{hemi_name} file not found: {.file {hemi}}")
  
    # nilearn.surface.load_surf_data()
    nilearn_surf_ext <- c("\\.mgh$", "\\.mgz$", "\\.gii$", "\\.nii$",
                          "\\.nii\\.gz$", "\\.npy$", "\\.txt$", "\\.csv$")
    
    if (!any(grepl(paste(nilearn_surf_ext, collapse = "|"), hemi, ignore.case = TRUE)))
        vw_message(c("!" = "{hemi_name}: unrecognised file extension in {.file {basename(hemi)}}.",
                     " " = "nilearn will attempt to load it anyway but things may get weird.",
                     ">" = "try reading it in yourself and providing a vector instead, or using 
                     one of the supported extensions (e.g. .mgh, .mgz, .csv... see `nilearn.surface.load_surf_data()`)"
    ))

    return(hemi)
  }

  hemi <- as.numeric(hemi)

  n_vert <- count_vertices(fs_template)

  if (length(hemi) != n_vert) {
      vw_message("!" = "{hemi_name} vector length ({.warn {length(hemi)}}) does not match 
      {fs_template} template ({n_vert}), I will try to subset it.")
  }

  hemi
}

# Session-level init flag
# Avoids querying reticulate::py (which is NULL before first Python call).
.vw_surf_env <- new.env(parent = emptyenv())
.vw_surf_env$ready <- FALSE

.vw_surf_init_py <- function() {
  if (isTRUE(.vw_surf_env$ready)) return(invisible(NULL))

  py_file <- system.file("python", "plot_vw_surf.py", package = "verywise")
  if (!nzchar(py_file))
    stop("verywise: could not find inst/python/plot_vw_surf.py - ",
         "try reinstalling the package.", call. = FALSE)
  
  py_code <- paste(readLines(py_file), collapse = "\n")

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    vw_error(c(
      "Surface plotting requires {.pkg reticulate}.",
      "i" = "Install it with {.code install.packages('reticulate')}."
    ))
  }

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

#   tryCatch(
#     reticulate::py_run_string(.VW_SURF_PY),
#     error = function(e) stop(
#       "verywise: failed to initialise the Python brain-plot renderer.\n",
#       "Ensure nilearn, matplotlib, and numpy are installed in the active ",
#       "reticulate Python environment.\n",
#       "Original error: ", conditionMessage(e),
#       call. = FALSE
#     )
#   )
  .vw_surf_env$ready <- TRUE
  invisible(NULL)
}


# Mesh-resolution helper (R side)
# Returns the FreeSurfer home path when the template is found locally,
# or NULL to signal Python to use the nilearn download+cache path (via fetch_surf_fsaverage)
resolve_mesh <- function(fs_template, fs_home, verbose = TRUE) {
  
  if (is.null(fs_home)) {
    fs_home <- Sys.getenv("FREESURFER_HOME")
  }
  
  if (nzchar(fs_home)) {
    surf_dir <- file.path(fs_home, "subjects", fs_template, "surf")
    if (dir.exists(surf_dir)) {
      vw_message("Using local FreeSurfer mesh from {fs_home}", verbose = verbose, type = 'note')
      return(fs_home)
    }
    vw_message("{fs_template} surface mesh not found in {.path $FREESURFER_HOME/subjects/}
         {cli::symbol$arrow_right} falling back to nilearn.")
  }

  # Warn only when the template is not yet cached
  cache <- file.path(path.expand("~"), "nilearn_data", fs_template)
  if (!dir.exists(cache) || length(list.files(cache, recursive = TRUE)) == 0L) {
     vw_message(c("i" = "Will download {fs_template} surface mesh. This may take up to a minute on the first run.",
                  " " = "Meshes will be then cached in {.path ~/nilearn_data/{fs_template}}",
                  " " = "Alternatively, provide a path to FreeSurfer (via `fs_home` argument or environment variables)"))
  }

  NULL 
}