#' @title
#' Plot vertex-wise coefficient maps on a 3D cortical surface
#'
#' @description
#' Visualize vertex-wise statistical coefficient maps (e.g., LMM fixed effects)
#' on the FreeSurfer fsaverage cortical surface using Python's
#' **`nilearn`** and **`plotly`** backends via `{reticulate}`.
#'
#' Optionally, outlines specific cortical regions of interest (ROIs)
#' using FreeSurfer annotation information (i.e. 36 regions Desikan-Killiany atlas).
#'
#' @param res_dir Character. Path to the directory containing FreeSurfer-style
#'   vertex-wise result files (`*.mgh`).
#' @param term Character. Name of the model term to visualize (matches entries
#'   in `stack_names.txt`).
#' @param hemi Character. Hemisphere to plot (`"lh"` or `"rh"`).
#' @param measure Character. Surface measure, e.g. `'area'`, `'thickness'`,
#'   `'volume'`. Defaults to `'area'`.
#' @param surface Character. Surface mesh to plot on (`"inflated"` or `"pial"`). Defaults
#'   to `"inflated"`.
#' @param outline_rois Optional character vector. ROI names to outline on the
#'   surface plot. You can use call `verywise::locate_roi()` to see the names of all
#'   available ROIs.
#' @param fs_home Character. Path to the FreeSurfer installation (`$FREESURFER_HOME`).
#' @param show_in_browser Logical (default = TRUE). Open the Figure in browser.
#'   This can be invoked manually as well using `$show()`.
#'
#' @details
#' The function:
#' 1. Locates the appropriate MGH files for a given model term and hemisphere.
#' 2. Loads them with **nibabel** via `{reticulate}`.
#' 3. Displays the coefficient map interactively on the fsaverage inflated surface
#'    using **`nilearn`**'s `plot_surf_stat_map()` (Plotly engine).
#' 4. Optionally outlines regions from the FreeSurfer annotation (via `verywise`).
#'
#' Python dependencies (`nibabel`, `nilearn`, `matplotlib`, `plotly`, `numpy`)
#' must be installed in the active `{reticulate}` environment.
#'
#' @return A `plotly.graph_objects.Figure` (Python object) that can be displayed
#'   or further modified via `reticulate`.
#'
#' @examplesIf rlang::is_installed("reticulate") && dir.exists("~/results/fs_results")
#' # Example usage (requires FreeSurfer fsaverage and precomputed MGH maps)
#' fs_home <- Sys.getenv("FREESURFER_HOME")
#' surf_fig <- plot_vw_map(
#'   res_dir = "~/results/fs_results",
#'   term = "age",
#'   hemi = "lh",
#'   measure = "area",
#'   outline_rois = c("entorhinal", "precuneus"),
#'   fs_home = fs_home)
#'
#' @export
#' 
plot_vw_map <- function(res_dir, term, 
                        hemi = c("lh", "rh"), measure = 'area', 
                        surface = c("inflated", "pial"),
                        outline_rois = NULL, fs_home = Sys.getenv("FREESURFER_HOME"), 
                        show_in_browser = TRUE
                      ) {

  # Match mesh arguments
  hemi <- match.arg(hemi)
  hemi_name <- if (hemi == "lh") "left" else "right"
  surface <- match.arg(surface)

  # Validate paths
  if (!dir.exists(res_dir))
    stop("Results directory does not exist: ", res_dir)
  if (!dir.exists(fs_home))
    stop("Freesurfer home not found: ", fs_home)

  # Retrieve fsaverage surface meshes
  surf_mesh <- file.path(fs_home, "subjects", "fsaverage", "surf", 
                         paste(hemi, surface, sep="."))

  # Find term stack index
  stack_file <- file.path(res_dir, "stack_names.txt")
  if (!file.exists(stack_file)) {
    stop(sprintf("Cannot find the `stack_names.txt` in `%s`.", res_dir),
           " Results folder is incorrect or corrupted.")
  }

  # Read existing file
  stack_ids <- utils::read.table(stack_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (!term %in% stack_ids$stack_name) {
    stop("Term '", term, "' not found in stack_names.txt")
  }
  # Extract stack
  stack <- paste0('stack', stack_ids[ stack_ids$stack_name == term, 'stack_number'])

  # Ensure Python deps available
  reticulate::py_require(c("nibabel", "nilearn", "matplotlib", "plotly"))

  # Import Python modules
  nb <- reticulate::import("nibabel")
  nl <- reticulate::import("nilearn.plotting")
  np <- reticulate::import("numpy")

  # Load coefficients .MGH file
  coef_map_file <- file.path(res_dir, paste(hemi, measure, stack, "coef.mgh", sep="."))
  coef_map <- nb$load(coef_map_file)$get_fdata()

  # Load cluster info
  ocn_map_file <- file.path(res_dir, paste(hemi, measure, stack, "cache.th30.abs.sig.ocn.mgh", sep = "."))
  ocn_map <- nb$load(ocn_map_file)$get_fdata()

  clusters <- table(ocn_map)
  print(clusters)
  
  surf_fig <- nl$plot_surf_stat_map(
    surf_mesh = surf_mesh,
    stat_map = coef_map,
    cmap = "coolwarm",
    colorbar = TRUE,
    hemi = hemi_name,
    engine = "plotly",
    vmin = min(coef_map), 
    vmax = max(coef_map),
    threshold = 0.0001,
    cbar_tick_format = ".3f",
    view = 'medial',
    title = paste('Effect of', term, 'on', measure, '-', hemi_name, 'hemisphere')
  )

  if (!is.null(outline_rois)) {
    roi_info <- locate_roi(hemi = hemi)
    roi_colors <- roi_info$roi_color[match(outline_rois, roi_info$roi_label)]

    roi_mask <- Reduce(`+`, lapply(seq_along(outline_rois), function(i) {
      mask_n <- locate_roi(outline_rois[i], hemi = hemi, verbose = FALSE)
      as.numeric(mask_n) * i
    }))

    surf_fig$add_contours(
      roi_map = np$array(roi_mask),
      levels = as.list(seq_along(outline_rois)),
      labels = as.list(outline_rois),
      lines = lapply(roi_colors, function (c) list(width = 6, color = c))
    )
  }

  # Launch browser 
  if (show_in_browser) surf_fig$show() 

  invisible(surf_fig)
}


# plot_vw_vec <- function(vec, vec_name = NULL, 
#                         hemi = c("lh", "rh"), measure = 'area', 
#                         surface = c("inflated", "pial"),
#                         outline_rois = NULL, fs_home = Sys.getenv("FREESURFER_HOME"), 
#                         show_in_browser = TRUE
#                       ) {

#   # Match mesh arguments
#   hemi <- match.arg(hemi)
#   hemi_name <- if (hemi == "lh") "left" else "right"
#   surface <- match.arg(surface)

#   if (!dir.exists(fs_home))
#     stop("Freesurfer home not found: ", fs_home)

#   # Retrieve fsaverage surface meshes
#   surf_mesh <- file.path(fs_home, "subjects", "fsaverage", "surf", 
#                          paste(hemi, surface, sep="."))

#   # Ensure Python deps available
#   reticulate::py_require(c("nibabel", "nilearn", "matplotlib", "plotly"))

#   # Import Python modules
#   nb <- reticulate::import("nibabel")
#   nl <- reticulate::import("nilearn.plotting")
#   np <- reticulate::import("numpy")
  
#   surf_fig <- nl$plot_surf_stat_map(
#     surf_mesh = surf_mesh,
#     stat_map = vec,
#     cmap = "coolwarm",
#     colorbar = TRUE,
#     hemi = hemi_name,
#     engine = "plotly",
#     vmin = min(vec), 
#     vmax = max(vec),
#     # threshold = 0.0001,
#     cbar_tick_format = ".3f",
#     view = 'medial',
#     title = paste('Value of', vec_name, 'on', measure, '-', hemi_name, 'hemisphere')
#   )

#   if (!is.null(outline_rois)) {
#     roi_info <- locate_roi(hemi = hemi)
#     roi_colors <- roi_info$roi_color[match(outline_rois, roi_info$roi_label)]

#     roi_mask <- Reduce(`+`, lapply(seq_along(outline_rois), function(i) {
#       mask_n <- locate_roi(outline_rois[i], hemi = hemi, verbose = FALSE)
#       as.numeric(mask_n) * i
#     }))

#     surf_fig$add_contours(
#       roi_map = np$array(roi_mask),
#       levels = as.list(seq_along(outline_rois)),
#       labels = as.list(outline_rois),
#       lines = lapply(roi_colors, function (c) list(width = 6, color = c))
#     )
#   }

#   # Launch browser 
#   if (show_in_browser) surf_fig$show() 

#   invisible(surf_fig)
# }

# #' @title Plot a vertex-wise map on brain surface
# #' 
# #' @description Projects a vector of vertex-wise values onto a FreeSurfer
# #'   template surface and saves a multi-view PNG. Uses \pkg{fsbrain} and
# #'   \pkg{rgl} (instead of calling Python via reticulate).
# #'
# #' @param surf_data Numeric vector of vertex-wise values (length = 163842 for fsaverage).
# #' @param filename Character. Output PNG path. Must end in \code{".png"}.
# #'   Default: \code{file.path(tempdir(), "vw_map.png")}.
# #' @param surface Character. Surface type: \code{"inflated"} (default),
# #'   \code{"pial"}, \code{"white"}, or \code{"smoothwm"}.
# #' @param cmap Character or character vector. A named palette from
# #'   \pkg{RColorBrewer} (e.g. \code{"RdBu"}) / \pkg{grDevices}
# #'   (e.g. \code{"viridis"}), or a vector of hex color codes for a
# #'   custom colormap. Default: \code{"RdBu"}.
# #' @param title Character or \code{NULL}. Optional plot title.
# #' @param subjects_dir Character. FreeSurfer subjects directory.
# #'   Defaults to the \code{SUBJECTS_DIR} environment variable.
# #' @param subject Character. Template subject. Default: \code{"fsaverage5"}.
# #' @param views Character vector. Camera views to tile in the output image.
# #'   Default: \code{c("lateral", "medial")}.
# #' @param clim Numeric vector of length 2 \code{c(min, max)} for the color
# #'   scale. If \code{NULL} (default) a symmetric range around 0 is used.
# #' @param thresh Numeric. Vertices with \code{abs(value) < thresh} are shown
# #'   in the background color (i.e. masked). Default: \code{0} (no masking).
# #'
# #' @return Invisibly returns \code{filename}.
# #' @seealso \code{\link{plot_vw_map}} for the Python-backed version.
# #' @export
# plot_vw_map_r <- function(surf_data,
#                            filename     = file.path(tempdir(), "vw_map.png"),
#                            surface      = "inflated",
#                            cmap         = "RdBu",
#                            title        = NULL,
#                            subjects_dir = Sys.getenv("SUBJECTS_DIR"),
#                            fs_template = "fsaverage5",
#                            views        = c("lateral", "medial"),
#                            clim         = NULL,
#                            thresh       = 0) {

#   # ---- dependency checks ----
#   if (!requireNamespace("fsbrain", quietly = TRUE))
#     stop("Install 'fsbrain': install.packages('fsbrain')", call. = FALSE)
#   if (!requireNamespace("rgl", quietly = TRUE))
#     stop("Install 'rgl': install.packages('rgl')", call. = FALSE)

#   # ---- input validation ----
#   if (!is.numeric(surf_data) || is.matrix(surf_data))
#     stop("`surf_data` must be a numeric vector.", call. = FALSE)
#   # if (!grepl("\\.png$", filename, ignore.case = TRUE))
#   #   stop("`filename` must have a .png extension.", call. = FALSE)
#   # if (subjects_dir == "")
#     # stop("Set `subjects_dir` or the SUBJECTS_DIR environment variable.",
#     #      call. = FALSE)

#   n_verts <- length(surf_data)

#   # ---- hemisphere split ----
#   # lh_data <- surf_data[seq_len(half)]
#   # rh_data <- surf_data[seq(half + 1L, n_verts)]

#   # ---- optional threshold masking ----
#   # if (thresh > 0) {
#   #   lh_data[abs(lh_data) < thresh] <- NA
#   #   rh_data[abs(rh_data) < thresh] <- NA
#   # }

#   # ---- color limits ----
#   if (is.null(clim)) {
#     abs_max <- max(abs(surf_data), na.rm = TRUE)
#     clim <- c(-abs_max, abs_max)
#   }

#   # ---- build palette function ----
#   diverging_rb <- c("RdBu", "RdYlBu", "Spectral", "PiYG",
#                     "PRGn", "BrBG", "RdGy", "PuOr")
#   if (length(cmap) == 1L && is.character(cmap)) {
#     if (cmap %in% diverging_rb) {
#       pal_fn <- function(n)
#         rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, cmap))(n))
#     } else {
#       pal_fn <- function(n) grDevices::hcl.colors(n, palette = cmap)
#     }
#   } else {
#     # custom hex vector
#     pal_fn <- grDevices::colorRampPalette(cmap)
#   }

#   # ---- render with fsbrain ----
#   cm <- fsbrain::vis.data.on.fsaverage(
#     vis_subject_id   = fs_template,
#     morph_data_lh    = surf_data,
#     morph_data_rh    = surf_data,
#     subjects_dir     = NULL,
#     surface          = "pial",
#     views            = 'si',
#     # makecmap_options = list(
#     #   colFn  = pal_fn,
#     #   n      = 256L,
#     #   range  = clim,
#     #   symm   = TRUE
#     # ),
#     draw_colorbar    = TRUE #
#     # title            = if (!is.null(title)) title else ""
#   )
  
#   rgl::rglwidget()
#   # ---- export PNG ----
#   dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
#   fsbrain::export(cm,
#                   img_only   = TRUE,
#                   output_img = filename,
#                   views      = views)

#   message("Saved: ", filename)
#   invisible(filename)
# }
