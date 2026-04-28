# # ── helpers ────────────────────────────────────────────────────────────────────

# .validate_morph <- function(hemi, fs_template) {
#   if (is.null(hemi)) return(NULL)
#   hemi <- as.numeric(hemi)
#   n    <- count_vertices(fs_template)
#   if (length(hemi) != n)
#     stop("Morphometry vector length (", length(hemi), ") does not match ",
#          fs_template, " vertex count (", n, ").")
#   hemi
# }

# .resolve_views <- function(views) {
#   presets <- c("t4", "t9", "si")
#   if (is.character(views) && length(views) == 1 && views %in% presets)
#     fsbrain::get.view.angle.names(angle_set = views)
#   else
#     views
# }

# # Monkey-patch rgl::rgl.snapshot → rgl::snapshot3d(webshot=TRUE) for the
# # duration of `expr`. Needed because fsbrain calls rgl.snapshot() internally.
# .with_webshot_snapshot <- function(expr) {
#   ns  <- asNamespace("rgl")
#   old <- get("rgl.snapshot", envir = ns)
#   on.exit(utils::assignInNamespace("rgl.snapshot", old, ns = "rgl"), add = TRUE)

#   utils::assignInNamespace("rgl.snapshot", ns = "rgl",
#                            function(filename, fmt = "png", top = TRUE) {
#                              rgl::snapshot3d(filename = filename, fmt = fmt, webshot = TRUE)
#                            }
#   )
#   force(expr)
# }

# # Maual static tiling??
# .plot_static_fsbrain <- function(cm, view_angles, colorbar_legend, output_file) {
#   tmp_dir <- tempfile()
#   dir.create(tmp_dir)

#   # Render each view angle as a separate PNG using webshot2 directly
#   imgs <- lapply(seq_along(view_angles), function(i) {
#     angle <- view_angles[[i]]
#     f     <- file.path(tmp_dir, sprintf("view_%02d.png", i))

#     fsbrain::vis.coloredmeshes.rotating(cm, angle = angle, output_img = f)
#     # ↑ fsbrain sets the viewpoint; we snapshot ourselves:
#     rgl::view3d(usermatrix = fsbrain::get.view.angle.matrix(angle))
#     rgl::snapshot3d(f, webshot = TRUE)
#     rgl::close3d()
#     magick::image_read(f)
#   })

#   # Stitch horizontally, add colourbar, save
#   combined <- magick::image_append(do.call(c, imgs))
#   magick::image_write(combined, output_file)
# }


# # ── main function ──────────────────────────────────────────────────────────────

# #' @title Plot vertex-wise maps on brain surfaces
# #'
# #' @description
# #' Render a vertex-wise scalar map on an fsaverage template surface using
# #' \pkg{fsbrain}. Works cross-platform without XQuartz:
# #' \itemize{
# #'   \item **Interactive** (no \code{output_file}): returns an
# #'         \code{\link[rgl]{rglwidget}} that opens in the RStudio Viewer,
# #'         Positron, or any browser. No display server required.
# #'   \item **Static PNG** (\code{output_file} supplied): renders via
# #'         \code{webshot2}/headless-Chrome on macOS and other platforms
# #'         without native OpenGL snapshot support; falls back to direct
# #'         OpenGL snapshot on Linux with a display.
# #' }
# #'
# #' @param lh,rh Numeric vector \strong{or} path to a \code{.mgh}/\code{.gii}
# #'   file, or \code{NULL}. Vertex-wise values for the left / right hemisphere.
# #'   At least one must be non-\code{NULL}.
# #' @param surface \code{"inflated"} (default), \code{"pial"}, or
# #'   \code{"white"}.
# #' @param cmap A \code{colorRampPalette}-style function.
# #'   Default: blue–white–red diverging palette.
# #' @param bg_map Background shading: \code{"sulc"} (default),
# #'   \code{"curv"}, or \code{"none"}.
# #' @param symm Logical. Use a symmetric (zero-centred) colour scale?
# #'   Default \code{TRUE}.
# #' @param views Camera angles. \code{"t4"} (4 standard views, default),
# #'   \code{"t9"} (9 views), \code{"si"} (superior/inferior), or a character
# #'   vector of individual angle names.
# #' @param colorbar_legend Character label for the colour bar, or \code{NULL}.
# #' @param output_file Path ending in \code{.png} for static export, or
# #'   \code{NULL} for interactive mode.
# #' @param fs_home FreeSurfer installation root. Defaults to
# #'   \env{FREESURFER_HOME}.
# #' @param fs_template fsaverage template name. Default \code{"fsaverage"}.
# #'
# #' @return In interactive mode, an \code{rglwidget} object (auto-prints as a
# #'   WebGL scene). In static mode, the \code{output_file} path invisibly.
# #'
# #' @export
# plot_vw_surf_r <- function(
#     lh              = NULL,
#     rh              = NULL,
#     surface         = c("inflated", "pial", "white"),
#     cmap            = grDevices::colorRampPalette(c("blue", "white", "red")),
#     bg_map          = c("sulc", "curv", "none"),
#     symm            = TRUE,
#     views           = "t4",
#     colorbar_legend = NULL,
#     output_file     = NULL,
#     fs_home         = Sys.getenv("FREESURFER_HOME"),
#     fs_template     = "fsaverage"
# ) {

#   ## ── validation ─────────────────────────────────────────────────────────────
#   if (is.null(lh) && is.null(rh))
#     stop("Supply at least one of `lh` or `rh`.")
#   if (!is.null(output_file) && !grepl("\\.png$", output_file, ignore.case = TRUE))
#     stop("`output_file` must end in '.png'.")
#   if (!dir.exists(file.path(fs_home, "subjects", fs_template)))
#     stop("Template '", fs_template, "' not found under ", fs_home, "/subjects/")

#   surface <- match.arg(surface)
#   bg_map  <- match.arg(bg_map)

#   lh <- .validate_morph(lh, fs_template)
#   rh <- .validate_morph(rh, fs_template)

#   subjects_dir  <- file.path(fs_home, "subjects")
#   cmap_opts     <- list(colFn = cmap, symm = symm)
#   view_angles   <- .resolve_views(views)
#   bg_arg        <- if (bg_map == "none") NULL else bg_map

#   ## ── always use null device: no XQuartz, no display server needed ──────────
#   withr::with_options(list(rgl.useNULL = TRUE), {

#     if (is.null(output_file)) {
#       ## ── interactive: WebGL widget → Viewer / browser ──────────────────────
#       cm <- fsbrain::vis.symmetric.data.on.subject(
#         subjects_dir,
#         fs_template,
#         morph_data_lh    = lh,
#         morph_data_rh    = rh,
#         surface          = surface,
#         bg               = bg_arg,
#         makecmap_options = cmap_opts
#         # no 'no_vis' rglaction: scene is rendered to the null device so
#         # rglwidget() can serialise it
#       )

#       widget <- rgl::rglwidget()
#       rgl::close3d()
#       return(widget)          # auto-prints as WebGL in any IDE or browser

#     } else {
#       ## ── static PNG: headless render + webshot2 snapshot ───────────────────
#       out_dir <- dirname(output_file)
#       if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

#       cm <- fsbrain::vis.symmetric.data.on.subject(
#         subjects_dir,
#         fs_template,
#         morph_data_lh    = lh,
#         morph_data_rh    = rh,
#         surface          = surface,
#         bg               = bg_arg,
#         rglactions       = list(no_vis = TRUE), # suppress window; keep scene
#         makecmap_options = cmap_opts
#       )

#       # webshot2 available → use it (macOS, HPC without Xvfb, etc.)
#       # absent  → fall back to native rgl.snapshot (Linux + Xvfb)
#       if (rlang::is_installed("webshot2")) {
#         .with_webshot_snapshot(
#           fsbrain::vis.export.from.coloredmeshes(
#             cm,
#             colorbar_legend  = colorbar_legend,
#             output_img       = output_file,
#             view_angles      = view_angles,
#             background_color = "white"
#           )
#         )
#       } else {
#         .plot_static_fsbrain(
#           cm,
#           colorbar_legend  = colorbar_legend,
#           output_img       = output_file,
#           view_angles      = view_angles,
#           background_color = "white"
#         )
#       }

#       message("\u2713  Saved to: ", output_file)
#       return(invisible(output_file))
#     }
#   })
# }