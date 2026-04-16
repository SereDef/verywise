# R/plot_vw_surf.R
# ──────────────────────────────────────────────────────────────────────────────
# Brain surface plotting via Python/nilearn.
# Fully headless — no XQuartz, rgl, or display server required.
# ──────────────────────────────────────────────────────────────────────────────

.VW_SURF_PY <- '
import io, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import matplotlib.image as mpimg
from nilearn import datasets as _nids
from nilearn import plotting as _nlplot
from nilearn.surface import load_surf_data as _lsd


# ── helpers ───────────────────────────────────────────────────────────────────

def _to_array(x):
    if x is None:
        return None
    if isinstance(x, str):
        return np.ravel(_lsd(x)).astype(np.float32)
    return np.asarray(x, dtype=np.float32).ravel()


def _auto_limits(arrays):
    vals = np.concatenate(
        [a[np.isfinite(a)] for a in arrays if a is not None]
    )
    mx = float(np.max(np.abs(vals))) if vals.size else 1.0
    return -mx, mx


def _get_mesh_paths(fs_template, surface, bg_map_type, fs_home):
    if fs_home:
        surf_dir = os.path.join(fs_home, "subjects", fs_template, "surf")
        if os.path.isdir(surf_dir):
            mlh = os.path.join(surf_dir, f"lh.{surface}")
            mrh = os.path.join(surf_dir, f"rh.{surface}")
            blh = os.path.join(surf_dir, f"lh.{bg_map_type}") if bg_map_type != "none" else None
            brh = os.path.join(surf_dir, f"rh.{bg_map_type}") if bg_map_type != "none" else None
            missing = [p for p in filter(None, [mlh, mrh, blh, brh])
                       if not os.path.exists(p)]
            if missing:
                raise FileNotFoundError(
                    "Expected FreeSurfer surface file(s) not found: " +
                    ", ".join(missing)
                )
            return mlh, mrh, blh, brh
    # Fall back to nilearn download + permanent cache
    fsavg    = _nids.fetch_surf_fsaverage(mesh=fs_template)
    surf_key = "infl" if surface == "inflated" else surface
    mlh = getattr(fsavg, f"{surf_key}_left")
    mrh = getattr(fsavg, f"{surf_key}_right")
    blh = getattr(fsavg, f"{bg_map_type}_left")  if bg_map_type != "none" else None
    brh = getattr(fsavg, f"{bg_map_type}_right") if bg_map_type != "none" else None
    return mlh, mrh, blh, brh


def _fig_to_arr(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
    buf.seek(0)
    arr = mpimg.imread(buf)
    buf.close()
    plt.close(fig)
    return arr


def _unwrap_plotly(fp):
    # nilearn returns PlotlySurfaceFigure; unwrap to raw plotly Figure
    return fp.figure if hasattr(fp, "figure") else fp


# ── static PNG renderer ───────────────────────────────────────────────────────

def vw_surf_static(lh, rh, surface, bg_map_type, cmap, vmin, vmax,
                   threshold, views, colorbar, colorbar_label,
                   title, output_file, dpi, fs_template, fs_home):

    lh, rh = _to_array(lh), _to_array(rh)
    if vmin is None or vmax is None:
        vmin, vmax = _auto_limits([lh, rh])

    mlh, mrh, blh, brh = _get_mesh_paths(fs_template, surface, bg_map_type, fs_home)
    hemis  = []
    if lh is not None: hemis.append(("left",  lh,  mlh, blh))
    if rh is not None: hemis.append(("right", rh,  mrh, brh))
    thresh = float(threshold) if threshold is not None else None

    panels = []
    for hname, data, mesh, bg in hemis:
        row = []
        for view in views:
            fp = _nlplot.plot_surf_stat_map(
                surf_mesh  = mesh,
                stat_map   = data,
                hemi       = hname,
                view       = view,
                cmap       = cmap,
                vmin       = vmin,
                vmax       = vmax,
                threshold  = thresh,
                bg_map     = bg,
                bg_on_data = bg is not None,
                colorbar   = False,
                engine     = "matplotlib",
            )
            row.append(_fig_to_arr(fp))
        panels.append(row)

    n_rows, n_cols = len(hemis), len(views)
    n_grid = n_cols + (1 if colorbar else 0)
    wrats  = [1.0] * n_cols + ([0.04] if colorbar else [])
    fig_w  = n_cols * 3.2 + (0.5 if colorbar else 0.0)
    fig_h  = n_rows * 2.5 + (0.6 if title else 0.15)

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white", dpi=dpi)
    gs  = gridspec.GridSpec(
        n_rows, n_grid,
        width_ratios = wrats,
        hspace = 0.02, wspace = 0.02,
        left   = 0.07,
        right  = 0.92 if colorbar else 0.99,
        top    = 0.92 if title else 0.98,
        bottom = 0.02,
    )
    for r, (hname, *_) in enumerate(hemis):
        for c, view in enumerate(views):
            ax = fig.add_subplot(gs[r, c])
            ax.imshow(panels[r][c])
            ax.axis("off")
            if c == 0:
                ax.set_ylabel(hname.capitalize(), fontsize=11,
                              fontweight="bold", rotation=90,
                              va="center", labelpad=6)
            if r == 0:
                ax.set_title(view.capitalize(), fontsize=10, pad=4)

    if colorbar:
        cax  = fig.add_subplot(gs[:, n_cols])
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        sm   = mcm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cb   = fig.colorbar(sm, cax=cax)
        if colorbar_label:
            cb.set_label(colorbar_label, fontsize=9)
        cb.ax.tick_params(labelsize=8)

    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", y=0.99)

    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    fig.savefig(output_file, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return output_file


# ── interactive HTML renderer ─────────────────────────────────────────────────

def vw_surf_interactive(lh, rh, surface, bg_map_type, cmap, vmin, vmax,
                        threshold, views, colorbar, colorbar_label,
                        title, output_html, fs_template, fs_home):
    # Lazy import — missing plotly never breaks the static path
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    lh, rh = _to_array(lh), _to_array(rh)
    if vmin is None or vmax is None:
        vmin, vmax = _auto_limits([lh, rh])

    mlh, mrh, blh, brh = _get_mesh_paths(fs_template, surface, bg_map_type, fs_home)
    hemis  = []
    if lh is not None: hemis.append(("left",  lh,  mlh, blh))
    if rh is not None: hemis.append(("right", rh,  mrh, brh))
    thresh = float(threshold) if threshold is not None else None

    n_rows, n_cols = len(hemis), len(views)
    specs  = [[{"type": "scene"}] * n_cols for _ in range(n_rows)]
    titles = [f"{h.capitalize()} - {v.capitalize()}"
              for h, *_ in hemis for v in views]

    fig_main = make_subplots(
        rows = n_rows, cols = n_cols,
        specs = specs,
        subplot_titles = titles,
        horizontal_spacing = 0.02,
        vertical_spacing   = 0.06,
    )

    first_bar = True
    for r, (hname, data, mesh, bg) in enumerate(hemis):
        for c, view in enumerate(views):
            fp  = _nlplot.plot_surf_stat_map(
                surf_mesh  = mesh,
                stat_map   = data,
                hemi       = hname,
                view       = view,
                cmap       = cmap,
                vmin       = vmin,
                vmax       = vmax,
                threshold  = thresh,
                bg_map     = bg,
                bg_on_data = bg is not None,
                colorbar   = False,
                engine     = "plotly",
            )
            raw = _unwrap_plotly(fp)
            # Plotly subplot scene keys: scene, scene2, scene3, ...
            idx = r * n_cols + c
            sk  = "scene" if idx == 0 else f"scene{idx + 1}"

            for trace in raw.data:
                if colorbar and first_bar:
                    trace.update(
                        showscale = True,
                        colorbar  = dict(
                            title      = dict(text=colorbar_label or "", side="right"),
                            thickness  = 15,
                            len        = 0.8,
                            x          = 1.02,
                            tickformat = ".2g",
                        ),
                    )
                    first_bar = False
                else:
                    trace.update(showscale=False)
                fig_main.add_trace(trace, row=r + 1, col=c + 1)

            if hasattr(raw.layout, "scene"):
                fig_main.update_layout(**{sk: raw.layout.scene.to_plotly_json()})

    fig_main.update_layout(
        title_text      = title or "",
        title_font_size = 14,
        paper_bgcolor   = "white",
        plot_bgcolor    = "white",
        height          = 380 * n_rows,
        width           = 400 * n_cols + (80 if colorbar else 0),
    )
    fig_main.write_html(output_html, include_plotlyjs="cdn",
                        full_html=True, auto_open=False)
    return output_html
'


# ── session-level init flag ────────────────────────────────────────────────────
# Avoids querying reticulate::py (which is NULL before first Python call).
.vw_surf_env       <- new.env(parent = emptyenv())
.vw_surf_env$ready <- FALSE

.vw_surf_init_py <- function() {
  if (isTRUE(.vw_surf_env$ready)) return(invisible(NULL))
  tryCatch(
    reticulate::py_run_string(.VW_SURF_PY),
    error = function(e) stop(
      "verywise: failed to initialise the Python brain-plot renderer.\n",
      "Ensure nilearn, matplotlib, and numpy are installed in the active ",
      "reticulate Python environment.\n",
      "Original error: ", conditionMessage(e),
      call. = FALSE
    )
  )
  .vw_surf_env$ready <- TRUE
  invisible(NULL)
}


# ── mesh-resolution helper (R side) ──────────────────────────────────────────
# Returns the FreeSurfer home path when the template is found locally,
# or NULL to signal Python to use the nilearn download+cache path.
.resolve_mesh <- function(fs_template) {
  fs_home <- Sys.getenv("FREESURFER_HOME")

  if (nzchar(fs_home)) {
    surf_dir <- file.path(fs_home, "subjects", fs_template, "surf")
    if (dir.exists(surf_dir)) {
      message("\u2714 Using local FreeSurfer mesh: ", surf_dir)
      return(fs_home)
    }
    message("\u2139  '", fs_template,
            "' not found under $FREESURFER_HOME/subjects/ — falling back to nilearn.")
  }

  # Warn only when the template is not yet cached
  cache <- file.path(path.expand("~"), "nilearn_data", fs_template)
  if (!dir.exists(cache) ||
      length(list.files(cache, recursive = TRUE)) == 0L) {
    message(
      "\u25ba Downloading '", fs_template, "' surface mesh (one-time setup).\n",
      "  Meshes will be cached in ~/nilearn_data/", fs_template, "/\n",
      "  This may take up to a minute on the first run."
    )
  }

  NULL  # Python will call fetch_surf_fsaverage
}


# ── main exported function ────────────────────────────────────────────────────

#' @title Plot vertex-wise maps on brain surfaces
#'
#' @description
#' Renders left and/or right hemisphere vertex-wise scalar maps on a standard
#' fsaverage surface using Python/nilearn. Fully headless — no XQuartz, rgl,
#' or display server required.
#'
#' Surface meshes are resolved in order:
#' \enumerate{
#'   \item \env{FREESURFER_HOME}/subjects/\code{fs_template}/surf/ — used
#'         directly if the directory exists (no download needed).
#'   \item nilearn automatic download, cached permanently in
#'         \file{~/nilearn_data/} (download only happens once per template).
#' }
#'
#' \strong{Interactive} (\code{output_file = NULL}): generates a
#' self-contained WebGL HTML file opened in the IDE Viewer or default
#' browser.\cr
#' \strong{Static} (\code{output_file} supplied): saves a tiled PNG via
#' matplotlib.
#'
#' @param lh Numeric vector, \code{.mgh}/\code{.gii} file path, or
#'   \code{NULL}.
#' @param rh Numeric vector, \code{.mgh}/\code{.gii} file path, or
#'   \code{NULL}. At least one of \code{lh}/\code{rh} must be supplied.
#' @param surface Surface mesh: \code{"inflated"} (default), \code{"pial"},
#'   or \code{"white"}.
#' @param cmap Matplotlib colormap name. Default \code{"RdBu_r"}.
#' @param bg_map Background shading: \code{"sulc"} (default),
#'   \code{"curv"}, or \code{"none"}.
#' @param vmin,vmax Numeric colour limits, or \code{NULL} for automatic
#'   symmetric scaling.
#' @param threshold Absolute-value masking threshold, or \code{NULL}
#'   (no masking).
#' @param views Character vector of camera angles — any subset of
#'   \code{"lateral"}, \code{"medial"}, \code{"dorsal"}, \code{"ventral"},
#'   \code{"anterior"}, \code{"posterior"}. Default
#'   \code{c("lateral", "medial")}.
#' @param colorbar Logical. Draw a shared colour bar? Default \code{TRUE}.
#' @param colorbar_label Character label for the colour bar axis, or
#'   \code{NULL}.
#' @param title Character figure title, or \code{NULL}.
#' @param output_file Path ending in \code{.png} for static export, or
#'   \code{NULL} for interactive HTML mode.
#' @param dpi Integer output resolution for static PNG. Default \code{150L}.
#' @param fs_template fsaverage template: \code{"fsaverage3"} through
#'   \code{"fsaverage"}. Default \code{"fsaverage"}.
#'
#' @return Invisibly: \code{output_file} path (static) or the temp HTML
#'   file path (interactive). Called primarily for the side-effect.
#'
#' @examples
#' \dontrun{
#' # Interactive — opens in Viewer / browser
#' plot_vw_surf(lh = my_coef, fs_template = "fsaverage5")
#'
#' # Static PNG, 4 views, 300 dpi
#' plot_vw_surf(
#'   lh          = lh_coef,
#'   rh          = rh_coef,
#'   cmap        = "RdBu_r",
#'   threshold   = 0.05,
#'   views       = c("lateral", "medial", "dorsal", "ventral"),
#'   title       = "Effect of age on cortical thickness",
#'   output_file = "figures/age_thickness.png",
#'   dpi         = 300L
#' )
#' }
#'
#' @export
plot_vw_surf <- function(
    lh              = NULL,
    rh              = NULL,
    surface         = c("inflated", "pial", "white"),
    cmap            = "RdBu_r",
    bg_map          = c("sulc", "curv", "none"),
    vmin            = NULL,
    vmax            = NULL,
    threshold       = NULL,
    views           = c("lateral", "medial"),
    colorbar        = TRUE,
    colorbar_label  = NULL,
    title           = NULL,
    output_file     = NULL,
    dpi             = 150L,
    fs_template     = "fsaverage"
) {

  ## ── input validation ───────────────────────────────────────────────────────
  if (is.null(lh) && is.null(rh))
    stop("At least one of `lh` or `rh` must be supplied.", call. = FALSE)

  if (!is.null(output_file) &&
      !grepl("\\.png$", output_file, ignore.case = TRUE))
    stop("`output_file` must end in '.png'.", call. = FALSE)

  surface <- match.arg(surface)
  bg_map  <- match.arg(bg_map)

  .valid_views <- c("lateral", "medial", "dorsal",
                    "ventral", "anterior", "posterior")
  bad_views <- setdiff(views, .valid_views)
  if (length(bad_views))
    stop("Invalid view(s): ", paste(bad_views, collapse = ", "),
         "\n  Choose from: ", paste(.valid_views, collapse = ", "),
         call. = FALSE)

  ## ── vertex-count validation ────────────────────────────────────────────────
  n_vert <- count_vertices(fs_template)

  .check_hemi <- function(x, label) {
    if (is.null(x) || is.character(x)) return(x)  # file path: skip check
    x <- as.numeric(x)
    if (length(x) != n_vert)
      stop(label, " vector length (", length(x), ") does not match ",
           fs_template, " vertex count (", n_vert, ").", call. = FALSE)
    x
  }
  lh <- .check_hemi(lh, "lh")
  rh <- .check_hemi(rh, "rh")

  ## ── mesh resolution + progress messaging ──────────────────────────────────
  fs_home_arg <- .resolve_mesh(fs_template)

  ## ── initialise Python renderer (once per session) ─────────────────────────
  reticulate::py_require(c("nilearn", "matplotlib", "numpy"))
  .vw_surf_init_py()

  ## ── coerce lh/rh for Python ───────────────────────────────────────────────
  .prep <- function(x)
    if (is.character(x)) x else if (!is.null(x)) as.numeric(x) else NULL

  common <- list(
    lh             = .prep(lh),
    rh             = .prep(rh),
    surface        = surface,
    bg_map_type    = bg_map,
    cmap           = cmap,
    vmin           = vmin,
    vmax           = vmax,
    threshold      = threshold,
    views          = as.list(views),
    colorbar       = colorbar,
    colorbar_label = colorbar_label,
    title          = title,
    fs_template    = fs_template,
    fs_home        = fs_home_arg        # NULL → Python None → nilearn download
  )

  ## ── dispatch ───────────────────────────────────────────────────────────────
  if (!is.null(output_file)) {

    out_dir <- dirname(normalizePath(output_file, mustWork = FALSE))
    if (!dir.exists(out_dir))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    do.call(reticulate::py$vw_surf_static,
            c(common, list(output_file = output_file,
                           dpi         = as.integer(dpi))))
    message("\u2714 Brain map saved to: ", output_file)
    return(invisible(output_file))

  } else {

    reticulate::py_require("plotly")
    tmp_html <- tempfile(fileext = ".html")

    do.call(reticulate::py$vw_surf_interactive,
            c(common, list(output_html = tmp_html)))

    viewer <- getOption("viewer", utils::browseURL)
    viewer(tmp_html)
    message("\u2714 Interactive brain map opened")
    return(invisible(tmp_html))
  }
}
