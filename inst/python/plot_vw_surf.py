import io
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import matplotlib.image as mpimg
from nilearn.datasets import fetch_surf_fsaverage
from nilearn.plotting import plot_surf 
from nilearn.surface import load_surf_data


# ── helpers ───────────────────────────────────────────────────────────────────

def _to_array(x):
    if x is None:
        return None
    if isinstance(x, str):
        arr = np.ravel(load_surf_data(x)).astype(np.float32)
    else:
        arr = np.asarray(x, dtype=np.float32).ravel()
    # R NA arrives as NaN in float32 — ensure they are proper IEEE NaN
    # (no-op for normal data, defensive for edge cases)
    arr[~np.isfinite(arr)] = np.nan
    return arr


def _auto_limits(arrays):
    vals = np.concatenate(
        [a[np.isfinite(a)] for a in arrays if a is not None]
    )
    mx = float(np.max(np.abs(vals))) if vals.size else 1.0
    return -mx, mx

def _match_mesh_size(data, target_len):
    """Helper to crop or pad an array to match the mesh vertex count."""

    if len(data) > target_len:
        return data[:target_len]

    elif len(data) < target_len:
        padded = np.full(target_len, np.nan)
        padded[:len(data)] = data
        return padded
        
    return data

def _fetch_surf_from_dir(hemi, surf_dir, surface, bg_map_type):
    
    mesh = os.path.join(surf_dir, f"{hemi}.{surface}")
    back = os.path.join(surf_dir, f"{hemi}.{bg_map_type}") if bg_map_type != "none" else None

    return mesh, back

def _fetch_surf_from_cache(hemi, fsavg, surface, bg_map_type):

    surf_key = "infl" if surface == "inflated" else surface

    mesh = getattr(fsavg, f"{surf_key}_{hemi}")
    back = getattr(fsavg, f"{bg_map_type}_{hemi}") if bg_map_type != "none" else None

    return mesh, back

def _get_mesh_paths(lh, rh, fs_template, surface, bg_map_type, fs_home):

    # Hardcoded vertex counts to avoid loading meshes when using standard templates
    FS_TEMPLATE_VERTICES = {
        "fsaverage": 163842,
        "fsaverage6": 40962,
        "fsaverage5": 10242,
        "fsaverage4": 2562,
        "fsaverage3": 642,
    }

    n_vertices = FS_TEMPLATE_VERTICES.get(fs_template)

    hemis  = []

    if fs_home:

        surf_dir = os.path.join(fs_home, "subjects", fs_template, "surf")
        
        if os.path.isdir(surf_dir):
            
            if lh is not None:
                lh = _match_mesh_size(lh, n_vertices)
                mlh, blh = _fetch_surf_from_dir('lh', surf_dir, surface, bg_map_type)
                hemis.append(("left", lh, mlh, blh))
        
            if rh is not None:
                rh = _match_mesh_size(rh, n_vertices)
                mrh, brh = _fetch_surf_from_dir('rh', surf_dir, surface, bg_map_type)
                hemis.append(("right", rh, mrh, brh))
            
            return hemis
    
    # Fall back to nilearn download + permanent cache
    fsavg = fetch_surf_fsaverage(mesh = fs_template)

    if lh is not None: 
        lh = _match_mesh_size(lh, n_vertices)
        mlh, blh = _fetch_surf_from_cache('left', fsavg, surface, bg_map_type)
        hemis.append(("left", lh, mlh, blh))

    if rh is not None:
        rh = _match_mesh_size(rh, n_vertices)
        mrh, brh = _fetch_surf_from_cache('right', fsavg, surface, bg_map_type)
        hemis.append(("right", rh, mrh, brh))

    return hemis


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


# --- Static PNG renderer ----------------------------------------------------

def vw_surf_static(lh, rh, surface, bg_map_type, cmap, vmin, vmax,
                   threshold, views, colorbar, colorbar_label,
                   title, output_file, dpi, fs_template, fs_home):

    lh, rh = _to_array(lh), _to_array(rh)
    
    if vmin is None or vmax is None:
        vmin, vmax = _auto_limits([lh, rh])

    hemis = _get_mesh_paths(lh, rh, fs_template, surface, bg_map_type, fs_home)

    thresh = float(threshold) if threshold is not None else None

    panels = []
    for hname, data, mesh, bg in hemis:
        row = []
        for view in views:
            fp = plot_surf(
                surf_mesh = mesh, 
                surf_map = data,
                hemi = hname,
                view = view, cmap = cmap, vmin = vmin, vmax = vmax,
                threshold = thresh,
                symmetric_cmap = False,
                avg_method = 'median',
                bg_map = bg, bg_on_data = False, alpha = 1.0,
                colorbar = False,
                engine = "matplotlib",
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


# --- Interactive HTML renderer ----------------------------------------------

def vw_surf_interactive(lh, rh, surface, bg_map_type, cmap, vmin, vmax,
                        threshold, views, colorbar, colorbar_label,
                        title, output_html, fs_template, fs_home):

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
            fp  = plot_surf(
                surf_mesh = mesh,
                surf_map  = data,
                hemi      = hname,
                view      = view,
                cmap      = cmap,
                vmin      = vmin,
                vmax      = vmax,
                threshold = thresh,
                symmetric_cmap = False,
                avg_method = 'median',
                bg_map    = bg,
                bg_on_data = False,
                colorbar  = True,
                engine    = "plotly",
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
