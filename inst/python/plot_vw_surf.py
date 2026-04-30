import os
import numpy as np
import matplotlib
matplotlib.use("Agg", force=True) # Force headless backend to prevent crashes
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from nilearn.datasets import fetch_surf_fsaverage
from nilearn.plotting import plot_surf 
from nilearn.surface import load_surf_data
from plotly.subplots import make_subplots
import warnings

# ── Constants ─────────────────────────────────────────────────────────────────

# Hardcoded vertex counts to avoid loading meshes when using standard templates
_FS_TEMPLATE_COUNTS = {
    "fsaverage":  163842,
    "fsaverage6": 40962,
    "fsaverage5": 10242,
    "fsaverage4": 2562,
    "fsaverage3": 642,
}

_CMAP_REG_NAME = "_verywise_cmap"

_DEFAULT_VIEWS  = ["lateral", "dorsal", "anterior", 
                   "medial", "ventral", "posterior"]
# ^ order matters: each entry becomes one column (bilateral hemi-views expand to 2).
# Default produces the 2×4 layout:
#   Left laternal | Right lateral | Dorsal  | Anterior
#   Left  medial  | Right  medial | Ventral | Posterior

# _NILEARN_CAMERAS = {
#     "left":      {"eye": {"x": -1.5,"y": 0,   "z": 0  }, "up": {"x": 0,  "y": 0, "z": 1}, "center": {"x": 0, "y": 0, "z": 0}},
#     "right":     {"eye": {"x": 1.5, "y": 0,   "z": 0  }, "up": {"x": 0,  "y": 0, "z": 1}, "center": {"x": 0, "y": 0, "z": 0}},
#     "dorsal":    {"eye": {"x": 0,   "y": 0,   "z": 1.5 },"up": {"x": -1, "y": 0, "z": 0}, "center": {"x": 0, "y": 0, "z": 0}},
#     "ventral":   {"eye": {"x": 0,   "y": 0,   "z": -1.5},"up": {"x": 1,  "y": 0, "z": 0}, "center": {"x": 0, "y": 0, "z": 0}},
#     "anterior":  {"eye": {"x": 0,   "y": 1.5, "z": 0  }, "up": {"x": 0,  "y": 0, "z": 1}, "center": {"x": 0, "y": 0, "z": 0}},
#     "posterior": {"eye": {"x": 0,   "y": -1.5,"z": 0  }, "up": {"x": 0,  "y": 0, "z": 1}, "center": {"x": 0, "y": 0, "z": 0}},
# }

_VIEW_DISTANCES = { # Smaller numbers ->  closer
    "lateral":   1.18,
    "medial":    1.00,  # concave surface, needs to be slightly closer
    "dorsal":    1.15,
    "ventral":   1.30,
    "anterior":  1.15,  # narrow axis, zoom in more
    "posterior": 1.15,
}

# ── Array helpers ─────────────────────────────────────────────────────────────

def _to_array(x):
    if x is None:
        return None

    if isinstance(x, str):
        arr = np.ravel(load_surf_data(x)).astype(np.float32)
    else:
        arr = np.asarray(x, dtype=np.float32).ravel()

    # Ensure they are proper IEEE NaN (R NA -> NaN in float32) 
    arr[~np.isfinite(arr)] = np.nan

    return arr

def _match_mesh_size(data, target_len):
    """Warn and crop/pad data array to match mesh vertex count."""

    if len(data) > target_len:
        return data[:target_len]

    elif len(data) < target_len:
        raise ValueError(f"'{target_len}' vertices is higher than the data resolution.")
        # padded = np.full(target_len, np.nan)
        # padded[:len(data)] = data
        # return padded
        
    return data


def _get_meshes(lh, rh, fs_template, surface, bg_map_type, fs_home, load_dk=False):
    """
    Returns list of (hemi_name, data, mesh_path, bg_path) for each present hemi.
    Handles both local FreeSurfer directory and nilearn cache.
    """

    n_vertices = _FS_TEMPLATE_COUNTS.get(fs_template)

    if n_vertices is None:
        raise ValueError(f"Unknown fs_template '{fs_template}'.")

    use_local = False
    if fs_home:
        surf_dir = os.path.join(fs_home, "subjects", fs_template, "surf")
        use_local = os.path.isdir(surf_dir)

        if load_dk:
            import nibabel as nib
        
    if not use_local:
        fsavg = fetch_surf_fsaverage(mesh=fs_template)
        surf_key = "infl" if surface == "inflated" else surface

        if load_dk:
            from nilearn.datasets import fetch_atlas_surf_destrieux
            destrieux = fetch_atlas_surf_destrieux()

    hemis  = []

    for hemi_name, data in [("left", lh), ("right", rh)]:

        if data is None:
            continue

        roi = None  # ← reset per hemi

        data = _match_mesh_size(data, n_vertices)
        short = hemi_name[0] + "h"   # "lh" / "rh"

        if use_local:
            mesh = os.path.join(surf_dir, f"{short}.{surface}")
            back = os.path.join(surf_dir, f"{short}.{bg_map_type}") if bg_map_type != "none" else None

            for p in filter(None, [mesh, back]):
                if not os.path.exists(p):
                    raise FileNotFoundError(f"Surface file not found: {p}")
            
            if load_dk:
                annot = os.path.join(surf_dir, "..", "label", f"{short}.aparc.annot")
                labels, _, names = nib.freesurfer.read_annot(annot)
                names = [n.decode() if isinstance(n, bytes) else n for n in names]
                roi = (labels, names)

        else:
            mesh = getattr(fsavg, f"{surf_key}_{hemi_name}")
            back = getattr(fsavg, f"{bg_map_type}_{hemi_name}") if bg_map_type != "none" else None

            if load_dk:
                labels = np.asarray(destrieux[f"map_{hemi_name}"], dtype=int)
                names = [n.decode() if isinstance(n, bytes) else str(n) for n in destrieux["labels"]]
                roi = (labels, names)
        
        back = _lighten_bg(back, darkness=0.3)

        hemis.append((hemi_name, data, mesh, back, roi))
    
    return hemis

def _build_surface_images(lh, rh, fs_template, surface, bg_map_type, fs_home, load_dk=False):
    """
    Parse mesh files ONCE, return (surf_img, bg_img, hemis, roi_dict).
    surf_img and bg_img are SurfaceImage objects carrying their own mesh.
    """
    from nilearn.surface.surface import PolyMesh
    from nilearn.surface import SurfaceImage
    import nibabel as nib

    n_vertices = _FS_TEMPLATE_COUNTS.get(fs_template)
    if n_vertices is None:
        raise ValueError(f"Unknown fs_template '{fs_template}'.")

    use_local = False
    if fs_home:
        surf_dir = os.path.join(fs_home, "subjects", fs_template, "surf")
        use_local = os.path.isdir(surf_dir)

    if not use_local:
        fsavg = fetch_surf_fsaverage(mesh=fs_template)
        surf_key = "infl" if surface == "inflated" else surface
        if load_dk:
            from nilearn.datasets import fetch_atlas_surf_destrieux
            destrieux = fetch_atlas_surf_destrieux()

    mesh_parts, data_parts, bg_parts, roi_dict = {}, {}, {}, {}
    hemi_present = []

    for hemi_name, data in [("left", lh), ("right", rh)]:
        if data is None:
            continue

        data = _match_mesh_size(data, n_vertices)
        short = hemi_name[0] + "h"

        if use_local:
            mesh_path = os.path.join(surf_dir, f"{short}.{surface}")
            bg_path   = os.path.join(surf_dir, f"{short}.{bg_map_type}") if bg_map_type != "none" else None
            for p in filter(None, [mesh_path, bg_path]):
                if not os.path.exists(p):
                    raise FileNotFoundError(f"Surface file not found: {p}")
            if load_dk:
                annot = os.path.join(surf_dir, "..", "label", f"{short}.aparc.annot")
                labels, _, names = nib.freesurfer.read_annot(annot)
                names = [n.decode() if isinstance(n, bytes) else n for n in names]
                roi_dict[hemi_name] = (labels, names)
        else:
            mesh_path = getattr(fsavg, f"{surf_key}_{hemi_name}")
            bg_path   = getattr(fsavg, f"{bg_map_type}_{hemi_name}") if bg_map_type != "none" else None
            if load_dk:
                labels = np.asarray(destrieux[f"map_{hemi_name}"], dtype=int)
                names  = [n.decode() if isinstance(n, bytes) else str(n) for n in destrieux["labels"]]
                roi_dict[hemi_name] = (labels, names)

        mesh_parts[hemi_name] = mesh_path
        data_parts[hemi_name] = data

        hemi_present.append(hemi_name)
        
        # Load and lighten background map once per hemisphere
        if bg_path is not None:
            bg_parts[hemi_name] = _lighten_bg(bg_path, darkness=0.3)

    # PolyMesh parses all geometry files exactly once here
    poly_mesh = PolyMesh(**mesh_parts)
    surf_img  = SurfaceImage(mesh=poly_mesh, data=data_parts)
    bg_img    = SurfaceImage(mesh=poly_mesh, data=bg_parts) if bg_parts else None

    hemis = 'both' if len(hemi_present) > 1 else hemi_present[0]

    return surf_img, bg_img, hemis, roi_dict

# ── Colormap helpers ──────────────────────────────────────────────────────────

def _wrap_cmap(name_or_obj):
    """Return a copy of the colormap with NaN rendered as fully transparent."""
    cmap = matplotlib.colormaps[name_or_obj] if isinstance(name_or_obj, str) \
           else name_or_obj
    c = cmap.copy()
    c.set_bad(alpha=0.0)
    return c


def _detect_discrete(vals, max_discrete=20):
    """True if all finite values look like integer labels (1-20 unique values)."""
    if vals.size == 0:
        return False
    unique = np.unique(vals)

    # Short-circuit before the allclose call (cheap check first)
    if len(unique) > max_discrete:
        return False

    return bool(np.allclose(vals, np.round(vals)))


def _make_diverging_cmap():
    """viridis (negative) → hot_r (positive) diverging colormap."""
    hot_r = matplotlib.colormaps['hot_r']
    viridis = matplotlib.colormaps['viridis']
    colors = [viridis(0.0), viridis(0.5), viridis(1.0),
              hot_r(0.5), hot_r(0.75), hot_r(1.0)]
    nodes  = [0.0, 0.25, 0.5, 0.6, 0.75, 1.0]
    return mcolors.LinearSegmentedColormap.from_list(
        "hot_r_viridis", list(zip(nodes, colors)))


def _make_continuous_cmap(vmin, vmax):
    """Pick colormap direction based on sign of the data range."""
    if vmin >= 0:       # all positive → viridis_r (dark=high)
        return _wrap_cmap('viridis_r')
    elif vmax <= 0:     # all negative → viridis (dark=low/most negative)
        return _wrap_cmap('viridis')
    else:               # diverging
        return _wrap_cmap(_make_diverging_cmap())


def _make_discrete_cmap(unique_vals):
    """One distinct color per integer label, NaN transparent."""
    n = len(unique_vals)
    palette = matplotlib.colormaps['tab20' if n > 10 else 'tab10']
    colors  = [palette(i / max(n - 1, 1)) for i in range(n)]
    return _wrap_cmap(mcolors.ListedColormap(colors, name='discrete'))


def _resolve_cmap_and_limits(data_arrays, cmap=None, vmin=None, vmax=None,
                             discrete=None):
    """
    Auto-detect continuous vs. discrete data and return
    (cmap_obj, plot_vmin, plot_vmax, norm).

    Parameters
    ----------
    cmap     : str | list of colors | None  (None = auto)
    discrete : True | False | None          (None = auto-detect)
    """

    vals = np.concatenate(
        [a[np.isfinite(a)] for a in data_arrays if a is not None])

    if vals.size == 0:
        empty_norm = mcolors.Normalize(0.0, 1.0)
        return _wrap_cmap('viridis'), 0.0, 1.0, empty_norm

    data_min = float(np.min(vals))
    data_max = float(np.max(vals))

    # Auto-detect discrete
    if discrete is None:
        discrete = _detect_discrete(vals)

    # ─── Build colormap ───
    if cmap is not None:
        # User override: list of hex/name colors → ListedColormap; string → named
        cmap_obj = _wrap_cmap(
            mcolors.ListedColormap(cmap) if isinstance(cmap, list)
            else matplotlib.colormaps[cmap]
        )
    elif discrete:
        unique_vals = np.unique(np.round(vals).astype(int))
        cmap_obj = _make_discrete_cmap(unique_vals)
    else:
        # Continuous: choose based on sign of auto limits (before user overrides)
        auto_vmin = vmin if vmin is not None else data_min
        auto_vmax = vmax if vmax is not None else data_max
        cmap_obj = _make_continuous_cmap(auto_vmin, auto_vmax)

    # ─── Compute limits & norm ───
    if discrete:
        unique_vals = np.unique(np.round(vals).astype(int))
        # Half-integer boundaries → each integer maps cleanly to one color bin
        plot_vmin = float(unique_vals[0])  - 0.5 if vmin is None else vmin
        plot_vmax = float(unique_vals[-1]) + 0.5 if vmax is None else vmax
        boundaries = np.arange(unique_vals[0] - 0.5, unique_vals[-1] + 1.5)
        norm = mcolors.BoundaryNorm(boundaries=boundaries,
                                    ncolors=len(unique_vals))

    else:
        plot_vmin = vmin if vmin is not None else data_min
        plot_vmax = vmax if vmax is not None else data_max

        # Enforce symmetry for diverging only when user specified neither bound
        if data_min < 0 < data_max and vmin is None and vmax is None:
            mx = max(abs(plot_vmin), abs(plot_vmax))
            plot_vmin, plot_vmax = -mx, mx

        norm = mcolors.Normalize(vmin=plot_vmin, vmax=plot_vmax)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Overwriting the cmap")
        matplotlib.colormaps.register(cmap_obj, name=_CMAP_REG_NAME, force=True)

    return cmap_obj, plot_vmin, plot_vmax, norm

def _lighten_bg(bg_path, darkness=0.3):
    """Compress bg_map dynamic range toward mid-grey (0.5)."""
    bg = load_surf_data(bg_path) if isinstance(bg_path, str) else bg_path

    if bg is None:
        return None

    arr = np.asarray(bg, dtype=np.float32)
    # normalize to [0,1] first
    lo, hi = np.nanmin(arr), np.nanmax(arr)
    if hi > lo:
        arr = (arr - lo) / (hi - lo)
    # blend toward 0.5 (mid-grey): darkness=1 → full range, darkness=0 → flat grey
    return 0.5 + (arr - 0.5) * darkness

# ── Figure helpers ────────────────────────────────────────────────────────────

def _unwrap_plotly(fp):
    """nilearn returns PlotlySurfaceFigure; unwrap to raw plotly Figure."""
    return fp.figure if hasattr(fp, "figure") else fp


def _hover_text(i, v, roi, x, y, z):
    if not np.isfinite(v):
        return ""
    line = f"Value: <b>{v:.4f}</b><br>"
    if roi is not None:
        lbl_idx = roi[0][i]   # roi = (labels_array, names_list)
        region  = roi[1][lbl_idx] if lbl_idx >= 0 else "unknown"
        line   += f"Region: <b>{region}</b><br>"
    line += f"x: {x:.2f}  y: {y:.2f}  z: {z:.2f}"
    return line


def _panel_layout(hemis, views='all', n_cols=None):
    """
    Build the subplot panel list from a views argument.

    Parameters
    ----------
    hemis : "left", "right", or "both"
    views : list of view names, e.g. ["lateral", "dorsal"].
            "lateral"/"medial" expand to left + right columns in
             bilateral mode; all other views map to hemi="both".
    n_cols : int or None. Auto-detected if None.

    Returns
    -------
    n_rows, n_cols, panels
        panels = list of (row, col, hemi, view, title_or_None)
    """
    import math

    if (views == 'all') or (len(views) == 6):
        views = _DEFAULT_VIEWS

    bilateral = hemis == "both"
   
    # Expand each view name → one or two (hemi, view) pairs
    expanded = []
    for v in views:
        if bilateral and v not in ["dorsal", "anterior", "ventral", "posterior"]:
            expanded.append(("left",  v))
            expanded.append(("right", v))
        else:
            expanded.append((hemis, v))

    n_panels = len(expanded)

    # Auto-detect n_cols: 4 for bilateral, 3 for single-hemi (capped at actual panel count)
    if n_cols is None:
        n_cols = min(4 if bilateral else 3, n_panels)
    n_rows = math.ceil(n_panels / n_cols)

    # Title logic:
    #   left/right → show "Left"/"Right" as column header only on first appearance
    #   both/single → always show the view name (each panel is a distinct view)
    col_hemi_seen = {}

    panels = []
    for i, (hemi, view) in enumerate(expanded):
        row, col = divmod(i, n_cols)

        if hemi in ("left", "right"):
            seen = col_hemi_seen.setdefault(col, set())
            panel_title = hemi.capitalize() if hemi not in seen else None
            seen.add(hemi)
        else:
            panel_title = view.capitalize()

        panels.append((row, col, hemi, view, panel_title))

    return n_rows, n_cols, panels


def _scene_key(row, col, n_cols):
    idx = row * n_cols + col
    return "scene" if idx == 0 else f"scene{idx + 1}"


def _copy_camera(raw_layout):
    """Extract only the camera dict from a nilearn plotly layout."""
    if hasattr(raw_layout, "scene") and raw_layout.scene.camera:
        return raw_layout.scene.camera.to_plotly_json()
    return {}

def _normalize_eye(camera_dict, view="lateral", distances=None):
    """
    Keep nilearn's view direction, apply per-view zoom distance,
    and reset center to origin (fixes off-center on inflated surface).
    """
    if distances is None:
        distances = _VIEW_DISTANCES

    distance = distances.get(view, 1.5)

    if not camera_dict:
        return {}

    eye = camera_dict.get("eye", {})
    x, y, z = eye.get("x", 0), eye.get("y", 0), eye.get("z", distance)
    mag = math.sqrt(x**2 + y**2 + z**2)
    if mag == 0:
        return camera_dict

    s = distance / mag
    return {
        **camera_dict,
        "eye":    {"x": x*s, "y": y*s, "z": z*s},
        "center": {"x": 0.0, "y": 0.0, "z": 0.0},  # ← critical for inflated
    }

# ──────────────────────────────────────────────────────────────────────────────
# ── Static PNG renderer ───────────────────────────────────────────────────────

# This static plotting, based on matplotlib. It works but it is much slower than 
# the kaleido based screenshots below 
# def vw_surf_static(lh, rh, surface, bg_map_type, vmin, vmax,
#                    threshold, views, colorbar, colorbar_label,
#                    title, output_file, dpi, fs_template, fs_home, cmap = _CMAP_REG_NAME):

#     lh, rh = _to_array(lh), _to_array(rh)
    
#     cmap_obj, vmin, vmax, norm = _resolve_cmap_and_limits(
#         [lh, rh], cmap=cmap, vmin=vmin, vmax=vmax)

#     # hemis = _get_meshes(lh, rh, fs_template, surface, bg_map_type, fs_home)
#     # Build the SurfaceImage objects (parses meshes only once)
#     surf_img, bg_img, hemis_present, _ = _build_surface_images(
#         lh, rh, fs_template, surface, bg_map_type, fs_home)

#     thresh = float(threshold) if threshold is not None else None

#     bilateral = hemis_present == {"left", "right"}

#     # Define panel layouts based on provided hemispheres
#     if bilateral:
#         n_rows, n_cols = 2, 4
#         panels = [
#             (0, 0, "left",  "lateral",   "Left"),
#             (0, 1, "right", "lateral",   "Right"),
#             (0, 2, "left",  "dorsal",    "Dorsal"),
#             (0, 3, "left",  "anterior",  "Anterior"),
#             (1, 0, "left",  "medial",    None),
#             (1, 1, "right", "medial",    None),
#             (1, 2, "left",  "ventral",   "Ventral"),
#             (1, 3, "left",  "posterior", "Posterior"),
#         ]
#     else:
#         n_rows, n_cols = 2, 3
#         h = next(iter(hemis_present))
#         panels = [
#             (0, 0, h, "lateral",   f"{h.capitalize()} hemisphere"),
#             (0, 1, h, "dorsal",    "Dorsal"),
#             (0, 2, h, "anterior",  "Anterior"),
#             (1, 0, h, "medial",    None),
#             (1, 1, h, "ventral",   "Ventral"),
#             (1, 2, h, "posterior", "Posterior"),
#         ]

#     n_grid = n_cols + (1 if colorbar else 0)
#     fig = plt.figure(
#         figsize=(n_cols * 3.5 + (0.5 if colorbar else 0.0),
#                  n_rows * 2.8 + (0.6 if title else 0.15)),
#         facecolor="white", dpi=dpi,
#     )
#     gs = gridspec.GridSpec(
#         n_rows, n_grid,
#         width_ratios=[1.0] * n_cols + ([0.04] if colorbar else []),
#         hspace=0.01, wspace=0.01,
#         left=0.02, right=0.92 if colorbar else 0.99,
#         top=0.92 if title else 0.98, bottom=0.02,
#     )

#     for row, col, hemi, view, col_title in panels:
#         ax = fig.add_subplot(gs[row, col], projection="3d")
        
#         # We pass surf_mesh=None because surf_map is a SurfaceImage
#         plot_surf(
#             surf_mesh=None, surf_map=surf_img, hemi=hemi, view=view,
#             cmap=cmap_obj, vmin=vmin, vmax=vmax, threshold=thresh,
#             symmetric_cmap=None, avg_method='mean',
#             bg_map=bg_img, bg_on_data=False, alpha=1.0,
#             colorbar=False, engine="matplotlib",
#             axes=ax, figure=fig,
#         )
        
#         ax.axis("off")
#         if col_title:
#             ax.set_title(col_title, fontsize=10, pad=4)

#     if colorbar:
#         cax = fig.add_subplot(gs[:, n_cols])
#         sm = matplotlib.cm.ScalarMappable(cmap=cmap_obj, norm=norm)
#         sm.set_array([])
#         cb = fig.colorbar(sm, cax=cax)
#         if isinstance(norm, mcolors.BoundaryNorm):
#             ticks = np.round(norm.boundaries[:-1] + 0.5).astype(int)
#             cb.set_ticks(ticks)
#             cb.set_ticklabels([str(v) for v in ticks])
#         if colorbar_label:
#             cb.set_label(colorbar_label, fontsize=9)
#         cb.ax.tick_params(labelsize=8)

#     if title:
#         fig.suptitle(title, fontsize=13, fontweight="bold", y=0.99)

#     os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
#     fig.savefig(output_file, dpi=dpi, bbox_inches="tight", facecolor="white")
#     plt.close(fig)

#     return output_file

# ──────────────────────────────────────────────────────────────────────────────
# ── Interactive HTML renderer ─────────────────────────────────────────────────

def vw_surf_interactive(lh, rh, surface, bg_map_type, vmin, vmax,
                        threshold, views, colorbar, colorbar_label,
                        title, output_html, fs_template, fs_home, cmap = _CMAP_REG_NAME):
 
    lh, rh = _to_array(lh), _to_array(rh)

    cmap_obj, vmin, vmax, norm = _resolve_cmap_and_limits(
        [lh, rh], cmap=cmap, vmin=vmin, vmax=vmax)

    surf_img, bg_img, hemis, rois = _build_surface_images(
        lh, rh, fs_template, surface, bg_map_type, fs_home, load_dk = True)

    thresh = float(threshold) if threshold is not None else None

    start_view = "lateral"
    # TODO update rois so they work for SurfaceImage
    data_map = {"left": lh, "right": rh, "both": lh}  # hover uses lh for bilateral

    fp  = plot_surf(
        surf_mesh=None, surf_map=surf_img,
        hemi=hemis, view=start_view,
        cmap=cmap, vmin=vmin, vmax=vmax,
        threshold=thresh, symmetric_cmap=None, avg_method=None,
        bg_map=bg_img, bg_on_data=True,
        colorbar=False, engine="plotly",
    )
    raw = _unwrap_plotly(fp)

    # hover text — use lh roi/data for left traces, rh for right
    roi  = rois.get("left" if hemi != "right" else "right", ({}, []))
    data = data_map[hemi]
    xs   = np.asarray(raw.data[0].x)
    ys   = np.asarray(raw.data[0].y)
    zs   = np.asarray(raw.data[0].z)

    text_vals = [
        _hover_text(i, v, roi, x, y, z)
        for i, (v, x, y, z) in enumerate(zip(data, xs, ys, zs))
    ]

    fig = go.Figure()
    for trace in raw.data:
        trace.update(
            text=text_vals,
            hovertemplate="%{text}<extra></extra>",
            showscale=False,
        )
        fig.add_trace(trace)

    # colorbar on stat trace
    if colorbar:
        fig.data[-1].update(
            showscale=True,
            colorscale=_plotly_colorscale(cmap_obj),
            cmin=vmin, cmax=vmax,
            colorbar=dict(
                title=dict(text=colorbar_label or "", side="right"),
                tickformat=".2g",
            ),
        )
    
    cam = _build_camera(hemi, start_view, surface=surface,
                        custom_presets=camera_presets)
    fig.update_layout(
        title_text=title or "", title_font_size=14,
        paper_bgcolor="white",
        hovermode="closest",
        width=600, height=550,
        scene=dict(
            camera=cam,
            aspectmode="cube",
            xaxis=dict(visible=False, showgrid=False, zeroline=False),
            yaxis=dict(visible=False, showgrid=False, zeroline=False),
            zaxis=dict(visible=False, showgrid=False, zeroline=False),
            bgcolor="white",
        ),
    )

    os.makedirs(os.path.dirname(os.path.abspath(output_html)), exist_ok=True)
    fig.write_html(output_html, include_plotlyjs="cdn",
                   full_html=True, auto_open=False)

    return output_html

def _mpl_to_plotly_colorscale(cmap, n=256):
    """Convert a matplotlib colormap to a Plotly colorscale list."""
    return [
        [i / (n - 1), f"rgba({int(r*255)},{int(g*255)},{int(b*255)},{a:.3f})"]
        for i, (r, g, b, a) in enumerate(cmap(np.linspace(0, 1, n)))
    ]

import math

# ─────────────────────────────────────────────────────────────────────────────
# Camera presets — the single source of truth for view angles and zoom.
#
# Tuning guide:
#   'eye'      — unit direction vector toward the camera. Changing sign flips
#                the view to the opposite side. Magnitude is normalised so only
#                direction matters here; zoom is set by 'distance'.
#   'distance' — eye vector is scaled to this length. Smaller = more zoomed in.
#                Start with 1.5 (lateral) and work inward for narrow-axis views.
#   'up'       — which world direction points "up" in the panel.
#                For lateral/medial/anterior/posterior: (0,0,1) = superior up.
#                For dorsal/ventral you must rotate: (-1,0,0) / (1,0,0).
#   'center'   — always (0,0,0). Changing this breaks rotation and centering.
# ─────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
# Surface-keyed camera presets.
#
# Structure: _CAMERA_PRESETS[surface][(hemi, view)]
#   surface  : "pial" | "inflated" | "white"
#   eye      : unit direction; only direction matters, magnitude is normalised
#   up       : which world axis points "up" in the panel
#   distance : eye-to-origin distance. smaller = more zoomed in.
#   center   : always omitted here — _build_camera hardcodes it to (0,0,0)
#
# For views that are hemi-agnostic (dorsal/ventral/anterior/posterior),
# "left"/"right" entries simply copy "both" — add them only if single-hemi
# use requires a different angle.
# ─────────────────────────────────────────────────────────────────────────────

def _bilateral_views(d):
    """Expand ("both", view) entries to also cover ("left", view) / ("right", view)."""
    extra = {}
    for (hemi, view), v in d.items():
        if hemi == "both":
            extra.setdefault(("left",  view), v)
            extra.setdefault(("right", view), v)
    return {**extra, **d}   # explicit entries win over auto-expanded ones


_CAMERA_PRESETS = {

    "pial": _bilateral_views({
        ("left",  "lateral"):   dict(eye=(-1.0, 0.0,  0.0), up=(0, 0, 1), distance=1.20),
        ("right", "lateral"):   dict(eye=( 1.0, 0.0,  0.0), up=(0, 0, 1), distance=1.20),
        ("left",  "medial"):    dict(eye=( 1.0, 0.0,  0.0), up=(0, 0, 1), distance=1.05),
        ("right", "medial"):    dict(eye=(-1.0, 0.0,  0.0), up=(0, 0, 1), distance=1.05),
        ("both",  "dorsal"):    dict(eye=( 0.0, 0.0,  1.0), up=(-1, 0, 0),distance=1.15),
        ("both",  "ventral"):   dict(eye=( 0.0, 0.0, -1.0), up=(1, 0, 0), distance=1.30),
        ("both",  "anterior"):  dict(eye=( 0.0, 1.0,  0.0), up=(0, 0, 1), distance=1.15),
        ("both",  "posterior"): dict(eye=( 0.0,-1.0,  0.0), up=(0, 0, 1), distance=1.15),
    }),

    "inflated": _bilateral_views({
        ("left",  "lateral"):   dict(eye=(-1.0, 0.0,  0.0), up=(0, 0, 1), distance=1.10),
        ("right", "lateral"):   dict(eye=( 1.0, 0.0,  0.0), up=(0, 0, 1), distance=1.10),
        ("left",  "medial"):    dict(eye=( 1.0, 0.0,  0.0), up=(0, 0, 1), distance=1.20),
        ("right", "medial"):    dict(eye=(-1.0, 0.0,  0.0), up=(0, 0, 1), distance=1.20),
        ("both",  "dorsal"):    dict(eye=(0.11, 0.0, 1.25), up=(-1, 0, 0),center=(0.15, 0, 0), distance=1.15),
        ("both",  "ventral"):   dict(eye=(0.25, 0.0, -1.3), up=(1, 0, 0), center=(0.18, 0, 0), distance=1.30),
        ("both",  "anterior"):  dict(eye=(0.16, 1.1, -0.1), up=(0, 0.2,1),center=(0.20, 0, 0), distance=1.15),
        ("both",  "posterior"): dict(eye=(0.20,-1.1, -0.0), up=(0, 0, 1), center=(0.18, 0, 0), distance=1.15),
    }),
  
    "white": _bilateral_views({      # white sits between pial and inflated
        ("left",  "lateral"):   dict(eye=(-1.0, 0.0,  0.0), up=(0, 0, 1), distance=1.20),
        ("right", "lateral"):   dict(eye=( 1.0, 0.0,  0.0), up=(0, 0, 1), distance=1.20),
        ("left",  "medial"):    dict(eye=( 1.0, 0.0,  0.0), up=(0, 0, 1), distance=1.05),
        ("right", "medial"):    dict(eye=(-1.0, 0.0,  0.0), up=(0, 0, 1), distance=1.05),
        ("both",  "dorsal"):    dict(eye=( 0.0, 0.0,  1.0), up=(-1, 0, 0),distance=1.18),
        ("both",  "ventral"):   dict(eye=( 0.0, 0.0, -1.0), up=(1, 0, 0), distance=1.32),
        ("both",  "anterior"):  dict(eye=( 0.0, 1.0,  0.0), up=(0, 0, 1), distance=1.15),
        ("both",  "posterior"): dict(eye=( 0.0,-1.0,  0.0), up=(0, 0, 1), distance=1.15),
    }),
}


def _build_camera(hemi, view, surface, custom_presets=None):
    """
    Resolve camera dict for a (hemi, view, surface) triple.
    Lookup order:
      1. custom_presets[surface][(hemi, view)]  — per-call override
      2. _CAMERA_PRESETS[surface][(hemi, view)] — surface-specific constant
      3. _CAMERA_PRESETS["inflated"][(hemi, view)] — last-resort fallback
    center is always locked to origin regardless of source.
    """
    surface_key = surface if surface in _CAMERA_PRESETS else "inflated"

    p = None
    for source in [custom_presets, _CAMERA_PRESETS]:
        if source and surface_key in source:
            p = source[surface_key].get((hemi, view))
            if p is not None:
                break

    if p is None:
        raise ValueError(
            f"No camera preset for surface={surface!r}, hemi={hemi!r}, view={view!r}. "
            f"Add it to _CAMERA_PRESETS[{surface_key!r}]."
        )

    ex, ey, ez = p["eye"]
    ux, uy, uz = p["up"]
    d   = p["distance"]
    mag = math.sqrt(ex**2 + ey**2 + ez**2) or 1.0
    s   = d / mag

    # center: use preset value if provided, else origin
    cx, cy, cz = p.get("center", (0.0, 0.0, 0.0))

    return dict(
        eye    = dict(x=ex*s, y=ey*s, z=ez*s),
        up     = dict(x=ux,   y=uy,   z=uz),
        center = dict(x=cx,   y=cy,   z=cz),
    )

# ─────────────────────────────────────────────────────────────────────────────

os.environ.setdefault("KALEIDO_CHROME_ARGS", "--headless --disable-gpu --no-sandbox")

def vw_surf_static_plotly(
        lh, rh, surface, views, bg_map_type,
        vmin, vmax, threshold, colorbar, colorbar_label, title,
        output_file, dpi, fs_template, fs_home,
        cmap=_CMAP_REG_NAME,
        camera_presets=None,
        cell_px=400):

    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    from nilearn.surface import load_surf_mesh

    lh, rh = _to_array(lh), _to_array(rh)
    cmap_obj, vmin, vmax, norm = _resolve_cmap_and_limits(
        [lh, rh], cmap=cmap, vmin=vmin, vmax=vmax)

    surf_img, bg_img, hemis, _ = _build_surface_images(
        lh, rh, fs_template, surface, bg_map_type, fs_home)

    thresh = float(threshold) if threshold is not None else None
    n_rows, n_cols, panels = _panel_layout(hemis, views)

    # ── 1. Uniform cubic scene bbox from mesh ─────────────────────────────────
    all_coords = np.vstack([
        load_surf_mesh(p)[0] for p in surf_img.mesh.parts.values()
    ])
    mid  = (all_coords.max(axis=0) + all_coords.min(axis=0)) / 2.0
    half = (all_coords.max(axis=0) - all_coords.min(axis=0)).max() / 2.0 * 1.05
    _ax  = lambda i: dict(visible=False, showgrid=False, zeroline=False,
                          range=[float(mid[i]-half), float(mid[i]+half)])

    # ── 2. Trace cache — ≤3 plot_surf() calls ────────────────────────────────
    CANONICAL   = {"left": "lateral", "right": "lateral", "both": "dorsal"}
    trace_cache = {}
    for hemi in {h for _, _, h, _, _ in panels}:
        fp  = plot_surf(
            surf_mesh=None, surf_map=surf_img,
            hemi=hemi, view=CANONICAL[hemi],
            cmap=_CMAP_REG_NAME, vmin=vmin, vmax=vmax,
            threshold=thresh, symmetric_cmap=None, avg_method=None,
            bg_map=bg_img, bg_on_data=False,
            colorbar=False, engine="plotly",
        )
        raw = _unwrap_plotly(fp)
        for tr in raw.data:
            tr.update(showscale=False, hoverinfo="skip")
        trace_cache[hemi] = list(raw.data)

    # ── 3. Single multi-scene figure — one Kaleido call ───────────────────────
    subplot_titles = [""] * (n_rows * n_cols)
    for row, col, _, _, lbl in panels:
        if lbl:
            subplot_titles[row * n_cols + col] = lbl

    fig = make_subplots(
        rows=n_rows, cols=n_cols,
        specs=[[{"type": "scene"}] * n_cols for _ in range(n_rows)],
        subplot_titles=subplot_titles,
        horizontal_spacing=0.01, vertical_spacing=0.04,
    )

    for row, col, hemi, view, _ in panels:
        for tr in trace_cache[hemi]:
            fig.add_trace(tr, row=row+1, col=col+1)

        sk = _scene_key(row, col, n_cols)
        fig.update_layout(**{sk: dict(
            aspectmode = "cube",
            xaxis  = _ax(0), yaxis = _ax(1), zaxis = _ax(2),
            bgcolor = "white",
            camera  = _build_camera(hemi, view,
                                    surface=surface,
                                    custom_presets=camera_presets),
        )})

    # ── 4. Colorbar via dummy trace ───────────────────────────────────────────
    if colorbar:
        fig.add_trace(go.Scatter3d(
            x=[None], y=[None], z=[None], mode="markers",
            marker=dict(
                size=0.001, color=[vmin, vmax],
                colorscale=_mpl_to_plotly_colorscale(cmap_obj),
                cmin=vmin, cmax=vmax, showscale=True,
                colorbar=dict(
                    title=dict(text=colorbar_label or "", side="right"),
                    thickness=15, len=0.8, x=1.02, tickformat=".2g",
                ),
            ),
            showlegend=False,
        ), row=1, col=1)

    margin_r = 100 if colorbar else 10
    fig.update_layout(
        title_text=title or "", title_font_size=14,
        paper_bgcolor="white", plot_bgcolor="white",
        height=cell_px * n_rows,
        width =cell_px * n_cols + margin_r,
        margin=dict(l=10, r=margin_r, t=50 if title else 30, b=10),
    )

    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    fig.write_image(output_file, scale=dpi / 72.0)
    return output_file


def vw_dump_camera_preset(hemi, view, surface, fs_template, fs_home=None,
                          lh=None, rh=None, bg_map_type="sulc"):

    import tempfile, json
    import plotly.graph_objects as go
    from nilearn.surface import load_surf_mesh

    lh_arr = _to_array(lh)
    rh_arr = _to_array(rh)

    if lh_arr is None and rh_arr is None:
        from nilearn.datasets import fetch_surf_fsaverage
        fsavg    = fetch_surf_fsaverage(mesh=fs_template)
        surf_key = "infl" if surface == "inflated" else surface
        coords, _ = load_surf_mesh(getattr(fsavg, f"{surf_key}_left"))
        n      = coords.shape[0]
        lh_arr = np.zeros(n, dtype=np.float32)
        rh_arr = np.zeros(n, dtype=np.float32)

    surf_img, bg_img, hemis, _ = _build_surface_images(
        lh_arr, rh_arr, fs_template, surface, bg_map_type, fs_home)

    CANONICAL = {"left": "lateral", "right": "lateral", "both": "dorsal"}
    fp  = plot_surf(
        surf_mesh=None, surf_map=surf_img,
        hemi=hemi, view=CANONICAL[hemi],
        cmap=_CMAP_REG_NAME, colorbar=False, engine="plotly",
        bg_map=bg_img, bg_on_data=True,   # ← sulc shading visible
        alpha=0.9,
    )
    raw = _unwrap_plotly(fp)

    all_coords = np.vstack([
        load_surf_mesh(p)[0] for p in surf_img.mesh.parts.values()
    ])
    mid  = (all_coords.max(axis=0) + all_coords.min(axis=0)) / 2.0
    half = (all_coords.max(axis=0) - all_coords.min(axis=0)).max() / 2.0 * 1.05
    _ax  = lambda i: dict(visible=False, showgrid=False, zeroline=False,
                          range=[float(mid[i]-half), float(mid[i]+half)])

    try:
        cam = _build_camera(hemi, view, surface=surface)
    except (ValueError, KeyError):
        cam = {}

    fig = go.Figure(data=list(raw.data))
    fig.update_layout(
        scene=dict(
            aspectmode="cube",
            xaxis=_ax(0), yaxis=_ax(1), zaxis=_ax(2),
            bgcolor="white", camera=cam,
        ),
        width=700, height=700, paper_bgcolor="white",
        margin=dict(l=10, r=10, t=80, b=10),
        title=dict(
            text=f"surface=<b>{surface}</b>  hemi=<b>{hemi}</b>  view=<b>{view}</b>",
            font=dict(size=13),
        ),
    )

    # Inject JS that shows live camera coords as a copyable overlay
    live_js = """
        <div id="cam-box" style="
            position:fixed; top:12px; right:12px; z-index:9999;
            background:rgba(20,20,20,0.82); color:#e8e8e8;
            font:13px/1.6 monospace; padding:12px 16px; border-radius:8px;
            min-width:340px; white-space:pre; user-select:all;
            box-shadow:0 2px 12px rgba(0,0,0,0.4);">
        Rotate or pan the brain to update…
        </div>
        <div style="position:fixed;bottom:12px;right:12px;z-index:9999;
            font:11px monospace;color:#888;">
            Click box → Ctrl/Cmd+C to copy
        </div>
        <script>
        (function() {
            function fmt(v) { return (v||0).toFixed(4); }
            function update() {
                var el = document.getElementsByClassName('js-plotly-plot')[0];
                if (!el) { setTimeout(update, 300); return; }
                el.on('plotly_relayout', function(e) {
                    var cam = el._fullLayout.scene.camera;
                    if (!cam) return;
                    var eye = cam.eye, up = cam.up, ctr = cam.center || {x:0,y:0,z:0};
                    var txt =
                        ' eye: (' + fmt(eye.x) + ', ' + fmt(eye.y) + ', ' + fmt(eye.z) + ')\n' +
                        '  up: (' + fmt(up.x)  + ', ' + fmt(up.y)  + ', ' + fmt(up.z)  + ')\n' +
                        ' ctr: (' + fmt(ctr.x) + ', ' + fmt(ctr.y) + ', ' + fmt(ctr.z) + ')\n\n' +
                        '# Paste into _CAMERA_PRESETS:\n' +
                        '("' + '""" + hemi + """' + '", "' + '""" + view + """' + '"): dict(\n' +
                        '    eye=(' + fmt(eye.x) + ', ' + fmt(eye.y) + ', ' + fmt(eye.z) + '),\n' +
                        '    up=('  + fmt(up.x)  + ', ' + fmt(up.y)  + ', ' + fmt(up.z)  + '),\n' +
                        '    center=(' + fmt(ctr.x) + ', ' + fmt(ctr.y) + ', ' + fmt(ctr.z) + '),\n' +
                        '    distance=' + fmt(Math.sqrt(eye.x**2+eye.y**2+eye.z**2)) + ',\n' +
                        '),';
                    document.getElementById('cam-box').innerText = txt;
                });
            }
            window.addEventListener('load', function() { setTimeout(update, 800); });
        })();
        </script>
        """

    # Write HTML manually so we can inject the overlay
    base_html = fig.to_html(include_plotlyjs="cdn", full_html=True)
    # Inject before </body>
    html_out  = base_html.replace("</body>", live_js + "\n</body>")

    f = tempfile.NamedTemporaryFile(suffix=".html", delete=False,
                                    mode="w", encoding="utf-8")
    f.write(html_out)
    f.close()

    print(f"\nPath: {f.name}")
    print("Open with:  utils::browseURL(reticulate::py$vw_dump_camera_preset(...))")
    return f.name


