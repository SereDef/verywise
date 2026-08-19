"""
verywise surface plotting
==========================
Static (PNG/Kaleido) and interactive (HTML/Plotly) renderers for vertex-wise
FreeSurfer/fsaverage surface data, built on nilearn's plot_surf.

Architecture
------------
VeryWiseMesh          : loads/validates lh & rh arrays, resolves mesh + bg map
                         paths, builds nilearn SurfaceImage objects, and always
                         loads DK (aparc) ROI labels per hemisphere.
VeryWiseColormap       : resolves colormap, colour limits, ticks, band mask and
                         ready-to-use Plotly colorbar/density traces.
VeryWiseLayout         : resolves the (rows, cols, panels) subplot grid from a
                         `hemis` + `views` spec.
VeryWiseCameraPresets  : class-level constant dict of camera presets + lookup.
VeryWiseSurfacePlotter : shared scaffolding (scene bbox, figure unwrap) used
                         by both renderers.
VeryWiseStaticPlotter  : builds its own per-hemi trace cache, ROI contour
                         overlays (via nilearn's add_contours), Kaleido export.
VeryWiseInteractivePlotter : builds its own per-view traces with rich
                         per-vertex hover text, HTML export.

Public API
----------
vw_surf_static_plotly(...)   -> PNG file path
vw_surf_interactive(...)     -> HTML file path
vw_dump_camera_preset(...)   -> standalone debug utility (unchanged)
"""

import os
import math
import numpy as np
import matplotlib
matplotlib.use("Agg", force=True)  # Force headless backend to prevent crashes
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator

from nilearn.plotting import plot_surf
from nilearn.surface import load_surf_data, load_surf_mesh
from plotly.subplots import make_subplots

# disable orjson that complains about big-endians in the mesh when writing html
import plotly.graph_objects as go
import plotly.io as pio
pio.json.config.default_engine = 'json'

# ═════════════════════════════════════════════════════════════════════════════
# Small helpers (pure, stateless, reused across classes) 
# ═════════════════════════════════════════════════════════════════════════════

def _to_array(x):
    """Coerce a path / array-like into a float32 1-D array with proper NaNs."""
    if x is None:
        return None
    if isinstance(x, str):
        arr = np.ravel(load_surf_data(x)).astype(np.float32)
    else:
        arr = np.asarray(x, dtype=np.float32).ravel()
    arr[~np.isfinite(arr)] = np.nan
    return arr

def _threshold(arr, mask):
    return np.where(np.asarray(mask, dtype=bool), arr, np.nan) \
        if mask is not None and arr is not None else arr


# ═════════════════════════════════════════════════════════════════════════════
# Data + geometry + ROI labels
# ═════════════════════════════════════════════════════════════════════════════

class VeryWiseMesh:
    """
    Resolves vertex data, mesh/background paths, builds nilearn SurfaceImage
    objects, and loads DK (Desikan-Killiany / aparc) ROI labels per hemisphere
    -- used for hover text (interactive) and contour overlays (static).

    Parameters
    ----------
    lh, rh : array-like, path, or None
        Vertex-wise data for each hemisphere.
    fs_template : str
        One of the keys in the templates supported by verywise.
    surface : str
        "pial" or "inflated"
    bg_map_type : str
        e.g. "sulc", "curv", or "none".
    fs_home : str or None
        FreeSurfer home dir; if it contains subjects/<fs_template>/surf, local
        files are used instead of the nilearn fsaverage cache.

    Attributes
    ----------
    lh, rh : np.ndarray or None
        Validated/cropped data arrays.
    hemis : "left", "right", or "both"
    surf_img, bg_img : nilearn SurfaceImage
    mesh_paths : dict hemi -> mesh path (local file or nilearn cache path)
    roi : dict hemi -> (labels_array, names_list)
    """

    def __init__(self, lh, rh, fs_template, surface, bg_map_type, fs_home=None):
        self.fs_template = fs_template
        self.surface = surface
        self.bg_map_type = bg_map_type
        self.fs_home = fs_home

        # Hardcoded vertex counts to avoid loading meshes when using standard templates
        _FS_TEMPLATE_COUNTS = {
            "fsaverage":  163842,
            "fsaverage6": 40962,
            "fsaverage5": 10242,
            "fsaverage4": 2562,
            "fsaverage3": 642,
        }

        self.n_vertices = _FS_TEMPLATE_COUNTS.get(fs_template)
        if self.n_vertices is None:
            raise ValueError(f"Unknown fs_template '{fs_template}'.")

        self.lh = _to_array(lh)
        self.rh = _to_array(rh)

        self._resolve_paths()
        self._load_roi_labels()
        self._build_surface_images()


    def _resolve_paths(self):
        self.use_local = False
        if self.fs_home:
            surf_dir = os.path.join(self.fs_home, "subjects", self.fs_template, "surf")
            self.use_local = os.path.isdir(surf_dir)
            self.surf_dir = surf_dir

        if not self.use_local:
            from nilearn.datasets import fetch_surf_fsaverage
            self._fsavg = fetch_surf_fsaverage(mesh=self.fs_template)
            self.surf_key = "infl" if self.surface == "inflated" else self.surface

        self.mesh_paths, self.bg_paths, self.hemi_present = {}, {}, []

        for hemi_name, data in [("left", self.lh), ("right", self.rh)]:
            if data is None:
                continue

            data = self._match_mesh_size(data, self.n_vertices)
            setattr(self, hemi_name, data)  # write dynamic attribute (self.left / self.right)
            short = hemi_name[0] + "h"

            if self.use_local:
                mesh_path = os.path.join(self.surf_dir, f"{short}.{self.surface}")
                bg_path = (os.path.join(self.surf_dir, f"{short}.{self.bg_map_type}")
                           if self.bg_map_type != "none" else None)
                for p in filter(None, [mesh_path, bg_path]):
                    if not os.path.exists(p):
                        raise FileNotFoundError(f"Surface file not found: {p}")
            else:
                mesh_path = getattr(self._fsavg, f"{self.surf_key}_{hemi_name}")
                bg_path = (str(getattr(self._fsavg, f"{self.bg_map_type}_{hemi_name}"))
                           if self.bg_map_type != "none" else None)

            self.mesh_paths[hemi_name] = mesh_path
            
            if bg_path is not None:
                self.bg_paths[hemi_name] = self._lighten_bg(bg_path, darkness=0.3)
            
            self.hemi_present.append(hemi_name)

        self.hemis = "both" if len(self.hemi_present) > 1 else self.hemi_present[0]

    def _load_roi_labels(self):
        """Retrieve roi labels and their colors. 
           TMP without local FreeSurfer, uses Destrieux instead of DK but colors are not matched
         - TODO: extract annot shipped with verywise (sysdata) as json and use that instead"""

        self.roi = {}
        self.roi_colors = {
            "unknown":                 "#190519",
            "bankssts":                "#196428",
            "caudalanteriorcingulate": "#7D64A0",
            "caudalmiddlefrontal":     "#641900",
            "corpuscallosum":          "#784632",
            "cuneus":                  "#DC1464",
            "entorhinal":              "#DC140A",
            "fusiform":                "#B4DC8C",
            "inferiorparietal":        "#DC3CDC",
            "inferiortemporal":        "#B42878",
            "isthmuscingulate":        "#8C148C",
            "lateraloccipital":        "#141E8C",
            "lateralorbitofrontal":    "#234B32",
            "lingual":                 "#E18C8C",
            "medialorbitofrontal":     "#C8234B",
            "middletemporal":          "#A06432",
            "parahippocampal":         "#14DC3C",
            "paracentral":             "#3CDC3C",
            "parsopercularis":         "#DCB48C",
            "parsorbitalis":           "#146432",
            "parstriangularis":        "#DC3C14",
            "pericalcarine":           "#78643C",
            "postcentral":             "#DC1414",
            "posteriorcingulate":      "#DCB4DC",
            "precentral":              "#3C14DC",
            "precuneus":               "#A08CB4",
            "rostralanteriorcingulate":"#50148C",
            "rostralmiddlefrontal":    "#4B327D",
            "superiorfrontal":         "#14DCA0",
            "superiorparietal":        "#14B48C",
            "superiortemporal":        "#8CDCDC",
            "supramarginal":           "#50A014",
            "frontalpole":             "#640064",
            "temporalpole":            "#4614AA",
            "transversetemporal":      "#9696C8",
            "insula":                  "#FFC020",
        }

        if self.use_local:
            import nibabel as nib
            for hemi_name in self.hemi_present:
                annot = os.path.join(self.surf_dir, "..", "label", f"{hemi_name[0]}h.aparc.annot")
                roi_map, _, names = nib.freesurfer.read_annot(annot)
                names = [n.decode() if isinstance(n, bytes) else n for n in names]
                # nib encodes "unknown" as level -1, everything else 0..len(names)-2
                roi_labels = {
                    name: (i, self.roi_colors.get(name, "#000000"))
                    for i, name in enumerate(names)
                }
                self.roi[hemi_name] = (roi_map, roi_labels)
        else:
            from nilearn.datasets import fetch_atlas_surf_destrieux
            destrieux = fetch_atlas_surf_destrieux()
            for hemi_name in self.hemi_present:
                labels = np.asarray(destrieux[f"map_{hemi_name}"], dtype=int)
                names = [n.decode() if isinstance(n, bytes) else str(n)
                         for n in destrieux["labels"]]
                self.roi[hemi_name] = (labels, names)

    def _build_surface_images(self):
        """SurfaceImage build (parse mesh geometry)"""
        from nilearn.surface.surface import PolyMesh
        from nilearn.surface import SurfaceImage

        data_parts = {h: getattr(self, h) for h in self.hemi_present}
        poly_mesh = PolyMesh(**self.mesh_paths) # parse all geometry files 

        # Data 
        self.surf_img = SurfaceImage(mesh=poly_mesh, data=data_parts)

        # Background
        self.bg_img = (SurfaceImage(mesh=poly_mesh, data=self.bg_paths)
                        if self.bg_paths else None)

        # the ROI atlas (to hand to add_contours()...?)
        roi_parts = {h: self.roi[h][0] for h in self.hemi_present}
        self.roi_img = SurfaceImage(mesh=poly_mesh, data=roi_parts)

    # --- Convenience ---------------------------------------------
    @staticmethod
    def _match_mesh_size(data, target_len):
        """Crop data array to match mesh vertex count (or raise if too short)."""
        if len(data) > target_len:
            return data[:target_len]
        if len(data) < target_len:
            raise ValueError(f"'{target_len}' vertices is higher than the data resolution.")
        return data

    @staticmethod
    def _lighten_bg(bg_path, darkness=0.3):
        """Compress bg_map dynamic range toward mid-grey (0.5)."""
        bg = load_surf_data(bg_path) if isinstance(bg_path, str) else bg_path
        if bg is None:
            return None
        arr = np.asarray(bg, dtype=np.float32)
        lo, hi = np.nanmin(arr), np.nanmax(arr)
        if hi > lo:
            arr = (arr - lo) / (hi - lo) # normalize to [0,1] first
        # blend toward 0.5 (mid-grey): darkness=1 -> full range, darkness=0 -> flat grey
        return 0.5 + (arr - 0.5) * darkness

    def scene_bbox(self):
        """Uniform cubic bounding box (mid, half) across all present hemis."""
        all_coords = np.vstack([load_surf_mesh(p)[0] for p in self.mesh_paths.values()])
        mid = (all_coords.max(axis=0) + all_coords.min(axis=0)) / 2.0
        half = (all_coords.max(axis=0) - all_coords.min(axis=0)).max() / 2.0 * 1.05
        return mid, half
    
    def roi_levels(self, roi_names):
        """Resolve DK region names to (levels, labels, colors). Assumes the same
        DK/aparc colortable order on both hemispheres."""
        _, roi_labels = self.roi[self.hemi_present[0]]
        levels, labels, colors = [], [], []
        for name in roi_names:
            if name in roi_labels:
                level, color = roi_labels[name]
                levels.append(level)
                labels.append(name)
                colors.append(color)
        return levels, labels, colors


# ═════════════════════════════════════════════════════════════════════════════
# Colour scale, ticks and colorbar/density traces
# ═════════════════════════════════════════════════════════════════════════════

class VeryWiseColormap:
    """
    Resolves a colormap + color limits + colorbar ticks + band mask from
    vertex-wise hemisphere data, with automatic continuous/discrete detection. 
    Also builds the Plotly traces for a colorbar strip and a density panels.

    Parameters
    ----------
    lh, rh : np.ndarray or None
        Hemisphere value arrays. NaNs are ignored.
    lh_mask, rh_mask : np.ndarray or None
        Boolean significance masks, same shape as lh/rh. None = all
        significant (no thresholding).
    cmap : str, list, Colormap, or None
    vmin, vmax : float or None
        User-specified colormap edges. When given, both the colormap norm
        and the panel axis range use these values. When None:
          - panel axis range (self.vmin/self.vmax) = min/max of all_vals
          - colormap norm edges (self.eff_min/self.eff_max) = min/max of
            sig_vals (the significant subset)
    max_discrete : int
    n_y, n_x : int
        Grid resolution for colorbar / density heatmaps.

    Attributes
    ----------
    cmap : Colormap
    vmin, vmax : float
        Panel axis extent (y_grid, ticks, colorbar/density range).
    sig_min, sig_max : float
        Colormap normalization edges.
    norm : matplotlib.colors.Normalize or BoundaryNorm
    discrete : bool
    tick_vals, tick_text : list
    y_grid : np.ndarray
    band_mask : np.ndarray of bool
    colorbar_trace, density_fill_trace, density_outline_trace : Plotly traces
    """

    def __init__(self, lh, rh=None, lh_mask=None, rh_mask=None, cmap=None,
                 vmin=None, vmax=None, max_discrete=20, n_y=1000, n_x=80):

        all_vals, sig_vals = self._masked_concat(lh, rh, lh_mask, rh_mask)

        # +/- Inf are already taken care of in _to_array(), set to NA
        self.data_min, self.data_max = float(np.nanmin(all_vals)), float(np.nanmax(all_vals))

        self.vmin = vmin if vmin is not None else self.data_min
        self.vmax = vmax if vmax is not None else self.data_max

        self.is_empty = np.all(np.isnan(sig_vals))

        if self.is_empty:
            self.cmap = self._wrap_nan_transparent(matplotlib.colormaps["binary"])
            self.eff_min, self.eff_max = self.vmin, self.vmax
            self.norm = mcolors.Normalize(vmin=self.vmin, vmax=self.vmax)
            self.is_discrete = False
            self.tick_vals, self.tick_text = self._colorbar_ticks(self.vmin, self.vmax)
            self.y_grid = np.linspace(self.vmin, self.vmax, n_y)
            self.band_mask = np.ones_like(self.y_grid, dtype=bool)
            self.colorscale = self._to_plotly_colorscale()
            self.colorbar_trace = self._build_colorbar_trace()
            self.density_fill_trace, self.density_outline_trace = self._build_density_traces(all_vals[np.isfinite(all_vals)], n_x)
            return

        sig_finite = sig_vals[np.isfinite(sig_vals)] # drop non-significant 
        self.sig_min, self.sig_max = float(sig_finite.min()), float(sig_finite.max())

        self.is_discrete = self._detect_discrete(sig_finite, max_discrete)

        self.cmap, self.eff_min, self.eff_max, self.norm = self._build_cmap(
            cmap, self.is_discrete, vmin, vmax, self.sig_min, self.sig_max)

        if self.is_discrete:
            self.tick_vals = list(range(1, int(self.sig_max) + 1))
            self.tick_text = [str(v) for v in self.tick_vals]
        else:
            self.tick_vals, self.tick_text = self._colorbar_ticks(self.vmin, self.vmax)

        # Now draw the colobar and the density
        self.y_grid = np.linspace(self.vmin, self.vmax, n_y)
        self.colorscale = self._to_plotly_colorscale()

        self.band_mask = self._resolve_band_mask(sig_finite, self.y_grid)

        self.colorbar_trace = self._build_colorbar_trace()
        self.density_fill_trace, self.density_outline_trace = \
            self._build_density_traces(all_vals[np.isfinite(all_vals)], n_x)
    
    
    @staticmethod
    def _masked_concat(lh, rh, lh_mask, rh_mask):
        pairs = [(a, m if m is not None else np.ones_like(a, dtype=bool))
                 for a, m in ((lh, lh_mask), (rh, rh_mask)) if a is not None]
        all_arrs = [a for a, m in pairs]
        masked_arrs = [np.where(np.asarray(m, dtype=bool), a, np.nan) for a, m in pairs]
        all_vals = np.concatenate(all_arrs) if all_arrs else np.array([])
        masked_vals = np.concatenate(masked_arrs) if masked_arrs else np.array([])
        return all_vals, masked_vals

    @staticmethod
    def _detect_discrete(vals, max_discrete=20):
        if vals.size == 0:
            return False
        unique = np.unique(vals)
        if len(unique) > max_discrete:
            return False
        return bool(np.allclose(vals, np.round(vals)))

    @classmethod
    def _build_cmap(self, cmap, is_discrete, vmin, vmax, sig_min, sig_max):

        if is_discrete:
            n_unique = int(sig_max)
            # Pick a color
            palette = cmap if cmap is not None else ("tab10" if n_unique <= 10 else "tab20")

            resolved_min = 0.5 if vmin is None else vmin
            resolved_max = n_unique + 0.5 if vmax is None else vmax

            boundaries = np.arange(1.5, resolved_max)
            norm = mcolors.BoundaryNorm(boundaries, ncolors=n_unique)

            return self._discrete_cmap(n_unique, palette), resolved_min, resolved_max, norm

        if cmap is not None:
            base = self._resolve_cmap_obj(cmap)
        else: 
            if sig_min >= 0:
                base = matplotlib.colormaps["hot_r"]
            elif sig_max <= 0:
                base = matplotlib.colormaps["viridis"]
            else:
                base = self._diverging_cmap()
        
        resolved_min = vmin if vmin is not None else sig_min
        resolved_max = vmax if vmax is not None else sig_max

        norm = mcolors.Normalize(vmin=resolved_min, vmax=resolved_max)

        return self._wrap_nan_transparent(base), resolved_min, resolved_max, norm
    
    # --- Convenience ---------------------------------------------
    @staticmethod
    def _resolve_cmap_obj(cmap):
        if isinstance(cmap, list):
            return mcolors.ListedColormap(cmap)
        if isinstance(cmap, mcolors.Colormap):
            return cmap
        return matplotlib.colormaps[cmap]
    
    @staticmethod
    def _wrap_nan_transparent(cmap_obj):
        c = cmap_obj.copy()
        c.set_bad(alpha=0.0)
        return c

    @classmethod
    def _discrete_cmap(self, n_unique, palette):
        palette_obj = self._resolve_cmap_obj(palette)
        if hasattr(palette_obj, "colors") and n_unique <= len(palette_obj.colors):
            colors = [palette_obj.colors[i] for i in range(n_unique)]
        else:
            colors = [palette_obj(i / max(n_unique - 1, 1)) for i in range(n_unique)]
        return self._wrap_nan_transparent(mcolors.ListedColormap(colors, name="discrete"))

    @staticmethod
    def _diverging_cmap():
        hot_r = matplotlib.colormaps["hot_r"]
        viridis = matplotlib.colormaps["viridis"]
        colors = [viridis(0.0), viridis(0.5), viridis(1.0),
                  hot_r(0.5), hot_r(0.75), hot_r(1.0)]
        nodes = [0.0, 0.25, 0.5, 0.6, 0.75, 1.0]
        return mcolors.LinearSegmentedColormap.from_list(
            "hot_r_viridis", list(zip(nodes, colors)))
    
    # ---- Ticks --------------------------------------------------

    @classmethod
    def _colorbar_ticks(self, vmin, vmax, n_ticks=10):
        data_range = vmax - vmin
        if data_range == 0:
            return [vmin], [self._fmt_tick(vmin)]

        is_integer = (float(vmin).is_integer() and float(vmax).is_integer()
                      and data_range <= 2 * n_ticks)
        locator = MaxNLocator(nbins=n_ticks - 1, integer=is_integer,
                               steps=[1, 2, 2.5, 5, 10])
        vals = locator.tick_values(vmin, vmax)
        vals = vals[(vals >= vmin - data_range * 1e-9) &
                    (vals <= vmax + data_range * 1e-9)]

        tol = data_range * 1e-9
        vals = np.where(np.abs(vals) < tol, 0.0, vals)

        must_include = [vmin, vmax] + ([0.0] if vmin < 0 < vmax else [])
        for v in must_include:
            if not np.any(np.abs(vals - v) < tol):
                vals = np.sort(np.append(vals, v))

        labels = [self._fmt_tick(v, vmin, vmax) for v in vals]
        return vals.tolist(), labels

    @staticmethod
    def _fmt_tick(v, vmin=None, vmax=None):
        if v == 0.0:
            return "0"
        data_range = abs(vmax - vmin) if (vmin is not None and vmax is not None) else abs(v)
        magnitude = abs(v)
        if magnitude != 0 and (magnitude >= 1e4 or magnitude < 1e-3):
            return f"{v:.2e}"
        if float(v).is_integer():
            return str(int(v))
        if data_range != 0:
            decimals = min(max(0, -int(np.floor(np.log10(data_range))) + 2), 4)
            return f"{v:.{decimals}f}".rstrip("0").rstrip(".")
        return f"{v:.3g}"

    # --- Plotting ------------------------------------------------

    def _to_plotly_colorscale(self, n=256):
        return [
            [i / (n - 1), f"rgba({int(r*255)},{int(g*255)},{int(b*255)},{a:.3f})"]
            for i, (r, g, b, a) in enumerate(self.cmap(np.linspace(0, 1, n)))
        ]

    def _resolve_band_mask(self, sig_vals, y_grid):
        """Blank the y_grid region where there's no significant values."""
        if sig_vals.size == 0:
            return np.ones_like(y_grid, dtype=bool)

        pos_sig = sig_vals[sig_vals > 0]
        neg_sig = sig_vals[sig_vals < 0]

        lo = neg_sig.max() if neg_sig.size > 0 else self.vmin
        hi = pos_sig.min() if pos_sig.size > 0 else self.vmax

        gap = (y_grid >= lo) & (y_grid <= hi)
        outside_raw = (y_grid < self.sig_min) | (y_grid > self.sig_max)
        return gap | outside_raw

    def _build_colorbar_trace(self):
        z_cb = np.tile(self.y_grid.reshape(-1, 1), (1, 2)).astype(float)
        if self.band_mask.any():
            z_cb[self.band_mask, :] = np.nan
        return go.Heatmap(
            x=[0, 1], y=self.y_grid, z=z_cb,
            zmin=self.eff_min, zmax=self.eff_max, colorscale=self.colorscale,
            showscale=False, hoverinfo="skip", xgap=0, ygap=0,
        )

    def _build_density_traces(self, vals, n_x):
        from scipy.stats import gaussian_kde

        kde = None
        if not self.is_discrete and vals.size > 1 and np.ptp(vals) > 0:
            try:
                kde = gaussian_kde(vals)
                density = kde(self.y_grid)
            except np.linalg.LinAlgError:
                kde = None
                density = np.zeros_like(self.y_grid)
        else:
            density = np.zeros_like(self.y_grid)

        dmax = density.max() if density.max() > 0 else 1.0
        density_norm = density / dmax

        outside = (self.y_grid < self.data_min) | (self.y_grid > self.data_max)

        density_fill = density_norm.copy()
        density_fill[outside] = 0.0

        density_line = density_norm.copy()
        density_line[outside] = np.nan

        y_line = self.y_grid.copy()
        if kde is not None:
            edge_vals = kde([self.data_min, self.data_max]) / dmax
            y_line = np.concatenate([[self.data_min], y_line, [self.data_max]])
            density_line = np.concatenate([[edge_vals[0]], density_line, [edge_vals[1]]])
            order = np.argsort(y_line)
            y_line, density_line = y_line[order], density_line[order]

        self.density_norm = density_fill

        x_grid = np.linspace(0, 1, n_x)
        z_kde = np.tile(self.y_grid.reshape(-1, 1), (1, n_x)).astype(float)
        z_kde[x_grid.reshape(1, -1) > density_fill.reshape(-1, 1)] = np.nan
        if self.band_mask.any():
            z_kde[self.band_mask, :] = np.nan

        outline_trace = go.Scatter(
            x=density_line, y=self.y_grid, mode="lines", line=dict(color="black", width=1.5),
            hoverinfo="skip", showlegend=False,
        )

        fill_trace = go.Heatmap(
            x=x_grid, y=self.y_grid, z=z_kde,
            zmin=self.eff_min, zmax=self.eff_max, colorscale=self.colorscale,
            showscale=False, hoverinfo="skip", xgap=0, ygap=0,
        )

        return fill_trace, outline_trace

    def zero_line_trace(self):
        if not (self.vmin < 0 < self.vmax):
            return None
        return go.Scatter(
            x=[0, 1], y=[0, 0], mode="lines",
            line=dict(color="grey", width=0.75, dash="dash"),
            hoverinfo="skip", showlegend=False,
        )

    def __repr__(self):
        if self.is_empty:
            return "VeryWiseColormap(empty=True, no finite data)"
        mode = "discrete" if self.is_discrete else "continuous"
        return (f"VeryWiseColormap(mode={mode}, axis=[{self.vmin:.3g}, {self.vmax:.3g}], "
                f"cmap=[{self.eff_min:.3g}, {self.eff_max:.3g}], n_ticks={len(self.tick_vals)})")


# ═════════════════════════════════════════════════════════════════════════════
# Panel-grid resolution and brain views 
# ═════════════════════════════════════════════════════════════════════════════

class VeryWiseLayout:
    """
    Resolves the subplot panel grid from a `hemis` + `views` spec, enforces 
    defaults for both modes (static or interactive)

    Attributes
    ----------
    n_rows, n_cols : int
    panels : list of (row, col, hemi, view, title_or_None)
    """
    DEFAULT_LAYOUT = {  # noqa: RUF012
        "static": {
            "views": ["lateral", "dorsal", "anterior",   # Left lateral | Right lateral | Dorsal  | Anterior
                      "medial", "ventral", "posterior"], # Left medial  | Right medial  | Ventral | Posterior
            "panel_size": (400, 440), # brain panel width, height in px
            "colorbar_width": 0.40,
            "colorbar_layout": [("spacer", 0.39), ("colorbar", 0.17), ("spacer", 0.03), ("density", 0.41)]
        },
        "interactive": {
            "views": ["lateral"],
            "panel_size": (500, 350),
            "colorbar_width": 0.20,
            "colorbar_layout": [("colorbar", 0.40), ("density", 0.60)]
        }
    }

    def __init__(self, hemis, views="all", n_cols=None,  mode="static"):
        self.mode = mode
        self.config = self.DEFAULT_LAYOUT[mode]

        if views == "all" or (isinstance(views, list) and len(views) == 6):
            views = self.config["views"] 

        if mode == "static":
            self._build_static_panels(hemis, views, n_cols)
        else:
            self._build_interactive_panels(hemis, views)

    def _build_static_panels(self, hemis, views, n_cols):
        bilateral = hemis == "both"

        # Expand each view name to one or two (hemi, view) pairs
        expanded = []
        for v in views:
            if bilateral and v not in ["dorsal", "anterior", "ventral", "posterior"]:
                expanded.extend([("left", v), ("right", v)])
            else:
                expanded.append((hemis, v))

        n_panels = len(expanded)

        # Auto-detect n_cols: 
        # 4 for bilateral, 3 for single-hemi (capped at actual panel count)
        self.n_cols = min(4 if bilateral else 3, n_panels) if n_cols is None else n_cols
        self.n_rows = math.ceil(n_panels / self.n_cols)

        # Title logic:
        #   left/right: show "Left"/"Right" as column header only on first appearance
        #   both/single: always show the view name (each panel is a distinct view)
        col_hemi_seen = set()
        self.panels = []
        for i, (hemi, view) in enumerate(expanded):
            row, col = divmod(i, self.n_cols)
            if hemi in ("left", "right"):
                key = f"{col}_{hemi}"
                panel_title = hemi.capitalize() if key not in col_hemi_seen else None
                col_hemi_seen.add(key)
            else:
                panel_title = view.capitalize()
            self.panels.append((row, col, hemi, view, panel_title))


    def _build_interactive_panels(self, hemis, views):
        hemi_order = ["left", "right"] if hemis == "both" else [hemis]
        self.n_rows, self.n_cols = len(hemi_order), len(views)
        self.panels = [
            (r, c, hemi, view, f"{hemi.capitalize()} hemisphere")
            for r, hemi in enumerate(hemi_order)
            for c, view in enumerate(views)
        ]

    def hemis_used(self):
        return {h for _, _, h, _, _ in self.panels}
    
    def build_colormap_panel(self, colorbar, colorbar_width):
        """Returns (col_widths, cbar_specs, cbar_col, dens_col)."""
        if not colorbar:
            return [], [], [], None, None
        
        colorbar_width = colorbar_width or self.config["colorbar_width"]
        tags = self.config["colorbar_layout"]

        col_widths = [colorbar_width * frac for _, frac in tags]
        cbar_col = self.n_cols + next(i for i, (t, _) in enumerate(tags) if t == "colorbar") + 1
        dens_col = self.n_cols + next(i for i, (t, _) in enumerate(tags) if t == "density") + 1

        # The first row spans all rows for the colorbar columns
        cbar_specs = [[{"type": "xy", "rowspan": self.n_rows}] * len(tags)]
        cbar_specs += [[None] * len(tags) for _ in range(1, self.n_rows)]

        return col_widths, cbar_specs, cbar_col, dens_col

    def build_figure(self, colorbar, colorbar_width, horizontal_spacing=0.000, vertical_spacing=0.04):
        """Builds the make_subplots figure with scene + colorbar/density columns."""
        cb_widths, cbar_specs, cbar_col, dens_col = self.build_colormap_panel(
            colorbar, colorbar_width)
        col_widths = [1.0] * self.n_cols + cb_widths
        total_cols = len(col_widths)

        specs = [
            [{"type": "scene"}] * self.n_cols + cbar_specs[r] if cbar_specs else [{"type": "scene"}] * self.n_cols
            for r in range(self.n_rows)
        ]

        subplot_titles = [None] * (self.n_rows * total_cols)
        for row, col, _, _, lbl in self.panels:
            if lbl:
                subplot_titles[row * total_cols + col] = lbl

        fig = make_subplots(
            rows=self.n_rows, cols=total_cols, specs=specs, column_widths=col_widths,
            subplot_titles=subplot_titles,
            horizontal_spacing=horizontal_spacing, vertical_spacing=vertical_spacing,
        )
        return fig, cbar_col, dens_col

    def figure_size(self, colorbar=True, colorbar_width=None, cell_px=None,
                     roi_names=None, legend_extra_px=None):
        """
        Resolve (width, height, margin_r) for this layout + colorbar + legend.

        cell_px : (width, height) tuple or single number, optional
            Overrides the mode's default panel size. A single number scales
            both width and height uniformly (height = cell_px * 1.1 for static,
            cell_px for interactive... actually just pass a tuple for full control).
        colorbar_width : float, optional
            Fraction of one panel's width reserved for colorbar/density.
            Defaults to 0.40 (static) / 0.20 (interactive) if not given.
        """
        default_w, default_h = self.config["panel_size"]

        if cell_px is None:
            w_unit, h_unit = default_w, default_h
        elif isinstance(cell_px, (tuple, list)):
            w_unit, h_unit = cell_px
        else:
            w_unit = cell_px
            h_unit = cell_px * (default_h / default_w)  # keep aspect ratio

        width = w_unit * self.n_cols
        height = h_unit * self.n_rows

        if colorbar:
            width += w_unit * (colorbar_width or self.config["colorbar_width"])

        if legend_extra_px is None:
            legend_extra_px = (80 + max(len(n) for n in roi_names) * 7) if roi_names else 0

        width += legend_extra_px
        margin_r = max(legend_extra_px, 10)

        return width, height, margin_r

    @staticmethod
    def scene_key(row, col, n_cols):
        idx = row * n_cols + col
        return "scene" if idx == 0 else f"scene{idx + 1}"


class VeryWiseCameraPresets:
    """
    Camera presets: manual settings for view angles and zoom.

    Tuning guide
    ------------
    'eye'      -- direction vector toward the camera. Rotates the brain.
    'distance' -- Smaller = more zoomed in (eye vector scaled to this length).
    'up'       -- For lateral/medial/anterior/posterior: (0,0,1) = superior up.
                  For dorsal/ventral: (-1,0,0) / (1,0,0).
    'center'   -- Shifts the brain -- always (0,0,0) except for inflated surfaces.

    This is a class-level constant registry (not instance-configurable);
    `resolve()` supports a per-call `custom_presets` override on top of it.
    """

    @staticmethod
    def _bilateral(d):
        """Expand ("both", view) entries to also cover ("left"/"right", view)."""
        extra = {}
        for (hemi, view), v in d.items():
            if hemi == "both":
                extra.setdefault(("left", view), v)
                extra.setdefault(("right", view), v)
        return {**extra, **d}

    PRESETS = { 
        "pial": _bilateral.__func__({
            ("left", "lateral"):  dict(eye=(-1.0, 0.0, 0.0), up=(0, 0, 1), distance=1.20),
            ("right", "lateral"): dict(eye=( 1.0, 0.0, 0.0), up=(0, 0, 1), distance=1.20),
            ("left", "medial"):   dict(eye=( 1.0, 0.0, 0.0), up=(0, 0, 1), distance=1.05),
            ("right", "medial"):  dict(eye=(-1.0, 0.0, 0.0), up=(0, 0, 1), distance=1.05),
            ("both", "dorsal"):   dict(eye=( 0.0, 0.0, 1.0), up=(-1, 0, 0), distance=1.15),
            ("both", "ventral"):  dict(eye=( 0.0, 0.0, -1.0), up=(1, 0, 0), distance=1.30),
            ("both", "anterior"): dict(eye=( 0.0, 1.0, 0.0), up=(0, 0, 1), distance=1.15),
            ("both", "posterior"):dict(eye=( 0.0,-1.0, 0.0), up=(0, 0, 1), distance=1.15),
        }),
        "inflated": _bilateral.__func__({
            ("left", "lateral"):  dict(eye=(-1.0, 0.0, 0.0), up=(0, 0, 1), distance=1.10),
            ("right", "lateral"): dict(eye=( 1.0, 0.0, 0.0), up=(0, 0, 1), distance=1.10),
            ("left", "medial"):   dict(eye=( 1.0, 0.0, 0.0), up=(0, 0, 1), distance=1.20),
            ("right", "medial"):  dict(eye=(-1.0, 0.0, 0.0), up=(0, 0, 1), distance=1.20),
            ("both", "dorsal"):   dict(eye=(0.11, 0.0, 1.25), up=(-1, 0, 0), center=(0.15, 0, 0), distance=1.15),
            ("both", "ventral"):  dict(eye=(0.25, 0.0, -1.3), up=(1, 0, 0),  center=(0.18, 0, 0), distance=1.30),
            ("both", "anterior"): dict(eye=(0.16, 1.1, -0.1), up=(0, 0.2, 1),center=(0.20, 0, 0), distance=1.15),
            ("both", "posterior"):dict(eye=(0.20,-1.1, -0.0), up=(0, 0, 1),  center=(0.18, 0, 0), distance=1.15),
        }),
        "white": _bilateral.__func__({  # white sits between pial and inflated
            ("left", "lateral"):  dict(eye=(-1.0, 0.0, 0.0), up=(0, 0, 1), distance=1.20),
            ("right", "lateral"): dict(eye=( 1.0, 0.0, 0.0), up=(0, 0, 1), distance=1.20),
            ("left", "medial"):   dict(eye=( 1.0, 0.0, 0.0), up=(0, 0, 1), distance=1.05),
            ("right", "medial"):  dict(eye=(-1.0, 0.0, 0.0), up=(0, 0, 1), distance=1.05),
            ("both", "dorsal"):   dict(eye=( 0.0, 0.0, 1.0), up=(-1, 0, 0), distance=1.18),
            ("both", "ventral"):  dict(eye=( 0.0, 0.0, -1.0), up=(1, 0, 0), distance=1.32),
            ("both", "anterior"): dict(eye=( 0.0, 1.0, 0.0), up=(0, 0, 1), distance=1.15),
            ("both", "posterior"):dict(eye=( 0.0,-1.0, 0.0), up=(0, 0, 1), distance=1.15),
        }),
    }

    @classmethod
    def resolve(self, hemi, view, surface, custom_presets=None):
        """
        Resolve camera dict for a (hemi, view, surface) triple.
        Lookup order:
          1. custom_presets[surface][(hemi, view)]  -- per-call override
          2. self.PRESETS[surface][(hemi, view)]    -- surface-specific constant
          3. self.PRESETS["inflated"][(hemi, view)] -- last-resort fallback
        center is always locked to origin regardless of source.
        """
        surface_key = surface if surface in self.PRESETS else "inflated"

        p = None
        for source in [custom_presets, self.PRESETS]:
            if source and surface_key in source:
                p = source[surface_key].get((hemi, view))
            if p is not None:
                break

        if p is None:
            raise ValueError(
                f"No camera preset for surface={surface!r}, hemi={hemi!r}, view={view!r}. "
                f"Add it to VeryWiseCameraPresets.PRESETS[{surface_key!r}]."
            )

        ex, ey, ez = p["eye"]
        ux, uy, uz = p["up"]
        d = p["distance"]
        mag = math.sqrt(ex**2 + ey**2 + ez**2) or 1.0
        s = d / mag
        cx, cy, cz = p.get("center", (0.0, 0.0, 0.0))

        return {"eye": dict(x=ex * s, y=ey * s, z=ez * s),
                "up": dict(x=ux, y=uy, z=uz),
                "center": dict(x=cx, y=cy, z=cz)}


# ═════════════════════════════════════════════════════════════════════════════
# VeryWiseSurfacePlotter : shared scaffolding for static and interactive render
# ═════════════════════════════════════════════════════════════════════════════

class VeryWiseSurfacePlotter:
    """
    Shared scaffolding used by both StaticSurfacePlotter and
    InteractiveSurfacePlotter: scene bbox / axis-range helpers and figure
    unwrap. Trace-building is intentionally NOT shared -- static reuses one
    cached trace set per hemi across multiple camera angles, while
    interactive needs fresh per-vertex hover strings per (hemi, view).
    """

    def __init__(self, mesh: "VeryWiseMesh", colormap: "VeryWiseColormap",
                 surface, camera_presets=None, roi_names=None):
        self.mesh = mesh
        self.colormap = colormap
        self.surface = surface
        self.camera_presets = camera_presets
        self.roi_names = roi_names
        self.mid, self.half = mesh.scene_bbox()

    def axis(self, i):
        return dict(visible=False, showgrid=False, zeroline=False,
                    range=[float(self.mid[i] - self.half), float(self.mid[i] + self.half)])

    def camera(self, hemi, view):
        return VeryWiseCameraPresets.resolve(
            hemi, view, surface=self.surface, custom_presets=self.camera_presets)

    @staticmethod
    def _to_plotly(fp):
        """Unwrap nilearn's PlotlySurfaceFigure -> raw plotly go.Figure."""
        return fp.figure if hasattr(fp, "figure") else fp

    def _add_colorbar_and_density(self, fig, cb_col, dens_col, colorbar_label):
        sc = self.colormap
        fig.add_trace(sc.colorbar_trace, row=1, col=cb_col)
        fig.add_trace(sc.density_fill_trace, row=1, col=dens_col)
        fig.add_trace(sc.density_outline_trace, row=1, col=dens_col)

        zero_trace = sc.zero_line_trace()
        if zero_trace is not None:
            fig.add_trace(zero_trace, row=1, col=dens_col)

        fig.update_xaxes(range=[0, 1], visible=True, showgrid=False, showticklabels=False,
                          ticks="", showline=True, linecolor="black", mirror=True,
                          fixedrange=True, row=1, col=cb_col)

        fig.update_yaxes(
            range=[sc.vmin, sc.vmax], side="left", showgrid=False, zeroline=False,
            tickmode="array", tickvals=sc.tick_vals, ticktext=sc.tick_text,
            title=dict(text=colorbar_label or "", standoff=2), ticks="outside",
            showline=True, linecolor="black", mirror=True, automargin=True, row=1, col=cb_col)

        fig.update_xaxes(range=[0, 1.05], visible=False, showgrid=False,
                          fixedrange=True, row=1, col=dens_col)
        fig.update_yaxes(range=[sc.vmin, sc.vmax], visible=False, showgrid=False,
                          row=1, col=dens_col)
    
    def _add_roi_contours(self, fp, hemi):
        """Overlay ROI boundary contours for the requested region names."""
        if not self.roi_names:
            return
        
        levels, labels, colors = self.mesh.roi_levels(self.roi_names)
        if not levels:
            return
        
        style_dict = [{"width": 5, "color": c} for c in colors]
        fp.add_contours(roi_map=self.mesh.roi_img, levels=levels, labels=labels,
                        lines=style_dict)

    def _add_roi_legend(self, fig):
        """Only show each ROI name once in the legend (rest get showlegend=False)."""
        n_roi = len(self.roi_names) if self.roi_names else 0
        count = 0
        for tr in fig.data:
            if getattr(tr, "name", None) is None:
                continue
            count += 1
            if count > n_roi:
                tr.update(showlegend=False)

# ═════════════════════════════════════════════════════════════════════════════
# Kaleido PNG renderer
# ═════════════════════════════════════════════════════════════════════════════

class VeryWiseStaticPlotter(VeryWiseSurfacePlotter):
    """
    Static multi-view PNG renderer. Builds one `plot_surf()` trace cache per
    hemisphere (reused across camera angles), optionally overlays ROI
    boundary contours (DK atlas) via nilearn's `add_contours`, then composes
    a single multi-scene Plotly figure exported through Kaleido.
    """

    CANONICAL_VIEW = {"left": "lateral", "right": "lateral", "both": "dorsal"}

    def __init__(self, mesh, colormap, surface, layout: "VeryWiseLayout",
                 camera_presets=None, roi_names=None):
        super().__init__(mesh, colormap, surface, camera_presets, roi_names)
        self.layout = layout
        self.trace_cache = self._build_trace_cache()

    def _build_trace_cache(self):
        sc, mesh = self.colormap, self.mesh
        cache = {}
        for hemi in self.layout.hemis_used():
            fp = plot_surf(
                surf_mesh=None, surf_map=mesh.surf_img, hemi=hemi,
                view=self.CANONICAL_VIEW[hemi],
                cmap=sc.cmap, vmin=sc.eff_min, vmax=sc.eff_max,
                threshold=None, symmetric_cmap=None, avg_method=None,
                bg_map=mesh.bg_img, bg_on_data=False, colorbar=False, engine="plotly",
            )

            self._add_roi_contours(fp, hemi)

            raw = self._to_plotly(fp)

            for tr in raw.data:
                tr.update(hoverinfo="skip")
                if "showscale" in tr:
                    tr.update(showscale=False)
            cache[hemi] = list(raw.data)
        return cache

    def render(self, colorbar=True, colorbar_label=None, title=None,
               output_file=None, dpi=150, cell_px=None, colorbar_width=None):

        layout = self.layout
        fig, cbar_col, dens_col = layout.build_figure(colorbar, colorbar_width)
        fig_width, fig_height, margin_r = self.layout.figure_size(
        colorbar=colorbar, colorbar_width=colorbar_width, cell_px=cell_px,
        roi_names=None)  # static legend is bottom-anchored, not right-side

        # Move ROI legend to the bottom 
        for ann in fig.layout.annotations:
            ann.update(y=ann.y - 0.02, yanchor="bottom")

        for row, col, hemi, view, _ in layout.panels:
            for tr in self.trace_cache[hemi]:
                fig.add_trace(tr, row=row + 1, col=col + 1)
            sk = layout.scene_key(row, col, layout.n_cols)
            fig.update_layout(**{sk: dict(
                aspectmode="cube", xaxis=self.axis(0), yaxis=self.axis(1), zaxis=self.axis(2),
                bgcolor="white", camera=self.camera(hemi, view))})

        if colorbar:
            self._add_colorbar_and_density(fig, cbar_col, dens_col, colorbar_label)

        # Only add ROI labels to legend once 
        self._add_roi_legend(fig)

        fig.update_layout(
            title=dict(text=title or "", font=dict(size=25, weight="bold")),
            paper_bgcolor="white", plot_bgcolor="white",
            height=fig_height,
            width=fig_width,
            margin=dict(l=10, r=margin_r, t=max(60 if title else 30, 45), b=10),
            legend=dict(orientation="h", yanchor="top", y=-0.05, xanchor="left", x=0.02, 
                        title=None, font=dict(size=18))
        )

        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        fig.write_image(output_file, scale=dpi / 72.0)
        return output_file


# ═════════════════════════════════════════════════════════════════════════════
# HTML renderer with per-vertex hover
# ═════════════════════════════════════════════════════════════════════════════

class VeryWiseInteractivePlotter(VeryWiseSurfacePlotter):
    """
    Interactive multi-view HTML renderer. Builds its own trace per (hemi,
    view) -- unlike the static renderer, hover text differs per vertex and
    per view, so nothing is cached across views.
    """

    def __init__(self, mesh, colormap, surface, layout: "VeryWiseLayout", 
                 camera_presets=None, roi_names=None):
        super().__init__(mesh, colormap, surface, camera_presets, roi_names)
        self.layout = layout
    
    @staticmethod
    def _hover_text(i, v, roi_map, level_to_name, x, y, z):
        """Per-vertex hover string for interactive plots."""
        if not np.isfinite(v):
            return ""
        line = f"Vertex index: {i + 1}<br>Value: <b>{v:.4f}</b><br>"
        if roi_map is not None:
            lbl_idx = roi_map[i]
            region = level_to_name.get(lbl_idx, "unknown")
            line += f"Region: <b>{region}</b><br>"
        # line += f"Position [x, y, z]: [{x:.2f}, {y:.2f}, {z:.2f}]"
        return line
    
    def render(self, colorbar=True, colorbar_label=None, title=None,
               output_html=None, colorbar_width=None, cell_px=None):
        mesh, sc, layout = self.mesh, self.colormap, self.layout
        fig, cbar_col, dens_col = layout.build_figure(
            colorbar, colorbar_width, horizontal_spacing=0.02)

        fig_width, fig_height, margin_r = self.layout.figure_size(
        colorbar=colorbar, colorbar_width=colorbar_width, cell_px=cell_px,
        roi_names=self.roi_names)

        for row, col, hemi, view, _ in layout.panels:
            data = getattr(mesh, hemi)
            mesh_path = mesh.mesh_paths[hemi]
            bg = mesh.bg_paths.get(hemi)
            # Pre-build the reverse lookup once per hemisphere, instead of reversing the dict
            # on every vertex call for hover 
            roi_map, roi_labels = mesh.roi.get(hemi, (None, {}))
            level_to_name = {level: name for name, (level, _color) in roi_labels.items()}

            fp = plot_surf(
                surf_mesh=mesh_path, surf_map=data, hemi=hemi, view=view,
                cmap=sc.cmap, vmin=sc.eff_min, vmax=sc.eff_max, threshold=None,
                symmetric_cmap=None, avg_method=None,
                bg_map=bg, bg_on_data=False, colorbar=False, engine="plotly")
            self._add_roi_contours(fp, hemi)
            raw = self._to_plotly(fp)

            # ----- Hover behavior ---------------------------------------------
            xs, ys, zs = (np.asarray(raw.data[0].x),
                          np.asarray(raw.data[0].y),
                          np.asarray(raw.data[0].z))

            text_vals = [self._hover_text(i, v, roi_map, level_to_name, x, y, z)
                for i, (v, x, y, z) in enumerate(zip(data, xs, ys, zs))]

            for trace in raw.data:
                if trace.type == "mesh3d":
                    trace.update(text=text_vals, hovertemplate="%{text}<extra></extra>",
                                hoverlabel=dict(bgcolor="silver", font=dict(color="black")),
                                showlegend=False)
                else:
                    trace.update(hoverinfo="skip") # ROI contour lines: skip vertex hover
                fig.add_trace(trace, row=row + 1, col=col + 1)

            sk = layout.scene_key(row, col, layout.n_cols)
            if hasattr(raw.layout, "scene") and raw.layout.scene.camera:
                fig.update_layout(**{sk: {"camera": raw.layout.scene.camera.to_plotly_json()}})

        self._add_roi_legend(fig)

        if colorbar:
            self._add_colorbar_and_density(fig, cbar_col, dens_col, colorbar_label)

        ax_style = dict(visible=False, showgrid=False, zeroline=False, showticklabels=False, showspikes=False)
        fig.update_scenes(xaxis=ax_style, yaxis=ax_style, zaxis=ax_style, 
            bgcolor="rgba(0,0,0,0)", aspectmode="data")

        fig.update_layout(
            title_text=title or "", title_font_size=14,
            paper_bgcolor="white", plot_bgcolor="white", hovermode="closest",
            height=fig_height,
            width=fig_width,
            legend=dict(orientation="v", yanchor="top", y=1, xanchor="left", x=1.02,
                        title=None, font=dict(size=12), tracegroupgap=4),
            margin=dict(l=10, r=margin_r, t=45, b=10),
        )

        os.makedirs(os.path.dirname(os.path.abspath(output_html)), exist_ok=True)
        fig.write_html(output_html, include_plotlyjs="cdn", full_html=True, auto_open=False)
        return output_html

# ═════════════════════════════════════════════════════════════════════════════
# ── Public API : thin wrapper functions ─────────────────────────────────────
# ═════════════════════════════════════════════════════════════════════════════


def vw_surf_static_plotly(lh, rh, lh_mask=None, rh_mask=None, vmin=None, vmax=None,
    views='all', surface='pial', bg_map_type='sulc', fs_template=None, fs_home=None,
    colorbar=True, colorbar_label=None, colorbar_width=0.40, cmap=None, 
    roi_names=None, title=None, output_file=None, dpi=150, cell_px=None):
    """
    Render a static multi-view PNG of vertex-wise surface data.
    """
    lh_arr, rh_arr = _to_array(lh), _to_array(rh)

    colormap = VeryWiseColormap(lh_arr, rh_arr, lh_mask=lh_mask, rh_mask=rh_mask,
                                 cmap=cmap, vmin=vmin, vmax=vmax)

    lh_plot = _threshold(lh_arr, lh_mask) 
    rh_plot = _threshold(rh_arr, rh_mask)

    mesh = VeryWiseMesh(lh_plot, rh_plot, fs_template, surface, bg_map_type, fs_home)
    layout = VeryWiseLayout(mesh.hemis, views, mode="static")

    plotter = VeryWiseStaticPlotter(mesh, colormap, surface, layout,
                                    roi_names=roi_names)

    return plotter.render(colorbar=colorbar, colorbar_label=colorbar_label, title=title,
                          output_file=output_file, dpi=dpi, cell_px=cell_px, colorbar_width=colorbar_width)


def vw_surf_interactive(lh, rh, lh_mask=None, rh_mask=None, vmin=None, vmax=None,
    views='all', surface='pial', bg_map_type='sulc', fs_template=None, fs_home=None,
    colorbar=True, colorbar_label=None, colorbar_width=0.20, cmap=None, 
    roi_names=None, title=None, output_html=None, cell_px=None):
    
    """Render an interactive multi-view HTML surface plot with per-vertex hover."""
    lh_arr, rh_arr = _to_array(lh), _to_array(rh)

    colormap = VeryWiseColormap(lh_arr, rh_arr, lh_mask=lh_mask, rh_mask=rh_mask,
                                cmap=cmap, vmin=vmin, vmax=vmax)

    lh_plot = _threshold(lh_arr, lh_mask) 
    rh_plot = _threshold(rh_arr, rh_mask) 
    
    mesh = VeryWiseMesh(lh_plot, rh_plot, fs_template, surface, bg_map_type, fs_home)
    layout = VeryWiseLayout(mesh.hemis, views, mode="interactive")

    plotter = VeryWiseInteractivePlotter(mesh, colormap, surface, layout, 
                                         roi_names=roi_names)

    return plotter.render(colorbar=colorbar, colorbar_label=colorbar_label, title=title,
                          output_html=output_html, cell_px=cell_px, colorbar_width=colorbar_width)
