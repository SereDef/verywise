import numpy as np
import sys
import time
from pathlib import Path

# sys.path.insert(0, str(Path(__file__).parent))
#sys.path.insert(0, str(Path().parent)) # interactive
pkg_python_dir = Path.home() / "Desktop" / "Packages" / "verywise" / "inst" / "python"
sys.path.insert(0, str(pkg_python_dir))

import patch_kaleido  # must run BEFORE plotly/kaleido render calls

from plot_vw_surf import vw_surf_static_plotly, vw_surf_interactive 

fs_template = 'fsaverage4'
pattern = 'gradient'
seed = 3108

def make_synthetic(n_vertices, pattern, seed=0):
    """Generate a plausible-looking vertex-wise map for visual testing."""
    rng = np.random.default_rng(seed)

    if pattern == "gradient":
        base = np.linspace(-3, 3, n_vertices)
        rng.shuffle(base)
        data = base + rng.normal(0, 0.4, n_vertices)

    elif pattern == "sparse_clusters":
        data = np.full(n_vertices, np.nan)
        for _ in range(5):
            idx = rng.integers(0, n_vertices)
            size = rng.integers(50, 400)
            span = slice(max(0, idx - size // 2), min(n_vertices, idx + size // 2))
            sign = rng.choice([-1, 1])
            data[span] = sign * rng.uniform(0.5, 3.0)

    elif pattern == "discrete":
        data = rng.integers(1, 8, n_vertices).astype(float)

    elif pattern == "all_positive":
        data = rng.uniform(0.1, 5.0, n_vertices)

    else:
        raise ValueError(f"unknown pattern: {pattern}")

    return data.astype(np.float32)

_FS_TEMPLATE_COUNTS = {
    "fsaverage":  163842,
    "fsaverage6": 40962,
    "fsaverage5": 10242,
    "fsaverage4": 2562,
    "fsaverage3": 642,}

n_vertices = _FS_TEMPLATE_COUNTS[fs_template]
lh = make_synthetic(n_vertices, pattern, seed=seed)
rh = make_synthetic(n_vertices, pattern, seed=seed + 1)

surface = "inflated"
fs_home = "/Applications/freesurfer/7.4.1"

lh_mask = ~((lh > -1) & (lh < 1))
rh_mask = ~((rh > -1) & (rh < 1))

t0 = time.time()

fig = vw_surf_interactive(lh=lh, rh=rh,
        surface=surface,
        # bg_map_type="sulc",
        # vmin=-4, vmax=4,
        lh_mask = lh_mask, 
        rh_mask = rh_mask,
        roi_names=["superiorfrontal", "rostralmiddlefrontal", "precentral", "posteriorcingulate"],
        colorbar=True,
        cmap=None,
        output_html= pkg_python_dir / 'test.html',
        colorbar_label='Betas',
        fs_template=fs_template,
        fs_home=fs_home,
    )


out_path = vw_surf_static_plotly(lh=lh, rh=rh,
        surface=surface,
        # views='all',
        # bg_map_type="sulc",
        vmin=-4, vmax=4.5,
        lh_mask = lh_mask, 
        rh_mask = rh_mask,
        # threshold=threshold,
        roi_names=["superiorfrontal", "rostralmiddlefrontal", "precentral", "posteriorcingulate"],
        colorbar=True,
        cmap=None,
        colorbar_label='Betas',
        title=f"[{pattern}] {fs_template} / {surface}",
        output_file= pkg_python_dir / 'test.png',
        dpi=150,
        fs_template=fs_template,
        fs_home=fs_home,
    )
elapsed = time.time() - t0
elapsed


# ──────────────────────────────────────────────────────────────────────────────
# --- Other (quick viz helper to get positions) -----------------------------

# def vw_dump_camera_preset(hemi, view, surface, fs_template, fs_home=None,
#                           lh=None, rh=None, bg_map_type="sulc"):

#     import tempfile
#     import plotly.graph_objects as go
#     from nilearn.surface import load_surf_mesh

#     lh_arr = _to_array(lh)
#     rh_arr = _to_array(rh)

#     if lh_arr is None and rh_arr is None:
#         from nilearn.datasets import fetch_surf_fsaverage
#         fsavg    = fetch_surf_fsaverage(mesh=fs_template)
#         surf_key = "infl" if surface == "inflated" else surface
#         coords, _ = load_surf_mesh(getattr(fsavg, f"{surf_key}_left"))
#         n      = coords.shape[0]
#         lh_arr = np.zeros(n, dtype=np.float32)
#         rh_arr = np.zeros(n, dtype=np.float32)

#     surf_img, bg_img, hemis, _ = _build_surface_images(
#         lh_arr, rh_arr, fs_template, surface, bg_map_type, fs_home)

#     CANONICAL = {"left": "lateral", "right": "lateral", "both": "dorsal"}
#     fp  = plot_surf(
#         surf_mesh=None, surf_map=surf_img,
#         hemi=hemi, view=CANONICAL[hemi],
#         cmap='viridis', colorbar=False, engine="plotly",
#         bg_map=bg_img, bg_on_data=True,   # ← sulc shading visible
#         alpha=0.9,
#     )
#     raw = _unwrap_plotly(fp)

#     all_coords = np.vstack([
#         load_surf_mesh(p)[0] for p in surf_img.mesh.parts.values()
#     ])
#     mid  = (all_coords.max(axis=0) + all_coords.min(axis=0)) / 2.0
#     half = (all_coords.max(axis=0) - all_coords.min(axis=0)).max() / 2.0 * 1.05
#     def _ax(i):
#         return dict(visible=False, showgrid=False, zeroline=False,
#                               range=[float(mid[i]-half), float(mid[i]+half)])

#     try:
#         cam = _build_camera(hemi, view, surface=surface)
#     except (ValueError, KeyError):
#         cam = {}

#     fig = go.Figure(data=list(raw.data))
#     fig.update_layout(
#         scene=dict(
#             aspectmode="cube",
#             xaxis=_ax(0), yaxis=_ax(1), zaxis=_ax(2),
#             bgcolor="white", camera=cam,
#         ),
#         width=700, height=700, paper_bgcolor="white",
#         margin=dict(l=10, r=10, t=80, b=10),
#         title=dict(
#             text=f"surface=<b>{surface}</b>  hemi=<b>{hemi}</b>  view=<b>{view}</b>",
#             font=dict(size=13),
#         ),
#     )

#     # Inject JS that shows live camera coords as a copyable overlay
#     live_js = """
#         <div id="cam-box" style="
#             position:fixed; top:12px; right:12px; z-index:9999;
#             background:rgba(20,20,20,0.82); color:#e8e8e8;
#             font:13px/1.6 monospace; padding:12px 16px; border-radius:8px;
#             min-width:340px; white-space:pre; user-select:all;
#             box-shadow:0 2px 12px rgba(0,0,0,0.4);">
#         Rotate or pan the brain to update…
#         </div>
#         <div style="position:fixed;bottom:12px;right:12px;z-index:9999;
#             font:11px monospace;color:#888;">
#             Click box → Ctrl/Cmd+C to copy
#         </div>
#         <script>
#         (function() {
#             function fmt(v) { return (v||0).toFixed(4); }
#             function update() {
#                 var el = document.getElementsByClassName('js-plotly-plot')[0];
#                 if (!el) { setTimeout(update, 300); return; }
#                 el.on('plotly_relayout', function(e) {
#                     var cam = el._fullLayout.scene.camera;
#                     if (!cam) return;
#                     var eye = cam.eye, up = cam.up, ctr = cam.center || {x:0,y:0,z:0};
#                     var txt =
#                         ' eye: (' + fmt(eye.x) + ', ' + fmt(eye.y) + ', ' + fmt(eye.z) + ')\n' +
#                         '  up: (' + fmt(up.x)  + ', ' + fmt(up.y)  + ', ' + fmt(up.z)  + ')\n' +
#                         ' ctr: (' + fmt(ctr.x) + ', ' + fmt(ctr.y) + ', ' + fmt(ctr.z) + ')\n\n' +
#                         '# Paste into _CAMERA_PRESETS:\n' +
#                         '("' + '""" + hemi + """' + '", "' + '""" + view + """' + '"): dict(\n' +
#                         '    eye=(' + fmt(eye.x) + ', ' + fmt(eye.y) + ', ' + fmt(eye.z) + '),\n' +
#                         '    up=('  + fmt(up.x)  + ', ' + fmt(up.y)  + ', ' + fmt(up.z)  + '),\n' +
#                         '    center=(' + fmt(ctr.x) + ', ' + fmt(ctr.y) + ', ' + fmt(ctr.z) + '),\n' +
#                         '    distance=' + fmt(Math.sqrt(eye.x**2+eye.y**2+eye.z**2)) + ',\n' +
#                         '),';
#                     document.getElementById('cam-box').innerText = txt;
#                 });
#             }
#             window.addEventListener('load', function() { setTimeout(update, 800); });
#         })();
#         </script>
#         """

#     # Write HTML manually so we can inject the overlay
#     base_html = fig.to_html(include_plotlyjs="cdn", full_html=True)
#     # Inject before </body>
#     html_out  = base_html.replace("</body>", live_js + "\n</body>")

#     f = tempfile.NamedTemporaryFile(suffix=".html", delete=False,
#                                     mode="w", encoding="utf-8")
#     f.write(html_out)
#     f.close()

#     print(f"\nPath: {f.name}")
#     print("Open with:  utils::browseURL(reticulate::py$vw_dump_camera_preset(...))")
#     return f.name
