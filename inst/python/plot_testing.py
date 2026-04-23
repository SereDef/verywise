import numpy as np

from nilearn.surface import SurfaceImage, PolyMesh
from nilearn.plotting import view_surf

from plot_vw_surf import _get_mesh_paths

n_vert = 642 # fs3
fs_template = 'fsaverage3'
surface = 'inflated'
bg_map_type = 'pial'
fs_home = '/Applications/freesurfer/7.4.1'

rng = np.random.default_rng(42)
lh = rng.standard_normal(n_vert) # n draws from N(0, 1)
rh = rng.standard_normal(n_vert) # n draws from N(0, 1)

nan_idx = rng.choice(n_vert, size=n_vert // 5, replace=False)
lh[nan_idx] = np.nan
rh[nan_idx] = np.nan

# b = rng.uniform(0, 1, size = n)      # n draws from U(0, 1)

mlh, mrh, blh, brh = _get_mesh_paths(fs_template, surface, bg_map_type, fs_home)


from nilearn.datasets import load_fsaverage


# Step 3: wrap geometry in PolyMesh
mesh = PolyMesh(left=mlh, right=mrh)

# Step 4: wrap geometry + values in SurfaceImage
# stat_img = SurfaceImage(
#     mesh = mesh,
#     data = {"left": lh, "right": rh}
# )

fsavg = load_fsaverage(fs_template)
stat_img = SurfaceImage(
    mesh = fsavg['inflated'], 
    data = {"left": lh, "right": rh}
)

mesh_dict = {}
if lh is not None: mesh_dict["left"] = fsavg['inflated'].parts['left']
if rh is not None: mesh_dict["right"] = fsavg['inflated'].parts['right']
mesh = PolyMesh(**mesh_dict)

stat_img = SurfaceImage(
    mesh = fsavg['inflated'], 
    data = {"left": lh, "right": None}
)

# bg = PolyMesh(left=blh, right=brh)

mesh_parts = {}
data_parts = {}
if lh is not None:
    mesh_parts["left"] = fsavg['inflated'].parts['left']
    data_parts["left"] = lh
if rh is not None:
    mesh_parts["right"] = fsavg['inflated'].parts['right']
    data_parts["right"] = rh

mesh = PolyMesh(**mesh_parts)
stat_img = SurfaceImage(mesh=mesh, data=data_parts)

# Set hemi accordingly
hemi = "both" if (lh is not None and rh is not None) else \
       ("left" if lh is not None else "right")

html = view_surf(
    surf_mesh = None,
    surf_map = stat_img,   # carries the mesh internally
    darkness=None,
    vmin = -2, vmax = 2, 
    hemi=hemi,
    threshold=0.05)

# html.open_in_browser()
html.save_as_html('/Users/Serena/Desktop/Packages/verywise/try.html')
