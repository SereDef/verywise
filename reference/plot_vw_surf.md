# Plot vertex-wise maps on brain surfaces

Renders left and/or right hemisphere vertex-wise scalar maps on a
standard fsaverage surface using Python/nilearn. Fully headless (no
XQuartz, rgl, or display server required).

Surface meshes are resolved in order:

1.  `FREESURFER_HOME`/subjects/`fs_template`/surf/ - used directly if
    the directory exists (no download needed).

2.  nilearn automatic download, cached permanently in `~/nilearn_data/`
    (download only happens once per template).

**Interactive** (`to_file = NULL`): generates a self-contained WebGL
HTML file opened in the IDE Viewer or default browser.  
**Static** (`to_file` supplied): saves a tiled PNG via matplotlib.

## Usage

``` r
plot_vw_surf(
  lh = NULL,
  rh = NULL,
  lh_mask = NULL,
  rh_mask = NULL,
  vmin = NULL,
  vmax = NULL,
  views = "all",
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
  cell_px = NULL
)
```

## Arguments

- lh, rh:

  Numeric vector,
  `.mgh`/``` .gii`` file path, or  ```NULL`. At least one of `lh`/`rh“
  must be supplied.

- lh_mask, rh_mask:

  Masking boolean maps, or `NULL` (no masking).

- vmin, vmax:

  Numeric colormap limits, or `NULL` for automatic scaling (minimum and
  maximum of the thresholded data).

- views:

  Character vector of camera angles - any subset of `"lateral"`,
  `"medial"`, `"dorsal"`, `"ventral"`, `"anterior"`, `"posterior"`.
  Default: `"all"`.

- surface:

  Surface mesh: `"pial"` (default) or `"inflated"`.

- bg_map:

  Background shading: `"sulc"` (default), `"curv"`, or `"none"`.

- fs_template:

  Resolution (fsaverage template). Must match the length of `lh`/`rh`.
  Default: `"fsaverage"`.

- fs_home:

  (optional) location of FreeSurfer home for templates.

- colorbar:

  Logical. Draw a shared colour bar? Default `TRUE`.

- colorbar_label:

  Character label for the colour bar axis, or `NULL`.

- colorbar_width:

  Fraction of one brain panel's width reserved for the colorbar/density
  strip, or `NULL` for the mode default (static: 0.40, interactive:
  0.20).

- cmap:

  Matplotlib colormap name.

- roi_outline:

  Character vector of DK/aparc region names to outline with contour
  lines (e.g. `c("superiorfrontal", "precentral")`), or `NULL` for no
  ROI overlay.

- title:

  Character figure title, or `NULL`.

- to_file:

  Path ending in `.png` for static export, or `NULL` for interactive
  HTML mode.

- dpi:

  Integer output resolution for static PNG. Default `150`.

- cell_px:

  Brain panel size in pixels: a single number, a length-2 vector
  `c(width, height)`, or `NULL` for the renderer's default (static:
  400x440, interactive: 500x350).

## Value

Invisibly: `to_file` path (static) or the temp HTML file path
(interactive). Called primarily for the side-effect.

## Examples

``` r
if (FALSE) { # \dontrun{
# interactive (opens in Viewer / browser) with file 
plot_vw_surf(lh = 'path/to/lh.coef.mgh', fs_template = "fsaverage5")

# static PNG (4 views, 300 dpi...) with vector objects
plot_vw_surf(
  lh = lh_coef,
  rh = rh_coef,
  threshold = 0.05,
  views = c("lateral", "medial", "dorsal", "ventral"),
  cmap = "RdBu_r",
  title = "Effect of age on cortical thickness",
  to_file = "figures/age_thickness.png",
  dpi = 300
)
} # }
```
