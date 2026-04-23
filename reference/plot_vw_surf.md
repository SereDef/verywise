# Plot vertex-wise maps on brain surfaces

Renders left and/or right hemisphere vertex-wise scalar maps on a
standard fsaverage surface using Python/nilearn. Fully headless — no
XQuartz, rgl, or display server required.

Surface meshes are resolved in order:

1.  `FREESURFER_HOME`/subjects/`fs_template`/surf/ — used directly if
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
  fs_template = "fsaverage",
  fs_home = NULL,
  surface = c("inflated", "pial", "white"),
  cmap = "RdBu_r",
  bg_map = c("sulc", "curv", "none"),
  vmin = NULL,
  vmax = NULL,
  threshold = NULL,
  views = "all",
  colorbar = TRUE,
  colorbar_label = NULL,
  title = NULL,
  to_file = NULL,
  dpi = 150L
)
```

## Arguments

- lh:

  Numeric vector, `.mgh`/`.gii` file path, or `NULL`.

- rh:

  Numeric vector, `.mgh`/`.gii` file path, or `NULL`. At least one of
  `lh`/`rh` must be supplied.

- fs_template:

  fsaverage template: `"fsaverage3"` through `"fsaverage"`. Must match
  the length of lh/rh. Default `"fsaverage"`.

- fs_home:

  (optional) location of FreeSurfer home for templates.

- surface:

  Surface mesh: `"inflated"` (default), `"pial"`, or `"white"`.

- cmap:

  Matplotlib colormap name. Default `"RdBu_r"`.

- bg_map:

  Background shading: `"sulc"` (default), `"curv"`, or `"none"`.

- vmin, vmax:

  Numeric colour limits, or `NULL` for automatic symmetric scaling.

- threshold:

  Absolute-value masking threshold, or `NULL` (no masking).

- views:

  Character vector of camera angles — any subset of `"lateral"`,
  `"medial"`, `"dorsal"`, `"ventral"`, `"anterior"`, `"posterior"`.
  Default: `"all"`.

- colorbar:

  Logical. Draw a shared colour bar? Default `TRUE`.

- colorbar_label:

  Character label for the colour bar axis, or `NULL`.

- title:

  Character figure title, or `NULL`.

- to_file:

  Path ending in `.png` for static export, or `NULL` for interactive
  HTML mode.

- dpi:

  Integer output resolution for static PNG. Default `150L`.

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
  cmap = "RdBu_r",
  threshold = 0.05,
  views = c("lateral", "medial", "dorsal", "ventral"),
  title = "Effect of age on cortical thickness",
  to_file = "figures/age_thickness.png",
  dpi = 300L
)
} # }
```
