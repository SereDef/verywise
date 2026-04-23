# Plot vertex-wise maps on brain surfaces

Render a vertex-wise scalar map on an fsaverage template surface using
fsbrain. Works cross-platform without XQuartz:

- **Interactive** (no `output_file`): returns an `rglwidget` that opens
  in the RStudio Viewer, Positron, or any browser. No display server
  required.

- **Static PNG** (`output_file` supplied): renders via
  `webshot2`/headless-Chrome on macOS and other platforms without native
  OpenGL snapshot support; falls back to direct OpenGL snapshot on Linux
  with a display.

## Usage

``` r
plot_vw_surf_r(
  lh = NULL,
  rh = NULL,
  surface = c("inflated", "pial", "white"),
  cmap = grDevices::colorRampPalette(c("blue", "white", "red")),
  bg_map = c("sulc", "curv", "none"),
  symm = TRUE,
  views = "t4",
  colorbar_legend = NULL,
  output_file = NULL,
  fs_home = Sys.getenv("FREESURFER_HOME"),
  fs_template = "fsaverage"
)
```

## Arguments

- lh, rh:

  Numeric vector **or** path to a `.mgh`/`.gii` file, or `NULL`.
  Vertex-wise values for the left / right hemisphere. At least one must
  be non-`NULL`.

- surface:

  `"inflated"` (default), `"pial"`, or `"white"`.

- cmap:

  A `colorRampPalette`-style function. Default: blue–white–red diverging
  palette.

- bg_map:

  Background shading: `"sulc"` (default), `"curv"`, or `"none"`.

- symm:

  Logical. Use a symmetric (zero-centred) colour scale? Default `TRUE`.

- views:

  Camera angles. `"t4"` (4 standard views, default), `"t9"` (9 views),
  `"si"` (superior/inferior), or a character vector of individual angle
  names.

- colorbar_legend:

  Character label for the colour bar, or `NULL`.

- output_file:

  Path ending in `.png` for static export, or `NULL` for interactive
  mode.

- fs_home:

  FreeSurfer installation root. Defaults to `FREESURFER_HOME`.

- fs_template:

  fsaverage template name. Default `"fsaverage"`.

## Value

In interactive mode, an `rglwidget` object (auto-prints as a WebGL
scene). In static mode, the `output_file` path invisibly.
