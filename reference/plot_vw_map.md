# Plot vertex-wise coefficient maps on a 3D cortical surface

Visualize vertex-wise statistical coefficient maps (e.g., LMM fixed
effects) on the FreeSurfer fsaverage cortical surface using Python's
**`nilearn`** and **`plotly`** backends via `{reticulate}`.

Optionally, outlines specific cortical regions of interest (ROIs) using
FreeSurfer annotation information (i.e. 36 regions Desikan-Killiany
atlas).

## Usage

``` r
plot_vw_map(
  res_dir,
  term,
  hemi = c("lh", "rh"),
  measure = "area",
  surface = c("inflated", "pial"),
  outline_rois = NULL,
  fs_home = Sys.getenv("FREESURFER_HOME"),
  show_in_browser = TRUE
)
```

## Arguments

- res_dir:

  Character. Path to the directory containing Freesurfer-style
  vertex-wise result files (`*.mgh`).

- term:

  Character. Name of the model term to visualize (matches entries in
  `stack_names.txt`).

- hemi:

  Character. Hemisphere to plot (`"lh"` or `"rh"`).

- measure:

  Character. Surface measure, e.g. `'area'`, `'thickness'`, `'volume'`.
  Defaults to `'area'`.

- surface:

  Character. Surface mesh to plot on (`"inflated"` or `"pial"`). Deaults
  to `"inflated"`.

- outline_rois:

  Optional character vector. ROI names to outline on the surface plot.
  You can use call
  [`verywise::locate_roi()`](https://seredef.github.io/verywise/reference/locate_roi.md)
  to see the names of all available ROIs.

- fs_home:

  Character. Path to the Freesurfer installation (`$FREESURFER_HOME`).

- show_in_browser:

  Logical (default = TRUE). Open the Figure in browser. This can be
  invoked manually as well using `$show()`.

## Value

A `plotly.graph_objects.Figure` (Python object) that can be displayed or
further modified via `reticulate`.

## Details

The function:

1.  Locates the appropriate MGH files for a given model term and
    hemisphere.

2.  Loads them with **nibabel** via `{reticulate}`.

3.  Displays the coefficient map interactively on the fsaverage inflated
    surface using **`nilearn`**’s `plot_surf_stat_map()` (Plotly
    engine).

4.  Optionally outlines regions from the Freesurfer annotation (via
    `verywise`).

Python dependencies (`nibabel`, `nilearn`, `matplotlib`, `plotly`,
`numpy`) must be installed in the active `{reticulate}` environment.

## Examples

``` r
if (FALSE) { # rlang::is_installed("reticulate") && dir.exists("~/results/fs_results")
# Example usage (requires FreeSurfer fsaverage and precomputed MGH maps)
fs_home <- Sys.getenv("FREESURFER_HOME")
surf_fig <- plot_vw_map(
  res_dir = "~/results/fs_results",
  term = "age",
  hemi = "lh",
  measure = "area",
  outline_rois = c("entorhinal", "precuneus"),
  fs_home = fs_home)
}
```
