# Plot vertex-wise coefficient maps on a 3D cortical surface

Locates the vertex-wise coefficient MGH files for a given model term and
surface measure, optionally applies cluster-wise significance (CWS)
masking, and renders the result on a standard fsaverage surface via
[`plot_vw_surf`](https://seredef.github.io/verywise/reference/plot_vw_surf.md).

Surface maps are loaded directly in R using
[`load.mgh()`](https://seredef.github.io/verywise/reference/load.mgh.md)
without any Python dependency at the file-loading stage. The
Python/nilearn rendering backend is invoked only through
[`plot_vw_surf`](https://seredef.github.io/verywise/reference/plot_vw_surf.md).

## Usage

``` r
plot_vw_map(
  res_dir,
  term,
  measure = "area",
  hemi = c("both", "lh", "rh"),
  surface = c("pial", "inflated"),
  threshold = "cws",
  ...
)
```

## Arguments

- res_dir:

  Character. Path to the directory containing FreeSurfer-style
  vertex-wise result files (`*.mgh`) and `stack_names.txt`.

- term:

  Character. Name of the model term to visualize (matched against
  entries in `stack_names.txt`).

- measure:

  Character. Surface measure to load, e.g. `"area"`, `"thickness"`,
  `"volume"`. Default `"area"`.

- hemi:

  Character. Which hemisphere(s) to plot: `"both"` (default), `"lh"`, or
  `"rh"`.

- surface:

  Character. Surface mesh: `"pial"` (default) or `"inflated"`.

- threshold:

  Controls vertex-level masking before plotting:

  `"cws"` (default)

  :   Cluster-wise significance masking. Loads the matching
      `*.cache.*.sig.ocn.mgh` file and sets all vertices not belonging
      to a significant cluster (OCN label `== 0`) to `NA`. If no OCN
      file is found, the unmasked coefficients are plotted with a
      warning.

  Numeric

  :   Passed directly to
      [`plot_vw_surf`](https://seredef.github.io/verywise/reference/plot_vw_surf.md)
      as an absolute-value threshold (vertices with
      `|value| < threshold` are hidden).

  `NULL`

  :   No masking; all vertices are rendered.

- ...:

  Additional arguments forwarded to
  [`plot_vw_surf`](https://seredef.github.io/verywise/reference/plot_vw_surf.md),
  e.g. `views`, `cmap`, `vmin`, `vmax`, `colorbar`, `colorbar_label`,
  `title`, `to_file`, `dpi`, `fs_home`, `fs_template`.

## Value

Invisibly: the output of
[`plot_vw_surf`](https://seredef.github.io/verywise/reference/plot_vw_surf.md)
— the temp HTML file path (interactive mode) or `to_file` path (static
PNG mode). Called primarily for its side-effect of opening or saving the
figure.

## See also

[`plot_vw_surf`](https://seredef.github.io/verywise/reference/plot_vw_surf.md),
[`plot_vw_diff`](https://seredef.github.io/verywise/reference/plot_vw_diff.md)

## Examples

``` r
if (FALSE) { # rlang::is_installed("reticulate") && dir.exists("~/results/fs_results")
# Both hemispheres, interctive plot, cluster-wise masking (default)
plot_vw_map(
  res_dir = "~/results/fs_results",
  term    = "age",
  measure = "thickness"
)

# Left hemisphere only, numeric threshold, save to PNG
plot_vw_map(
  res_dir   = "~/results/fs_results",
  term      = "age",
  hemi      = "lh",
  threshold = 0.05,
  to_file   = "figures/age_lh.png",
  dpi       = 300L
)

# No masking, custom colour map and title
plot_vw_map(
  res_dir   = "~/results/fs_results",
  term      = "sex",
  threshold = NULL,
  cmap      = "RdBu_r",
  title     = "Sex difference in cortical area"
)
}
```
