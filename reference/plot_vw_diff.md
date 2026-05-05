# Plot the vertex-wise difference between two scalar maps

Computes the vertex-wise difference (`a - b`) between two surface maps
(numeric vectors or MGH/GII file paths) for each hemisphere, then
renders the result on a standard fsaverage surface via
[`plot_vw_surf`](https://seredef.github.io/verywise/reference/plot_vw_surf.md).

Useful for contrasting two conditions, time-points, groups, or model
terms without pre-computing the difference outside R.

## Usage

``` r
plot_vw_diff(
  lh_a = NULL,
  lh_b = NULL,
  rh_a = NULL,
  rh_b = NULL,
  label_a = "a",
  label_b = "b",
  ...
)
```

## Arguments

- lh_a, lh_b:

  Left-hemisphere map: numeric vector or path to an MGH/GII file. Both
  must have the same length / vertex count. Pass `NULL` to omit the left
  hemisphere entirely.

- rh_a, rh_b:

  Right-hemisphere map: same as above for the right hemisphere. Pass
  `NULL` to omit the right hemisphere.

- label_a, label_b:

  Short character labels used in the default figure title (e.g.
  `"group A"`, `"group B"`). Ignored when `title` is supplied via `...`.

- ...:

  Additional arguments forwarded to
  [`plot_vw_surf`](https://seredef.github.io/verywise/reference/plot_vw_surf.md)
  (e.g. `surface`, `views`, `cmap`, `vmin`, `vmax`, `threshold`,
  `colorbar`, `colorbar_label`, `title`, `to_file`, `dpi`, `fs_home`,
  `fs_template`).

## Value

Invisibly: the output of
[`plot_vw_surf`](https://seredef.github.io/verywise/reference/plot_vw_surf.md)
(temp HTML path or `to_file` path).

## Examples

``` r
if (FALSE) { # rlang::is_installed("reticulate") && dir.exists("path/to/model1/")
# Two numeric vectors
plot_vw_diff(
  lh_a = lh_term1, lh_b = lh_term2,
  rh_a = rh_term1, rh_b = rh_term2,
  label_a = "Term 1", label_b = "Term 2",
  cmap = "RdBu_r"
)

# Two MGH files, left hemisphere only
plot_vw_diff(
  lh_a = "path/to/model1/lh.area.aic.mgh",
  lh_b = "path/to/model2/lh.area.aic.mgh",
  label_a = "model1", label_b = "model2"
)
}
```
