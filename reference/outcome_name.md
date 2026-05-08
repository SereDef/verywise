# Format a Human-Readable Outcome Label

Constructs a descriptive label for a FreeSurfer surface-based
morphometric outcome by combining hemisphere and measure names. Intended
for use in plot titles, model output tables, and report generation.

## Usage

``` r
outcome_name(hemi, measure, sep = " - ")
```

## Arguments

- hemi:

  A character string specifying the hemisphere. Must be one of `"lh"`
  (left hemisphere) or `"rh"` (right hemisphere).

- measure:

  A character string specifying the FreeSurfer surface measure. Must be
  one of:

  `"thickness"`

  :   Cortical thickness.

  `"area"`

  :   Cortical surface area (white surface).

  `"area.pial"`

  :   Cortical surface area (pial surface).

  `"curv"`

  :   Mean curvature.

  `"jacobian_white"`

  :   Jacobian determinant of the white surface deformation field,
      reflecting distortion applied during spherical registration
      (`-surfreg`).

  `"pial"`

  :   Pial surface coordinates.

  `"pial_lgi"`

  :   Local gyrification index (pial surface).

  `"sulc"`

  :   Sulcal depth.

  `"volume"`

  :   Cortical gray matter volume.

  `"w_g.pct"`

  :   White/gray matter intensity ratio.

  `"white.H"`

  :   Mean curvature of the white surface.

  `"white.K"`

  :   Gaussian curvature of the white surface.

- sep:

  A character string used to separate the hemisphere label from the
  measure name. Defaults to `" - "`.

## Value

A single character string combining the hemisphere and measure labels,
e.g. `"Left hemisphere - Cortical thickness"`. Returns a string with
`NA` as the measure component if `measure` does not match any known
value (due to [`switch()`](https://rdrr.io/r/base/switch.html) returning
`NULL`, coerced by `paste0`).

## Details

Hemisphere labels map `"lh"` → `"Left"` and `"rh"` → `"Right"`. Measure
labels follow FreeSurfer's naming conventions as produced by
`mris_preproc` and `mrisurf_thickness`. The `sep` argument allows
flexible formatting for different output contexts (e.g., `sep = ": "`
for axis labels, `sep = "\n"` for multi-line plot titles).

## See also

[`count_vertices`](https://seredef.github.io/verywise/reference/count_vertices.md)
for template resolution lookup.

## Examples

``` r
outcome_name("lh", "thickness")
#> [1] "Left hemisphere - Cortical thickness"
# "Left hemisphere - Cortical thickness"

outcome_name("rh", "sulc", sep = ": ")
#> [1] "Right hemisphere: Sulcal depth"
# "Right hemisphere: Sulcal depth"
```
