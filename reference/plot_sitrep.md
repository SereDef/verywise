# Print a plotting environment sit-rep

Runs a full diagnostic of the Python/Kaleido/Xvfb environment used for
static brain map export. Useful for debugging `Error 525` and related
Kaleido/Chromium issues on HPC clusters.

## Usage

``` r
plot_sitrep()
```

## Value

Invisible `NULL`. Output is printed to the console.

## Details

Mirrors the same Python setup as
[`plot_vw_surf`](https://seredef.github.io/verywise/reference/plot_vw_surf.md):
calls `py_require()` to ensure the uv environment is ready, then
`.vw_surf_init_py()` to load `patch_kaleido.py`, and finally runs the
Python sitrep. If Python cannot be initialised at all, a partial R-side
sitrep is printed instead.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_sitrep()
} # }
```
