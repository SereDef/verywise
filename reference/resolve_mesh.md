# Resolve FreeSurfer surface mesh location

Determines whether to use a local FreeSurfer installation for surface
meshes or fall back to `nilearn`'s automatic download and caching
mechanism.

## Usage

``` r
resolve_mesh(fs_template, fs_home, verbose = TRUE)
```

## Arguments

- fs_template:

  A character string specifying the template (e.g., `"fsaverage"`).

- fs_home:

  A character string specifying the local `FREESURFER_HOME` path, or
  `NULL` to read from environment variables.

- verbose:

  Logical. Whether to print informational messages. Default `TRUE`.

## Value

The FreeSurfer home path when the template is found locally, or `NULL`
to signal Python to use the `nilearn` download+cache path (via
fetch_surf_fsaverage)
