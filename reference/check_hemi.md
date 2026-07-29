# Validate and load hemisphere data

Checks if the input is a valid file path or numeric vector. If a `.mgh`
file path is provided, it automatically loads it. Validates that the
number of vertices matches the specified FreeSurfer template.

## Usage

``` r
check_hemi(hemi, fs_template)
```

## Arguments

- hemi:

  A numeric vector, a character string representing a file path (e.g.,
  `.mgh`, `.gii`), or `NULL`.

- fs_template:

  A character string specifying the FreeSurfer template (e.g.,
  `"fsaverage5"`, `"fsaverage"`).

## Value

A numeric vector representing the surface data, or `NULL` if the input
is `NULL`.
