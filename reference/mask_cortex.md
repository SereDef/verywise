# Clean out vertices that are not on the cortex

Clean the super-subject matrix of the vertices that have 0 value or are
not located on the cortex, according to FreeSurfer `*cortical.mask.mgh`
file distributed in `inst/extdata`. This function is used in
[`build_supersubject`](https://seredef.github.io/verywise/reference/build_supersubject.md).

## Usage

``` r
mask_cortex(hemi = "lh", fs_template = "fsaverage")
```

## Arguments

- hemi:

  : (default = "lh") hemisphere ("lh" or "rh")

- fs_template:

  : (default = "fsaverage")

## Value

A (large) logical vector for cortical (or technical, i.e. vertex = 0)
mask.
