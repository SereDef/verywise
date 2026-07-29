# Validate a logical surface mask

Ensures the provided mask is a logical vector and matches the vertex
resolution of the specified FreeSurfer template.

## Usage

``` r
check_mask(mask, fs_template)
```

## Arguments

- mask:

  A logical vector used for thresholding, or `NULL`.

- fs_template:

  A character string specifying the FreeSurfer template.

## Value

A logical vector, or `NULL` if the input is `NULL`.
