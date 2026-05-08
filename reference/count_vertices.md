# Count Vertices for a FreeSurfer Template Surface

Returns the number of vertices per hemisphere for a given FreeSurfer
fsaverage surface template. Templates are icosphere subdivisions of the
fsaverage standard surface, with vertex counts decreasing at lower
resolutions.

## Usage

``` r
count_vertices(fs_template)
```

## Arguments

- fs_template:

  A character string specifying the FreeSurfer template. Must be one of:

  `"fsaverage"`

  :   Full-resolution standard template (163,842 vertices per
      hemisphere).

  `"fsaverage6"`

  :   Medium-resolution template (40,962 vertices per hemisphere).

  `"fsaverage5"`

  :   Low-resolution template (10,242 vertices per hemisphere).

  `"fsaverage4"`

  :   Very low-resolution template (2,562 vertices per hemisphere).

  `"fsaverage3"`

  :   Minimal-resolution template (642 vertices per hemisphere).

  `"fsmicro"`

  :   Micro test template (100 vertices per hemisphere).

## Value

An integer giving the number of vertices per hemisphere for the
specified template. Returns `NULL` invisibly if `fs_template` does not
match any known template.

## Examples

``` r
count_vertices("fsaverage")   # 163842
#> [1] 163842
count_vertices("fsaverage5")  # 10242
#> [1] 10242
count_vertices("fsmicro")     # 100
#> [1] 100
```
