# @title Save logical cortex mask from FreeSurfer cortex.label files

@title Save logical cortex mask from FreeSurfer cortex.label files

## Usage

``` r
create_cortex_mask(
  freesurfer_home = "",
  fsaverage_template = "fsaverage",
  hemi = "lh",
  outp_path = getwd()
)
```

## Arguments

- freesurfer_home:

  : path to freesurfer home

- fsaverage_template:

  : one of "fsaverage" (163842 vertices), "fsaverage6" (40962),
  "fsaverage5" (10242), "fsaverage4" (2562), "fsaverage3" (642)

- hemi:

  : "lh" of "rh"

- outp_path:

  : path where to save the logical mask
