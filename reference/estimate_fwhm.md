# Estimate full-width half maximum (FWHM)

This function is a wrapper for the FreeSurfer command `mris_fwhm`. It
estimates the full-width half maximum (FWHM) of the brain data based on
the residuals of your model.

## Usage

``` r
estimate_fwhm(
  result_path,
  hemi,
  mask = NULL,
  fs_template = "fsaverage",
  verbose = FALSE
)
```

## Arguments

- result_path:

  : path where model residuals are stored and where results are
  generated: this includes two files: "hemi.measure.fwhm.dat" and
  "hemi.measure.finamMask.mgh".

- hemi:

  : hemisphere.

- mask:

  : apply a mask (default = cortical mask + excluding problematic
  vertices)

- fs_template:

  : (default = "fsaverage") data template.

- verbose:

  : (default = FALSE) verbosity.

## Value

An integer FWHM value
