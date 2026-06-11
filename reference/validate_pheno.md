# Validate a user-supplied phenotype data frame

Checks that `pheno` is a `data.frame` containing all columns required by
`roi_associations` + a `folder_id` column. If `pheno = NULL` and
`roi_associations` is empty, a minimal bookkeeping-only phenotype is
built from `data_structure`.

## Usage

``` r
validate_pheno(pheno, data_structure, roi_associations)
```

## Arguments

- pheno:

  A `data.frame` in long format, or `NULL`.

- data_structure:

  Named list passed to
  [`build_minimal_pheno()`](https://seredef.github.io/verywise/reference/build_minimal_pheno.md)
  when `pheno = NULL`.

- roi_associations:

  Named list of numeric vectors as accepted by
  [`simulate_freesurfer_data()`](https://seredef.github.io/verywise/reference/simulate_freesurfer_data.md).

## Value

The validated (or freshly built) `data.frame`.
