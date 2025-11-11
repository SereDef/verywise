# Subset an existing supersubject matrix by matching folder IDs

This function subsets a
[FBM](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
supersubject matrix that was created using
[`build_supersubject()`](https://seredef.github.io/verywise/reference/build_supersubject.md),
retaining only the rows that correspond to a set of folder (or row) IDs.
It reads row names from an associated `.csv` file, checks for missing
IDs, and writes logs if any are not found. The new subsetted matrix can
be saved for future use.

## Usage

``` r
subset_supersubject(
  supsubj_dir,
  supsubj_file,
  folder_ids,
  error_cutoff = 20,
  new_supsubj_dir,
  n_cores = 1,
  save_rds = FALSE,
  verbose = TRUE
)
```

## Arguments

- supsubj_dir:

  Character string indicating the path to the directory containing the
  supersubject files (i.e the supersubject matrix itself as a `.rds`
  file, and the associated `.bk` and `.rownames.csv` files).

- supsubj_file:

  Character string indicating the name of the supersubject `.rds` file.
  Must follow the naming pattern
  `"<hemi>.<measure>.<fs_template>.supersubject.rds"`.

- folder_ids:

  Character vector of folder IDs to retain in the new ss matrix. This
  should also be a column in the phenotype dataset.

- error_cutoff:

  Integer indicating the maximum number of missing IDs that is allowed
  before the function throws an error. If the number of missing IDs is
  ` <= error_cutoff `, a warning is issued instead. Default: 20.

- new_supsubj_dir:

  Character string indicating the path to the directory where the new
  supersubject files should be stored (either temporarily or permanently
  if ` save_rds == TRUE `. Created if it does not exist.

- n_cores:

  Integer indicating the number of cores to use for parallel processing.

- save_rds:

  Logical. If `TRUE`, the new ss is also saved to a `.rds` file inside
  `new_supsubj_dir`.

- verbose:

  Logical. Default: `TRUE`.

## Value

A [FBM](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
object containing the subsetted supersubject matrix.

## Details

The function performs the following steps:

1.  Reads row names from the supersubject's `.csv` file.

2.  Checks whether all `folder_ids` exist in the supersubject.

3.  Logs missing IDs to `issues.log` in `new_supsubj_dir`.

4.  If the number of missing IDs exceeds `error_cutoff`, stops with an
    error.

5.  Creates a new FBM with only the matching rows, writing it blockwise
    to avoid excessive RAM usage.

6.  Writes the filtered row names to `ss.rownames.csv` in
    `new_supsubj_dir`.

## See also

[build_supersubject](https://seredef.github.io/verywise/reference/build_supersubject.md)
