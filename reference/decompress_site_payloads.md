# Load and harmonise site payload archives

Extracts each site's `.tar.gz` archive, loads the static sufficient
statistics (`.static.rds`) and attaches the memory-mapped FBM objects
(\\X^\top Y\\ and partial sums of \\Y\\), then transposes the resulting
list-of-sites into a list-of-components ready for the aggregation step.

## Usage

``` r
decompress_site_payloads(site_names, input_dir, hemi, measure)
```

## Arguments

- site_names:

  Character vector of site identifiers.

- input_dir:

  Character. Directory containing `<site>.tar.gz` archives.

- hemi:

  Character. Hemisphere (`"lh"` or `"rh"`).

- measure:

  Character. Surface measure (derived from the model formula LHS).

## Value

A named list with one element per payload component. Key elements:

- `n_obs`:

  Named integer vector of per-site sample sizes.

- `terms`:

  Character vector of fixed-effect term names (shared).

- `fs_template`:

  FreeSurfer template name (shared).

- `XtX`:

  List of \\K\\ matrices \\X_k^\top X_k\\ (\\p \times p\\).

- `X1`:

  List of \\K\\ vectors \\X_k^\top \mathbf{1}\_{n_k}\\ (length \\p\\).

- `XtY`:

  List of \\K\\ FBMs \\X_k^\top Y_k\\ (\\p \times V\\).

- `psumsY`:

  List of \\K\\ FBMs with row 1 = \\\mathbf{1}^\top Y_k\\, row 2 =
  \\Y_k^\top Y_k\\ (\\2 \times V\\).

## Details

Components that must be **identical** across sites (e.g. `terms`,
`fs_template`) are validated and reduced to a single value; an error is
thrown if any site disagrees. Scalar components (e.g. `n_obs`) are
unlisted into named vectors. Matrix and FBM components are kept as lists
of length \\K\\.
