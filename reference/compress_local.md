# Compress a site's results into a tar.gz archive

Creates a gzip-compressed tar archive containing a sufficient statistics
from the local site. The archive is written to the output directory and
can be safely shared among insititutions

## Usage

``` r
compress_local(outp_dir, site_name, then_rm_files = TRUE)
```

## Arguments

- outp_dir:

  Character scalar. Directory containing the local output files to be
  packaged (e.g., FBM `.rds` + `.bk` pairs and other `.rds` objects).

- site_name:

  Character scalar. Site identifier used as the archive basename (output
  file will be `{site_name}.tar.gz`).

- then_rm_files:

  Logical. Remove input files after compressing them? Default: TRUE.

## Value

Invisibly returns `NULL`. Called for its side effect of creating the
`{site_name}.tar.gz` file.

## Author

Serena Defina, 2026.

## Examples

``` r
if (FALSE) { # \dontrun{
# Suppose outp_dir contains:
#   XtY.rds, XtY.bk, psumsY.rds, psums.bk, static.rds
compress_local(outp_dir = "path/to/site_outputs", site_name = "siteA")
# Creates "siteA.tar.gz" in the current working directory.
} # }
```
