# Move key result files from one directory to another (for sharing and visualization)

The function searches through a `verywise` results directory for all
files needed for visualization: i.e. clusters, coefficient maps, and
"stack names". It then copies them to a destination directory while
preserving any sub-directory structure in the source results folder.

## Usage

``` r
move_result_files(from_dir, to_dir)
```

## Arguments

- from_dir:

  Character string indicating the path to the "source" results directory
  (i.e. `outp_dir` in the analysis call)

- to_dir:

  Character string indicating the path to the directory where matching
  files will be copied. Note: any required sub-directories will be
  created automatically.

## Value

Invisibly returns a logical vector indicating whether each file was
successfully copied. Side effects: `.mgh` files are copied to `to_dir`.

## Details

Files are matched using these regular expression patterns:

- `"[a-z]{1}h.[a-z]+.stack[0-9]{1,2}.cache.th30.abs.sig.ocn.mgh"`:
  cluster files

- `"[a-z]{1}h.[a-z]+.stack[0-9]{1,2}.coef.mgh"`: coefficient (beta) maps

- `"stack_names.txt"`: stack names

## Examples

``` r
if (FALSE) { # \dontrun{
# Move matching files from "results" to "plots"
move_result_files("path/to/results", "path/to/plots")
} # }
```
