# List sub-directories till depth `n`

Small recursive utility to list directories given complex folder
structure. This is used in
[`build_supersubject`](https://seredef.github.io/verywise/reference/build_supersubject.md)
for listing all sub-directories until the correct level.

## Usage

``` r
list.dirs.till(path, n)
```

## Arguments

- path:

  : the path where to begin search

- n:

  : the depth of the recursive process (i.e., how many levels of
  sub-directories to list)

## Value

A list of directory/sub-directory names.

## Author

Serena Defina, 2024.
