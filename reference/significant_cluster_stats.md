# Calculate Significant Cluster Statistics

This function computes statistics (mean or median) for significant
clusters based on the supersubject data and results directories. It
retrieves the clusters from an MGH file, extracts relevant measurements
for each cluster, and computes either the mean or median for each
subject in the supersubject data.

## Usage

``` r
significant_cluster_stats(stat, ss_dir, res_dir, term, measure, hemi)
```

## Arguments

- stat:

  A character string specifying the statistic to compute. Either "mean"
  or "median".

- ss_dir:

  A character string indicating the directory containing the
  super-subject data matrix.

- res_dir:

  A character string indicating the directory containing the result
  data, including the stack information and the significant clusters.

- term:

  A character string indicating the term or "stack name" to find in the
  stack file.

- measure:

  A character string indicating the measure used (e.g., "thickness",
  "volume").

- hemi:

  A character string indicating the hemisphere. One of "lh" for left
  hemisphere or "rh" for right hemisphere.

## Value

A data frame with the computed statistics for each cluster and subject.
The data frame contains one column for each computed statistic, indexed
by the subject's folder ID.

## Examples

``` r
if (FALSE) { # dir.exists("path/to/ss_dir")
result <- significant_cluster_stats(
  stat = "mean", 
  ss_dir = "path/to/ss_dir", 
  res_dir = "path/to/res_dir", 
  term = "age", 
  measure = "thickness", 
  hemi = "lh"
)
}
```
