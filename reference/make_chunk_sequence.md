# Define chunks of vertices for analyses

Split a vector of vertex indices into chunks for memory-efficient
parallel processing. Each chunk is annotated with attributes that
facilitate progress reporting (chunk index and 25/50/75% milestones).

## Usage

``` r
make_chunk_sequence(iv, chunk_size = 1000)
```

## Arguments

- iv:

  Integer vector of vertex indices to use.

- chunk_size:

  Integer; target size of each chunk (default = 1000).

## Value

A list of integer vectors, each containing vertex indices and progress
marks attributes.

## Author

Serena Defina, 2024.
