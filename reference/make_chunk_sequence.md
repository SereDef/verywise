# Define chunks of vertices for analyses

This function takes a vector of values representing the dimension of the
vertex-wise data to use and splits it into chunks for memory-efficient
parallel processing.

## Usage

``` r
make_chunk_sequence(iv, chunk_size = 1000)
```

## Arguments

- iv:

  : indices of the vertex-wise brain data to use.

- chunk_size:

  : (default = 1000) how big are the chunks

## Value

A list with n_chunks elements. Each element is a vector of vertex
positions.

## Author

Serena Defina, 2024.
