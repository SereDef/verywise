# Load "phenotype" file into R based on its extension

Takes a path to a file and reads it into memory.

## Usage

``` r
load_pheno_file(pheno_str, verbose = TRUE)
```

## Arguments

- pheno_str:

  A character string specifying the file path.

- verbose:

  (default = TRUE).

## Value

A data frame containing the loaded data

## Details

The function checks if the file exists, identifies the file extension,
and loads the data accordingly. Supported file types: .csv, .rds, .sav,
and .txt.
