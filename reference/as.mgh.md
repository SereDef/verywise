# Create an object structured like an MGH object

Create an object structured like an MGH object

## Usage

``` r
as.mgh(
  x,
  version = 1L,
  ndim1 = as.integer(length(x)),
  ndim2 = 1L,
  ndim3 = 1L,
  nframes = 1L,
  type = 3L,
  dof = 0L
)
```

## Arguments

- x:

  The vertex-wise values

- version:

  Version (default = 1)

- ndim1:

  Width / 1st dimension (default = length(x))

- ndim2:

  Height / 2nd dimension (default = 1)

- ndim3:

  Depth / 3rd dimension (default = 1)

- nframes:

  Number of scalar components (default = 1)

- type:

  Data type, can be UCHAR (0), SHORT (4), INT (1) or FLOAT (3) (default
  = 3)

- dof:

  Degrees of freedom (default = 0)

## Value

A mgh object with data and various header elements

## Author

Sander Lamballais, 2018.
