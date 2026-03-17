# Build a block-diagonal matrix from a list of square matrices.

Build a block-diagonal matrix from a list of square matrices.

## Usage

``` r
block_diag(blocks, dims)
```

## Arguments

- blocks:

  List of matrices (scalars are promoted to 1x1 matrices).

- dims:

  Integer vector of block dimensions.

## Value

A square matrix with blocks placed along the diagonal.
