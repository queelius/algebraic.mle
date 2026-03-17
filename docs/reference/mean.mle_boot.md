# Mean of bootstrap replicates.

Returns the empirical mean of the bootstrap replicates (which includes
bootstrap bias). For the bias-corrected point estimate, use
[`params()`](https://queelius.github.io/algebraic.dist/reference/params.html).

## Usage

``` r
# S3 method for class 'mle_boot'
mean(x, ...)
```

## Arguments

- x:

  An `mle_boot` object.

- ...:

  Additional arguments (not used).

## Value

Numeric vector of column means of bootstrap replicates.
