# PDF of the asymptotic distribution of an MLE.

Returns a closure computing the probability density function of the
asymptotic normal (univariate) or multivariate normal (multivariate)
distribution of the MLE.

## Usage

``` r
# S3 method for class 'mle'
density(x, ...)
```

## Arguments

- x:

  An `mle` object.

- ...:

  Additional arguments (not used).

## Value

A function computing the PDF at given points.
