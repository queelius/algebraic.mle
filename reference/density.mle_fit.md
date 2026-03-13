# PDF of the asymptotic distribution of an MLE.

Returns a closure computing the probability density function of the
asymptotic normal (univariate) or multivariate normal (multivariate)
distribution of the MLE.

## Usage

``` r
# S3 method for class 'mle_fit'
density(x, ...)
```

## Arguments

- x:

  An `mle_fit` object.

- ...:

  Additional arguments (not used).

## Value

A function computing the PDF at given points.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5), sigma = matrix(0.1))
f <- density(fit)
f(5)
#> [1] 1.261566
```
