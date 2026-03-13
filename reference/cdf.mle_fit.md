# CDF of the asymptotic distribution of an MLE.

CDF of the asymptotic distribution of an MLE.

## Usage

``` r
# S3 method for class 'mle_fit'
cdf(x, ...)
```

## Arguments

- x:

  An `mle_fit` object.

- ...:

  Additional arguments (not used).

## Value

A function computing the CDF at given points.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5), sigma = matrix(0.1))
F <- cdf(fit)
F(5)
#> [1] 0.5
```
