# Quantile function of the asymptotic distribution of an MLE.

Quantile function of the asymptotic distribution of an MLE.

## Usage

``` r
# S3 method for class 'mle_fit'
inv_cdf(x, ...)
```

## Arguments

- x:

  An `mle_fit` object.

- ...:

  Additional arguments (not used).

## Value

A function computing quantiles for given probabilities.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5), sigma = matrix(0.1))
q <- inv_cdf(fit)
q(0.975)
#> [1] 5.619795
```
