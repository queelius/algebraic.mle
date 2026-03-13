# Computes the variance-covariance matrix of \`mle_fit\` object.

Computes the variance-covariance matrix of \`mle_fit\` object.

## Usage

``` r
# S3 method for class 'mle_fit'
vcov(object, ...)
```

## Arguments

- object:

  the \`mle_fit\` object to obtain the variance-covariance of

- ...:

  additional arguments to pass (not used)

## Value

the variance-covariance matrix

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
vcov(fit)
#>      [,1] [,2]
#> [1,] 0.04 0.00
#> [2,] 0.00 0.32
```
