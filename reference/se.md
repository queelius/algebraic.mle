# Generic method for obtaining the standard errors of an estimator.

Generic method for obtaining the standard errors of an estimator.

## Usage

``` r
se(x, ...)
```

## Arguments

- x:

  the estimator

- ...:

  additional arguments to pass

## Value

Vector of standard errors for each parameter.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
se(fit)
#> [1] 0.2000000 0.5656854
```
