# Dimension (number of parameters) of an MLE.

Dimension (number of parameters) of an MLE.

## Usage

``` r
# S3 method for class 'mle_fit'
dim(x)
```

## Arguments

- x:

  An `mle_fit` object.

## Value

Integer; the number of parameters.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
dim(fit)
#> [1] 2
```
