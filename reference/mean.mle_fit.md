# Mean of the asymptotic distribution of an MLE.

Returns the parameter estimates, which are the mean of the asymptotic
distribution.

## Usage

``` r
# S3 method for class 'mle_fit'
mean(x, ...)
```

## Arguments

- x:

  An `mle_fit` object.

- ...:

  Additional arguments (not used).

## Value

Numeric vector of parameter estimates.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
mean(fit)
#>     mu sigma2 
#>      5      4 
```
