# Extract coefficients from an `mle_fit` object.

Delegates to
[`params()`](https://queelius.github.io/algebraic.dist/reference/params.html).

## Usage

``` r
# S3 method for class 'mle_fit'
coef(object, ...)
```

## Arguments

- object:

  the `mle_fit` object

- ...:

  additional arguments (not used)

## Value

Named numeric vector of parameter estimates.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
coef(fit)
#>     mu sigma2 
#>      5      4 
```
