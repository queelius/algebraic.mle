# Log-likelihood of an `mle_fit` object.

Returns a `"logLik"` object with `df` (number of parameters) and `nobs`
attributes, enabling [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) via `stats::AIC.default`.

## Usage

``` r
# S3 method for class 'mle_fit'
logLik(object, ...)
```

## Arguments

- object:

  the `mle_fit` object

- ...:

  additional arguments (not used)

## Value

A `"logLik"` object.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
logLik(fit)
#> 'log Lik.' -120 (df=2)
```
