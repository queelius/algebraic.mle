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
