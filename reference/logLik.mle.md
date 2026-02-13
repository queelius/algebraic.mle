# Extract log-likelihood from an \`mle\` object.

Provides compatibility with the standard R
[`logLik()`](https://rdrr.io/r/stats/logLik.html) generic. The returned
object has class `"logLik"` with attributes `df` (number of estimated
parameters) and `nobs` (number of observations), enabling automatic
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) support.

## Usage

``` r
# S3 method for class 'mle'
logLik(object, ...)
```

## Arguments

- object:

  the \`mle\` object

- ...:

  additional arguments (not used)

## Value

An object of class `"logLik"`, or `NULL` if the log-likelihood is not
available.
