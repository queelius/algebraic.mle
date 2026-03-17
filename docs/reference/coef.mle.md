# Extract coefficients from an \`mle\` object.

Provides compatibility with the standard R
[`coef()`](https://rdrr.io/r/stats/coef.html) generic. Returns the same
values as `params(x)`.

## Usage

``` r
# S3 method for class 'mle'
coef(object, ...)
```

## Arguments

- object:

  the \`mle\` object

- ...:

  additional arguments (not used)

## Value

Named numeric vector of parameter estimates.
