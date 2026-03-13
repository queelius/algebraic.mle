# Computes the variance-covariance matrix of \`boot\` object. Note: This impelements the \`vcov\` method defined in \`stats\`.

Computes the variance-covariance matrix of \`boot\` object. Note: This
impelements the \`vcov\` method defined in \`stats\`.

## Usage

``` r
# S3 method for class 'mle_fit_boot'
vcov(object, ...)
```

## Arguments

- object:

  the \`boot\` object to obtain the variance-covariance of

- ...:

  additional arguments to pass into \`stats::cov\`

## Value

The variance-covariance matrix estimated from bootstrap replicates.

## Examples

``` r
set.seed(1)
b <- boot::boot(rexp(50, 2), function(d, i) 1/mean(d[i]), R = 99)
fit <- mle_boot(b)
vcov(fit)
#>            [,1]
#> [1,] 0.08252172
```
