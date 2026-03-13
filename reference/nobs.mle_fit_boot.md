# Method for obtaining the number of observations in the sample used by an \`mle_fit_boot\`.

Method for obtaining the number of observations in the sample used by an
\`mle_fit_boot\`.

## Usage

``` r
# S3 method for class 'mle_fit_boot'
nobs(object, ...)
```

## Arguments

- object:

  the \`mle_fit_boot\` object to obtain the number of observations for

- ...:

  additional arguments to pass (not used)

## Value

Integer number of observations in the original sample.

## Examples

``` r
set.seed(1)
b <- boot::boot(rexp(50, 2), function(d, i) 1/mean(d[i]), R = 99)
fit <- mle_boot(b)
nobs(fit)
#> [1] 50
```
