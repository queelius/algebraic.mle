# Method for obtaining the parameters of an \`boot\` object.

Method for obtaining the parameters of an \`boot\` object.

## Usage

``` r
# S3 method for class 'mle_fit_boot'
params(x)
```

## Arguments

- x:

  the \`boot\` object to obtain the parameters of.

## Value

Numeric vector of parameter estimates (the original MLE).

## Examples

``` r
set.seed(1)
b <- boot::boot(rexp(50, 2), function(d, i) 1/mean(d[i]), R = 99)
fit <- mle_boot(b)
params(fit)
#> [1] 2.032434
```
