# Method for obtaining the number of parameters of an \`boot\` object.

Method for obtaining the number of parameters of an \`boot\` object.

## Usage

``` r
# S3 method for class 'mle_fit_boot'
nparams(x)
```

## Arguments

- x:

  the \`boot\` object to obtain the number of parameters of

## Value

Integer number of parameters.

## Examples

``` r
set.seed(1)
b <- boot::boot(rexp(50, 2), function(d, i) 1/mean(d[i]), R = 99)
fit <- mle_boot(b)
nparams(fit)
#> [1] 1
```
