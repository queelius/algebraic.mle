# Dimension (number of parameters) of a bootstrap MLE.

Dimension (number of parameters) of a bootstrap MLE.

## Usage

``` r
# S3 method for class 'mle_fit_boot'
dim(x)
```

## Arguments

- x:

  An `mle_fit_boot` object.

## Value

Integer; the number of parameters.

## Examples

``` r
set.seed(1)
b <- boot::boot(rexp(50, 2), function(d, i) 1/mean(d[i]), R = 99)
fit <- mle_boot(b)
dim(fit)
#> [1] 1
```
