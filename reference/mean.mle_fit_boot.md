# Mean of bootstrap replicates.

Returns the empirical mean of the bootstrap replicates (which includes
bootstrap bias). For the bias-corrected point estimate, use
[`params()`](https://queelius.github.io/algebraic.dist/reference/params.html).

## Usage

``` r
# S3 method for class 'mle_fit_boot'
mean(x, ...)
```

## Arguments

- x:

  An `mle_fit_boot` object.

- ...:

  Additional arguments (not used).

## Value

Numeric vector of column means of bootstrap replicates.

## Examples

``` r
set.seed(1)
b <- boot::boot(rexp(50, 2), function(d, i) 1/mean(d[i]), R = 99)
fit <- mle_boot(b)
mean(fit)
#> [1] 2.055982
```
