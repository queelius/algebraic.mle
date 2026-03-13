# Computes the estimate of the bias of a \`mle_fit_boot\` object.

Computes the estimate of the bias of a \`mle_fit_boot\` object.

## Usage

``` r
# S3 method for class 'mle_fit_boot'
bias(x, theta = NULL, ...)
```

## Arguments

- x:

  the \`mle_fit_boot\` object to compute the bias of.

- theta:

  true parameter value (not used for \`mle_fit_boot\`).

- ...:

  pass additional arguments (not used)

## Value

Numeric vector of estimated bias (mean of bootstrap replicates minus
original estimate).

## Examples

``` r
set.seed(1)
b <- boot::boot(rexp(50, 2), function(d, i) 1/mean(d[i]), R = 99)
fit <- mle_boot(b)
bias(fit)
#> [1] 0.0235476
```
