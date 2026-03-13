# Computes the estimate of the MSE of a \`boot\` object.

Computes the estimate of the MSE of a \`boot\` object.

## Usage

``` r
# S3 method for class 'mle_fit_boot'
mse(x, theta = NULL, ...)
```

## Arguments

- x:

  the \`boot\` object to compute the MSE of.

- theta:

  true parameter value (not used for \`mle_boot\`)

- ...:

  pass additional arguments into \`vcov\`

## Value

The MSE matrix estimated from bootstrap variance and bias.

## Examples

``` r
set.seed(1)
b <- boot::boot(rexp(50, 2), function(d, i) 1/mean(d[i]), R = 99)
fit <- mle_boot(b)
mse(fit)
#>            [,1]
#> [1,] 0.08307621
```
