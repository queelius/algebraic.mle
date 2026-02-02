# Determine if an object \`x\` is an \`mle\` object.

Determine if an object \`x\` is an \`mle\` object.

## Usage

``` r
is_mle(x)
```

## Arguments

- x:

  the object to test

## Value

Logical TRUE if `x` is an `mle` object, FALSE otherwise.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5), sigma = matrix(0.1))
is_mle(fit)       # TRUE
#> [1] TRUE
is_mle(list(a=1)) # FALSE
#> [1] FALSE
```
