# Determine if an object is an \`mle_boot\` object.

Determine if an object is an \`mle_boot\` object.

## Usage

``` r
is_mle_boot(x)
```

## Arguments

- x:

  the object to test

## Value

Logical TRUE if `x` is an `mle_boot` object, FALSE otherwise.

## Examples

``` r
# Create a simple mle object (not bootstrap)
fit_mle <- mle(theta.hat = 5, sigma = matrix(0.1))
is_mle_boot(fit_mle)  # FALSE
#> [1] FALSE

# Bootstrap example would return TRUE
```
