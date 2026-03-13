# Computes the score of an \`mle_fit\` object (score evaluated at the MLE).

If reguarlity conditions are satisfied, it should be zero (or
approximately, if rounding errors occur).

## Usage

``` r
# S3 method for class 'mle_fit'
score_val(x, ...)
```

## Arguments

- x:

  the \`mle_fit\` object to compute the score of.

- ...:

  additional arguments to pass (not used)

## Value

The score vector evaluated at the MLE, or NULL if not available.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), score = c(0.001, -0.002),
  nobs = 100L)
score_val(fit)
#> [1]  0.001 -0.002
```
