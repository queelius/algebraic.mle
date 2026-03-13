# Function for obtaining an estimate of the standard error of the MLE object \`x\`.

Function for obtaining an estimate of the standard error of the MLE
object \`x\`.

## Usage

``` r
# S3 method for class 'mle_fit'
se(x, se.matrix = FALSE, ...)
```

## Arguments

- x:

  the MLE object

- se.matrix:

  if \`TRUE\`, return the square root of the variance-covariance

- ...:

  additional arguments to pass (not used)

## Value

Vector of standard errors, or matrix if `se.matrix = TRUE`, or NULL if
variance-covariance is not available.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
se(fit)
#> [1] 0.2000000 0.5656854
```
