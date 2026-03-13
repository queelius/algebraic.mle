# Function for printing a \`summary\` object for an \`mle_fit\` object.

Function for printing a \`summary\` object for an \`mle_fit\` object.

## Usage

``` r
# S3 method for class 'summary_mle_fit'
print(x, ...)
```

## Arguments

- x:

  the \`summary_mle_fit\` object

- ...:

  pass additional arguments

## Value

Invisibly returns `x`.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
print(summary(fit))
#> Maximum likelihood estimator of type mle_fit is normally distributed.
#> The estimates of the parameters are given by:
#>     mu sigma2 
#>      5      4 
#> The standard error is  0.2 0.5656854 .
#> The asymptotic 95% confidence interval of the parameters are given by:
#>            2.5%    97.5%
#> mu     4.608007 5.391993
#> sigma2 2.891277 5.108723
#> The MSE of the individual components in a multivariate estimator is:
#>      [,1] [,2]
#> [1,] 0.04 0.00
#> [2,] 0.00 0.32
#> The log-likelihood is  -120 .
#> The AIC is  244 .
```
