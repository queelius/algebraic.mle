# Generic method for computing the predictive confidence interval given an estimator object \`x\`.

Generic method for computing the predictive confidence interval given an
estimator object \`x\`.

## Usage

``` r
pred(x, samp = NULL, alpha = 0.05, ...)
```

## Arguments

- x:

  the estimator object

- samp:

  a sampler for random variable that is parameterized by mle \`x\`

- alpha:

  (1-alpha)/2 confidence interval

- ...:

  additional arguments to pass

## Value

Matrix of predictive confidence intervals.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5), sigma = matrix(0.1), nobs = 100L)
# \donttest{
pred(fit, samp = function(n, theta) rnorm(n, theta[1], 1))
#>          mean    lower    upper
#> [1,] 5.000321 2.962705 7.027373
# }
```
