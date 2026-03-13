# Generic method for determining the orthogonal parameters of an estimator.

Generic method for determining the orthogonal parameters of an
estimator.

## Usage

``` r
orthogonal(x, tol, ...)
```

## Arguments

- x:

  the estimator

- tol:

  the tolerance for determining if a number is close enough to zero

- ...:

  additional arguments to pass

## Value

Logical vector or matrix indicating which parameters are orthogonal.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), info = solve(diag(c(0.04, 0.32))),
  loglike = -120, nobs = 100L)
orthogonal(fit)
#>       [,1]  [,2]
#> [1,] FALSE  TRUE
#> [2,]  TRUE FALSE
```
