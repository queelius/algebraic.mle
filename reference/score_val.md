# Generic method for computing the score of an estimator object (gradient of its log-likelihood function evaluated at the MLE).

Generic method for computing the score of an estimator object (gradient
of its log-likelihood function evaluated at the MLE).

## Usage

``` r
score_val(x, ...)
```

## Arguments

- x:

  the object to compute the score of.

- ...:

  pass additional arguments

## Value

The score vector evaluated at the MLE.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), score = c(0.001, -0.002),
  loglike = -120, nobs = 100L)
score_val(fit)
#> [1]  0.001 -0.002
```
