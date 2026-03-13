# Generic method for computing the bias of an estimator object.

Generic method for computing the bias of an estimator object.

## Usage

``` r
bias(x, theta = NULL, ...)
```

## Arguments

- x:

  the object to compute the bias of.

- theta:

  true parameter value. usually, this is unknown (NULL), in which case
  we estimate the bias

- ...:

  pass additional arguments

## Value

The bias of the estimator. The return type depends on the specific
method.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
bias(fit, theta = c(mu = 5, sigma2 = 4))
#> [1] 0 0
```
