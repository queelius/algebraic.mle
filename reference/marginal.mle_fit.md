# Method for obtaining the marginal distribution of an MLE that is based on asymptotic assumptions:

\`x ~ MVN(params(x), inv(H)(x))\`

## Usage

``` r
# S3 method for class 'mle_fit'
marginal(x, indices)
```

## Arguments

- x:

  The distribution object.

- indices:

  The indices of the marginal distribution to obtain.

## Value

An `mle_fit` object representing the marginal distribution for the
selected parameter indices.

## Details

where H is the (observed or expecation) Fisher information matrix.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
marginal(fit, 1)
#> Maximum likelihood estimator of type mle_fit is normally distributed.
#> The estimates of the parameters are given by:
#> mu 
#>  5 
#> The standard error is  0.2 .
#> The asymptotic 95% confidence interval of the parameters are given by:
#>        2.5%    97.5%
#> mu 4.608007 5.391993
#> The MSE of the estimator is  0.04 .
```
