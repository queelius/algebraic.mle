# Compose independent MLEs into a joint MLE.

Given two or more independent MLEs with disjoint parameter sets,
produces a joint MLE with block-diagonal variance-covariance structure.

## Usage

``` r
joint(x, ...)

# S3 method for class 'mle'
joint(x, ...)
```

## Arguments

- x:

  An `mle` object.

- ...:

  Additional `mle` objects to join.

## Value

An `mle` object representing the joint MLE.

## Details

The joint MLE has:

- `theta.hat`: concatenation of all parameter vectors

- `sigma`: block-diagonal from individual vcov matrices

- `loglike`: sum of log-likelihoods (when all available)

- `info`: block-diagonal from individual FIMs (when all available)

- `score`: concatenation of score vectors (when all available)

- `nobs`: NULL (different experiments have no shared sample size)

## Examples

``` r
# Two independent experiments
fit_rate <- mle(theta.hat = c(lambda = 2.1), sigma = matrix(0.04), nobs = 50L)
fit_shape <- mle(theta.hat = c(k = 1.5, s = 3.2),
                 sigma = matrix(c(0.1, 0.02, 0.02, 0.3), 2, 2), nobs = 100L)

# Joint MLE: 3 params, block-diagonal covariance
j <- joint(fit_rate, fit_shape)
params(j)   # c(lambda = 2.1, k = 1.5, s = 3.2)
#> lambda      k      s 
#>    2.1    1.5    3.2 
vcov(j)     # 3x3 block-diagonal
#>      [,1] [,2] [,3]
#> [1,] 0.04 0.00 0.00
#> [2,] 0.00 0.10 0.02
#> [3,] 0.00 0.02 0.30

# Existing algebra works on the joint:
marginal(j, 2:3)   # recover shape params
#> Maximum likelihood estimator of type mle is normally distributed.
#> The estimates of the parameters are given by:
#>   k   s 
#> 1.5 3.2 
#> The standard error is  0.3162278 0.5477226 .
#> The asymptotic 95% confidence interval of the parameters are given by:
#>       2.5%    97.5%
#> k 0.880205 2.119795
#> s 2.126484 4.273516
#> The MSE of the individual components in a multivariate estimator is:
#>      [,1] [,2]
#> [1,] 0.10 0.02
#> [2,] 0.02 0.30
as_dist(j)          # MVN for distribution algebra
#> Multivariate normal distribution (3 dimensions) 
#>   mu:
#> [1] 2.1 1.5 3.2
#>   sigma:
#>      [,1] [,2] [,3]
#> [1,] 0.04 0.00 0.00
#> [2,] 0.00 0.10 0.02
#> [3,] 0.00 0.02 0.30
```
