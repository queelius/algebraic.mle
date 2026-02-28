# Combine independent MLEs for the same parameter.

Given multiple independent MLEs that estimate the same parameter
\\\theta\\, produces an optimally weighted combination using
inverse-variance (Fisher information) weighting.

## Usage

``` r
combine(x, ...)

# S3 method for class 'list'
combine(x, ...)

# S3 method for class 'mle'
combine(x, ...)
```

## Arguments

- x:

  An `mle` object, or a list of `mle` objects.

- ...:

  Additional `mle` objects to combine.

## Value

An `mle` object representing the optimally weighted combination.

## Details

The combined estimator has:

- `theta.hat`: \\(\sum I_i)^{-1} \sum I_i \hat\theta_i\\

- `sigma`: \\(\sum I_i)^{-1}\\

- `info`: \\\sum I_i\\

- `nobs`: sum of individual sample sizes

When the Fisher information matrix is not directly available but the
variance-covariance matrix is, the FIM is computed as `ginv(vcov)`.

For the legacy interface that accepts a list, see
[`mle_weighted`](https://queelius.github.io/algebraic.mle/reference/mle_weighted.md).

## See also

[`mle_weighted`](https://queelius.github.io/algebraic.mle/reference/mle_weighted.md),
[`joint`](https://queelius.github.io/algebraic.mle/reference/joint.md)

## Examples

``` r
# Three independent estimates of the same rate
fit1 <- mle(theta.hat = c(lambda = 2.1), sigma = matrix(0.04), nobs = 50L)
fit2 <- mle(theta.hat = c(lambda = 1.9), sigma = matrix(0.02), nobs = 100L)
fit3 <- mle(theta.hat = c(lambda = 2.0), sigma = matrix(0.03), nobs = 70L)

comb <- combine(fit1, fit2, fit3)
params(comb)
#>   lambda 
#> 1.976923 
se(comb)
#> [1] 0.09607689
```
