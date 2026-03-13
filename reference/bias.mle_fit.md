# Computes the bias of an \`mle_fit\` object assuming the large sample approximation is valid and the MLE regularity conditions are satisfied. In this case, the bias is zero (or zero vector).

This is not a good estimate of the bias in general, but it's arguably
better than returning \`NULL\`.

## Usage

``` r
# S3 method for class 'mle_fit'
bias(x, theta = NULL, ...)
```

## Arguments

- x:

  the \`mle_fit\` object to compute the bias of.

- theta:

  true parameter value. normally, unknown (NULL), in which case we
  estimate the bias (say, using bootstrap)

- ...:

  additional arguments to pass

## Value

Numeric vector of zeros (asymptotic bias is zero under regularity
conditions).

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
bias(fit)
#> [1] 0 0
```
