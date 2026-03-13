# Method for obtaining the number of observations in the sample used by an \`mle_fit\`.

Method for obtaining the number of observations in the sample used by an
\`mle_fit\`.

## Usage

``` r
# S3 method for class 'mle_fit'
nobs(object, ...)
```

## Arguments

- object:

  the \`mle_fit\` object to obtain the number of observations for

- ...:

  additional arguments to pass (not used)

## Value

Integer number of observations, or NULL if not available.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
nobs(fit)
#> [1] 100
```
