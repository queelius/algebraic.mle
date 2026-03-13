# Method for obtaining the number of parameters of an \`mle_fit\` object.

Method for obtaining the number of parameters of an \`mle_fit\` object.

## Usage

``` r
# S3 method for class 'mle_fit'
nparams(x)
```

## Arguments

- x:

  the \`mle_fit\` object to obtain the number of parameters of

## Value

Integer number of parameters.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
nparams(fit)
#> [1] 2
```
