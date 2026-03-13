# Method for obtaining the parameters of an \`mle_fit\` object.

Method for obtaining the parameters of an \`mle_fit\` object.

## Usage

``` r
# S3 method for class 'mle_fit'
params(x)
```

## Arguments

- x:

  the \`mle_fit\` object to obtain the parameters of

## Value

Numeric vector of parameter estimates.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
params(fit)
#>     mu sigma2 
#>      5      4 
```
