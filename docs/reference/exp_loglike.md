# Log-Likelihood Function Generator for the Exponential Distribution

This function returns a function of the rate parameter that computes the
log-likelihood of the given data assuming an Exponential distribution.

## Usage

``` r
exp_loglike(x)
```

## Arguments

- x:

  a numeric vector representing a sample of observations.

## Value

Returns a function that computes the log-likelihood of \`x\` given a
rate parameter.

## Examples

``` r
if (FALSE) {
  obs <- rexp(100, rate = 0.5)
  loglike <- exp_loglike(obs)
  plot(loglike, from = 0, to = 2)
}
```
