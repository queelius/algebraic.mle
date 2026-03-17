# Maximum Likelihood Estimation (MLE) for the Exponential Distribution

This function computes the Maximum Likelihood Estimation (MLE) of the
rate parameter for a given sample that is assumed to be independently
and identically distributed (i.i.d) and drawn from the Exponential
distribution.

## Usage

``` r
mle_exp(x, keep_obs = FALSE)
```

## Arguments

- x:

  a numeric vector representing a sample of observations.

- keep_obs:

  logical. If TRUE, the observations are stored with the \`mle\` object.
  Default is FALSE.

## Value

Returns an \`mle\` object that contains:

## Examples

``` r
if (FALSE) {
  obs <- rexp(100, rate = 0.5)
  rate.hat <- mle_exp(obs)
  summary(rate.hat)
}
```
