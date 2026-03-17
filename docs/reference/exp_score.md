# Score Function Generator for the Exponential Distribution

This function returns a function of the rate parameter that computes the
score (derivative of log-likelihood) of the given data assuming an
Exponential distribution.

## Usage

``` r
exp_score(x)
```

## Arguments

- x:

  a numeric vector representing a sample of observations.

## Value

Returns a function that computes the score of \`x\` given a rate
parameter.
