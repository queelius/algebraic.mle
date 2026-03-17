# A log-likelihood function generator given data \`x\` for the weibull distribution.

The returned log-likelihood function takes a single vector \`par\` of
size \`2\` (at least) where the first component is the shape parameter
\`k\` and the second component is the scale parameter \`lambda\` such
that E\[X\] = lambda \* GAMMA(1 + 1/k).

## Usage

``` r
weibull_loglike(x)
```

## Arguments

- x:

  data

## Value

log-likelihood function for the weibull distribution with shape
parameter \`k\` and scale parameter \`lambda\` given data \`x\`.
(Accepts a parameter vector \`par\` of size 2.)
