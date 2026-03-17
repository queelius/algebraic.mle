# log-likelihood function generator given data \`x\` for the weibull distribution.

The returned log-likelihood function takes a single vector \`theta\` of
size \`2\` (at least) where the first component is the shape parameter
\`k\` and the second component is the scale parameter \`lambda\`.

## Usage

``` r
weibull_shape_scale_loglike(x)
```

## Arguments

- x:

  data

## Details

It can be used in statistical models or optimization algorithms to
estimate the parameters of the Weibull distribution.
