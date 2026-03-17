# MLE of the rate parameter when we assume the sample is i.i.d. and drawn from the exponential distribution.

Of course, the draws are unlikely to be exponential, but the exponential
may be a sufficiently good model. Hypothesis testing, such as relative
likelihoods, can be used to assess the appropriateness of the
exponential model to the data.

## Usage

``` r
mle_exp_rate(x, keep_obs = F)
```

## Arguments

- x:

  a sample of observations

- keep_obs:

  store the observations with the \`mle\` object, default is \`F\`

## Value

an \`mle\` object.
