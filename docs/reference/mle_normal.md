# MLE of the (mu,var) parameter vector when we assume the sample is i.i.d. and drawn from the normal distribution.

Of course, the draws are unlikely to be normal, but the normal
distribution is often a good model. Hypothesis testing, such as relative
likelihoods, can be used to assess the appropriateness of the normal
model to the data.

## Usage

``` r
mle_normal(x, keep_obs = TRUE)
```

## Arguments

- x:

  a sample of observations

- keep_obs:

  whether to store observations in \`mle_normal\` object

## Value

an \`mle\` object.
