# Log-Likelihood Function Generator for the Weibull Distribution with a shape and scale parameters that are a function of predictors in \`df\`,

\`resp(df) ~ weibull(shape(df, beta), scale(df, beta))\`

## Usage

``` r
weibull_conditional_shape_scale_loglike(df, resp, shape, scale)
```

## Arguments

- df:

  a numeric vector representing a sample of observations.

- resp:

  a function that returns the response variable given a data frame.

- shape:

  a function that returns the shape given a data frame and a parameter
  vector \`beta\`

- scale:

  a function that returns the scale given a data frame and a parameter
  vector \`beta\`

## Value

Returns a function that computes the conditional log-likelihood

## Details

where \`beta\` is a vector of parameters.

The function returned by this function is suitable for use with
\`optim\` or \`nlm\` to find the maximum likelihood estimate of
\`beta\`.

Note that the expected value given the data is \`scale(...) \* gamma(1 +
1/shape(...))\`.
