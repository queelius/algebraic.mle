# Log-Likelihood Function Generator for the Exponential Distribution with a rate parameter that is a function of predictors in \`df\`,

\`resp(df) ~ exp(rate(df, beta))\`

## Usage

``` r
exp_conditional_loglike(df, resp, rate)
```

## Arguments

- df:

  a numeric vector representing a sample of observations.

- resp:

  a function that returns the response variable given a row from a data
  frame.

- rate:

  a function that returns the rate given a row from a data frame and a
  parameter vector \`beta\`

## Value

Returns a function that computes the conditional log-likelihood

## Details

where \`beta\` is a vector of parameters.

The function returned by this function is suitable for use with
\`optim\` or \`nlm\` to find the maximum likelihood estimate of
\`beta\`.

Note that the expected value is \`1/rate(...)\`.
