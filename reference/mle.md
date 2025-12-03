# Constructor for making \`mle\` objects, which provides a common interface for maximum likelihood estimators.

This MLE makes the asymptotic assumption by default. Other MLEs, like
\`mle_boot\`, may not make this assumption.

## Usage

``` r
mle(
  theta.hat,
  loglike = NULL,
  score = NULL,
  sigma = NULL,
  info = NULL,
  obs = NULL,
  nobs = NULL,
  superclasses = NULL
)
```

## Arguments

- theta.hat:

  the MLE

- loglike:

  the log-likelihood of \`theta.hat\` given the data

- score:

  the score function evaluated at \`theta.hat\`

- sigma:

  the variance-covariance matrix of \`theta.hat\` given that data

- info:

  the information matrix of \`theta.hat\` given the data

- obs:

  observation (sample) data

- nobs:

  number of observations in \`obs\`

- superclasses:

  class (or classes) with \`mle\` as base
