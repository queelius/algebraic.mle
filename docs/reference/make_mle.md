# `make_mle` makes an `mle` object.

`make_mle` makes an `mle` object.

## Usage

``` r
make_mle(
  theta.hat,
  loglike = NULL,
  score = NULL,
  sigma = NULL,
  info = NULL,
  obs = NULL,
  sample_size = NULL,
  superclasses = NULL
)
```

## Arguments

- theta.hat:

  the MLE

- loglike:

  the log-likelihood of `theta.hat` given the data

- score:

  the score function evaluated at `theta.hat`

- sigma:

  the variance-covariance matrix of `theta.hat` given that data

- info:

  the information matrix of `theta.hat` given the data

- obs:

  observation (sample) data

- sample_size:

  number of observations in `obs`

- superclasses:

  class hierarchy, with `mle` as base
