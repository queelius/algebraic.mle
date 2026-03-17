# A function for computing the sampling distribution of a statistic of the MLE's sampling distribution using the Bootstrap method.

A function for computing the sampling distribution of a statistic of the
MLE's sampling distribution using the Bootstrap method.

## Usage

``` r
mle_boot_loglike(
  mle,
  loglike.gen,
  data = NULL,
  R = NULL,
  method = mle_newton_raphson,
  ...
)
```

## Arguments

- mle:

  a fitted `mle` object.

- loglike.gen:

  a generator for the log-likelihood function; it accepts observations
  and constructs the log-likelihood function

- data:

  data for generating MLEs for the bootstrap resampling.

- R:

  bootstrap replicates

- method:

  method for solving the MLE, defaults to numerically solving the root
  of the gradient of the log-likelihood using Newton-raphson.

- ...:

  additional arguments to pass.
