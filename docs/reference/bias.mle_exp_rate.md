# Computes the bias of an \`exp_mle\` object (exponential mle).

An unbiased estimator of the rate parameter of the exponential
distribution is given by: \`1/(nobs(x)-1)\*bias(x)\`, where \`x\` is an
\`mle_exp_rate\` object.

## Usage

``` r
# S3 method for mle_exp_rate
bias(x, par = NULL, ...)
```

## Arguments

- x:

  the \`mle_exp_rate\` object to compute the bias of.

- par:

  the true rate parameter value. Usually, rate is not known, and so we
  estimate the bias

- ...:

  pass additional arguments
