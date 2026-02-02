# Computes the estimate of the bias of a \`mle_boot\` object.

Computes the estimate of the bias of a \`mle_boot\` object.

## Usage

``` r
# S3 method for class 'mle_boot'
bias(x, theta = NULL, ...)
```

## Arguments

- x:

  the \`mle_boot\` object to compute the bias of.

- theta:

  true parameter value (not used for \`mle_boot\`).

- ...:

  pass additional arguments (not used)

## Value

Numeric vector of estimated bias (mean of bootstrap replicates minus
original estimate).
