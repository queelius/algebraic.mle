# Computes the bias of an \`mle\` object assuming the large sample approximation is valid and the MLE regularity conditions are satisfied. In this case, the bias is zero (or zero vector).

This is not a good estimate of the bias in general, but it's arguably
better than returning \`NULL\`.

## Usage

``` r
# S3 method for class 'mle'
bias(x, theta = NULL, ...)
```

## Arguments

- x:

  the \`mle\` object to compute the bias of.

- theta:

  true parameter value. normally, unknown (NULL), in which case we
  estimate the bias (say, using bootstrap)

- ...:

  additional arguments to pass

## Value

Numeric vector of zeros (asymptotic bias is zero under regularity
conditions).
