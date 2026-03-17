# Computes the bias of an \`mle_normal\` object.

This function computes the bias of an \`mle_normal\` object. It is a
method for the generic function \`bias\`.

## Usage

``` r
# S3 method for mle_normal
bias(x, par = NULL, ...)
```

## Arguments

- x:

  An \`mle_normal\` object (subclass of \`mle\`) to compute the bias of.

- par:

  The true parameter value. If this is unknown (NULL), the bias is
  estimated.

- ...:

  Additional arguments (currently unused).

## Value

A numeric vector of length 2, the bias of the \`mle_normal\` estimator.

## See also

[`bias`](https://queelius.github.io/algebraic.mle/reference/bias.md) for
the generic function.
