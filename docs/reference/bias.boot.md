# Computes the estimate of the bias of a `boot` object.

Computes the estimate of the bias of a `boot` object.

## Usage

``` r
# S3 method for boot
bias(x, par = NULL, ...)
```

## Arguments

- x:

  the `boot` object to compute the bias of.

- par:

  if the true parameter value is known, you may provide it; otherwise we
  use the MLE of `par`.

- ...:

  pass additional arguments
