# Computes the estimate of the MSE of a `boot` object.

Computes the estimate of the MSE of a `boot` object.

## Usage

``` r
# S3 method for boot
mse(x, par = NULL, ...)
```

## Arguments

- x:

  the `boot` object to compute the MSE of.

- par:

  if the true parameter value is known, you may provide it; otherwise we
  use the MLE of `par`.

- ...:

  pass additional arguments
