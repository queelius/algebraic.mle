# Generic method for computing the bias of an estimator object.

Generic method for computing the bias of an estimator object.

## Usage

``` r
bias(x, theta, ...)
```

## Arguments

- x:

  the object to compute the bias of.

- theta:

  true parameter value. usually, this is unknown (NULL), in which case
  we estimate the bias

- ...:

  pass additional arguments

## Value

The bias of the estimator. The return type depends on the specific
method.
