# Convert an MLE to its asymptotic distribution.

Internal helper that creates a
[`normal`](https://queelius.github.io/algebraic.dist/reference/normal.html)
(univariate) or
[`mvn`](https://queelius.github.io/algebraic.dist/reference/mvn.html)
(multivariate) distribution from an `mle_fit` object using the
asymptotic normality of the MLE.

## Usage

``` r
mle_to_dist(x)
```

## Arguments

- x:

  An `mle_fit` object with a non-NULL variance-covariance matrix.

## Value

A `normal` or `mvn` distribution object.
