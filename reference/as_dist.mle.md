# Convert an MLE to a distribution object.

Converts an `mle` object to its asymptotic
[`normal`](https://queelius.github.io/algebraic.dist/reference/normal.html)
or [`mvn`](https://queelius.github.io/algebraic.dist/reference/mvn.html)
distribution. The MLE must have a variance-covariance matrix available.

## Usage

``` r
# S3 method for class 'mle'
as_dist(x, ...)
```

## Arguments

- x:

  An `mle` object.

- ...:

  Additional arguments (not used).

## Value

A `normal` (univariate) or `mvn` (multivariate) distribution.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5), sigma = matrix(0.25))
d <- as_dist(fit)
mean(d)   # 5
#> [1] 5
vcov(d)   # 0.25
#> [1] 0.25

fit2 <- mle(theta.hat = c(a = 1, b = 2), sigma = diag(c(0.1, 0.2)))
d2 <- as_dist(fit2)
mean(d2)  # c(1, 2)
#> [1] 1 2
```
