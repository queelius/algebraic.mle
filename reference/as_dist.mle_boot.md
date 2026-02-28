# Convert a bootstrap MLE to an empirical distribution.

Converts an `mle_boot` object to an
[`empirical_dist`](https://queelius.github.io/algebraic.dist/reference/empirical_dist.html)
built from the bootstrap replicates.

## Usage

``` r
# S3 method for class 'mle_boot'
as_dist(x, ...)
```

## Arguments

- x:

  An `mle_boot` object.

- ...:

  Additional arguments (not used).

## Value

An `empirical_dist` object.

## Examples

``` r
set.seed(123)
x <- rexp(50, rate = 2)
rate_mle <- function(data, indices) 1 / mean(data[indices])
boot_result <- boot::boot(data = x, statistic = rate_mle, R = 200)
fit <- mle_boot(boot_result)
d <- as_dist(fit)
mean(d)
#> [1] 1.817156
```
