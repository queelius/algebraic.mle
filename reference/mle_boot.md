# Bootstrapped MLE

Sometimes, the large sample asymptotic theory of MLEs is not applicable.
In such cases, we can use the bootstrap to estimate the sampling
distribution of the MLE.

## Usage

``` r
mle_boot(x)
```

## Arguments

- x:

  the \`boot\` return value

## Value

An `mle_boot` object (wrapper for `boot` object).

## Details

This takes an approach similiar to the \`mle_numerical\` object, which
is a wrapper for a \`stats::optim\` return value, or something that is
compatible with the \`optim\` return value. Here, we take a \`boot\`
object, which is the sampling distribution of an MLE, and wrap it in an
\`mle_boot\` object and then provide a number of methods for the
\`mle_boot\` object that satisfies the concept of an \`mle\` object.

Look up the \`boot\` package for more information on the bootstrap.

## Examples

``` r
# Bootstrap MLE for mean of exponential distribution
set.seed(123)
x <- rexp(50, rate = 2)

# Statistic function: MLE of rate parameter
rate_mle <- function(data, indices) {
  d <- data[indices]
  1 / mean(d)  # MLE of rate is 1/mean
}

# Run bootstrap
boot_result <- boot::boot(data = x, statistic = rate_mle, R = 200)

# Wrap in mle_boot
fit <- mle_boot(boot_result)
params(fit)
#> [1] 1.769331
bias(fit)
#> [1] 0.04782478
confint(fit)
#>            2.5%    97.5%
#> param1 1.196127 2.246885
```
