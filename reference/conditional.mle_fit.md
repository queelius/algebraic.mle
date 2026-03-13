# Conditional distribution from an MLE.

Computes the conditional distribution of a subset of parameters given
observed values for other parameters. Uses the closed-form Schur
complement for the multivariate normal. Returns a `dist` object (not an
`mle_fit`), since conditioning loses the MLE context.

## Usage

``` r
# S3 method for class 'mle_fit'
conditional(x, P = NULL, ..., given_indices = NULL, given_values = NULL)
```

## Arguments

- x:

  An `mle_fit` object with at least 2 parameters.

- P:

  Optional predicate function for Monte Carlo fallback.

- ...:

  Additional arguments forwarded to `P`.

- given_indices:

  Integer vector of conditioned parameter indices.

- given_values:

  Numeric vector of observed values.

## Value

A `normal`, `mvn`, or `empirical_dist` object.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
conditional(fit, given_indices = 2, given_values = 4)
#> Normal distribution (mu = 5, var = 0.04) 
```
