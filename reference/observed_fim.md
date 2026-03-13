# Generic method for computing the observed FIM of an \`mle_fit\` object.

Fisher information is a way of measuring the amount of information that
an observable random variable \`X\` carries about an unknown parameter
\`theta\` upon which the probability of \`X\` depends.

## Usage

``` r
observed_fim(x, ...)
```

## Arguments

- x:

  the object to obtain the fisher information of

- ...:

  additional arguments to pass

## Value

The observed Fisher Information Matrix.

## Details

The inverse of the Fisher information matrix is the variance-covariance
of the MLE for \`theta\`.

Some MLE objects do not have an observed FIM, e.g., if the MLE's
sampling distribution was bootstrapped.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), info = solve(diag(c(0.04, 0.32))),
  loglike = -120, nobs = 100L)
observed_fim(fit)
#>      [,1]  [,2]
#> [1,]   25 0.000
#> [2,]    0 3.125
```
