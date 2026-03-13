# Function for obtaining the observed FIM of an \`mle_fit\` object.

Function for obtaining the observed FIM of an \`mle_fit\` object.

## Usage

``` r
# S3 method for class 'mle_fit'
observed_fim(x, ...)
```

## Arguments

- x:

  the \`mle_fit\` object to obtain the FIM of.

- ...:

  pass additional arguments

## Value

The observed Fisher Information Matrix, or NULL if not available.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), info = solve(diag(c(0.04, 0.32))),
  nobs = 100L)
observed_fim(fit)
#>      [,1]  [,2]
#> [1,]   25 0.000
#> [2,]    0 3.125
```
