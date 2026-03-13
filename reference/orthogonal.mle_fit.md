# Method for determining the orthogonal components of an \`mle_fit\` object \`x\`.

Method for determining the orthogonal components of an \`mle_fit\`
object \`x\`.

## Usage

``` r
# S3 method for class 'mle_fit'
orthogonal(x, tol = sqrt(.Machine$double.eps), ...)
```

## Arguments

- x:

  the \`mle_fit\` object

- tol:

  the tolerance for determining if a number is close enough to zero

- ...:

  pass additional arguments

## Value

Logical matrix indicating which off-diagonal FIM elements are
approximately zero (orthogonal parameters), or NULL if FIM unavailable.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), info = solve(diag(c(0.04, 0.32))),
  nobs = 100L)
orthogonal(fit)
#>       [,1]  [,2]
#> [1,] FALSE  TRUE
#> [2,]  TRUE FALSE
```
