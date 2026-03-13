# Method for sampling from an \`mle_fit\` object.

It creates a sampler for the \`mle_fit\` object. It returns a function
that accepts a single parameter \`n\` denoting the number of samples to
draw from the \`mle_fit\` object.

## Usage

``` r
# S3 method for class 'mle_fit'
sampler(x, ...)
```

## Arguments

- x:

  the \`mle_fit\` object to create sampler for

- ...:

  additional arguments to pass

## Value

A function that takes parameter `n` and returns `n` samples from the
asymptotic distribution of the MLE.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
s <- sampler(fit)
head(s(10))
#>          [,1]     [,2]
#> [1,] 4.857919 4.145315
#> [2,] 4.950662 3.803400
#> [3,] 4.809676 3.974528
#> [4,] 4.843019 3.056470
#> [5,] 4.923955 4.519863
#> [6,] 4.884931 4.343917
```
