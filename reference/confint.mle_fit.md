# Function to compute the confidence intervals of \`mle_fit\` objects.

Function to compute the confidence intervals of \`mle_fit\` objects.

## Usage

``` r
# S3 method for class 'mle_fit'
confint(object, parm = NULL, level = 0.95, use_t_dist = FALSE, ...)
```

## Arguments

- object:

  the \`mle_fit\` object to compute the confidence intervals for

- parm:

  the parameters to compute the confidence intervals for (not used)

- level:

  confidence level, defaults to 0.95 (alpha=.05)

- use_t_dist:

  logical, whether to use the t-distribution to compute the confidence
  intervals.

- ...:

  additional arguments to pass

## Value

Matrix of confidence intervals with columns for lower and upper bounds.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), loglike = -120, nobs = 100L)
confint(fit)
#>            2.5%    97.5%
#> mu     4.608007 5.391993
#> sigma2 2.891277 5.108723
```
