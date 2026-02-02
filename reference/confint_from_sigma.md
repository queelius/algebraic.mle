# Function to compute the confidence intervals from a variance-covariance matrix

Function to compute the confidence intervals from a variance-covariance
matrix

## Usage

``` r
confint_from_sigma(sigma, theta, level = 0.95)
```

## Arguments

- sigma:

  either the variance-covariance matrix or the vector of variances of
  the parameter estimator

- theta:

  the point estimate

- level:

  confidence level, defaults to 0.95 (alpha=.05)

## Value

Matrix of confidence intervals with rows for each parameter and columns
for lower and upper bounds.

## Examples

``` r
# Compute CI for a bivariate parameter
theta <- c(mu = 5.2, sigma2 = 4.1)
vcov_matrix <- diag(c(0.1, 0.5))  # Variance of estimators

confint_from_sigma(vcov_matrix, theta)
#>            2.5%    97.5%
#> mu     4.580205 5.819795
#> sigma2 2.714096 5.485904
confint_from_sigma(vcov_matrix, theta, level = 0.99)
#>            0.5%    99.5%
#> mu     4.385451 6.014549
#> sigma2 2.278614 5.921386
```
