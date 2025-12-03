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
