# Method for obtaining the observations used by the \`mle_fit\` object \`x\`.

Method for obtaining the observations used by the \`mle_fit\` object
\`x\`.

## Usage

``` r
# S3 method for class 'mle_fit'
obs(x)
```

## Arguments

- x:

  the \`mle_fit\` object to obtain the number of observations for

## Value

The observation data used to fit the MLE, or NULL if not stored.

## Examples

``` r
fit <- mle(theta.hat = c(mu = 5, sigma2 = 4),
  sigma = diag(c(0.04, 0.32)), obs = rnorm(100), nobs = 100L)
head(obs(fit))
#> [1] -0.2266080 -0.5233127  0.1305302  0.2645772  0.3476103  0.3018569
```
