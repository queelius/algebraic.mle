# Constructor for making \`mle\` objects, which provides a common interface for maximum likelihood estimators.

This MLE makes the asymptotic assumption by default. Other MLEs, like
\`mle_boot\`, may not make this assumption.

## Usage

``` r
mle(
  theta.hat,
  loglike = NULL,
  score = NULL,
  sigma = NULL,
  info = NULL,
  obs = NULL,
  nobs = NULL,
  superclasses = NULL
)
```

## Arguments

- theta.hat:

  the MLE

- loglike:

  the log-likelihood of \`theta.hat\` given the data

- score:

  the score function evaluated at \`theta.hat\`

- sigma:

  the variance-covariance matrix of \`theta.hat\` given that data

- info:

  the information matrix of \`theta.hat\` given the data

- obs:

  observation (sample) data

- nobs:

  number of observations in \`obs\`

- superclasses:

  class (or classes) with \`mle\` as base

## Value

An object of class `mle`.

## Examples

``` r
# MLE for normal distribution (mean and variance)
set.seed(123)
x <- rnorm(100, mean = 5, sd = 2)
n <- length(x)
mu_hat <- mean(x)
var_hat <- mean((x - mu_hat)^2)  # MLE of variance

# Asymptotic variance-covariance of MLE
# For normal: Var(mu_hat) = sigma^2/n, Var(var_hat) = 2*sigma^4/n
sigma_matrix <- diag(c(var_hat/n, 2*var_hat^2/n))

fit <- mle(
  theta.hat = c(mu = mu_hat, var = var_hat),
  sigma = sigma_matrix,
  loglike = sum(dnorm(x, mu_hat, sqrt(var_hat), log = TRUE)),
  nobs = n
)

params(fit)
#>       mu      var 
#> 5.180812 3.299602 
vcov(fit)
#>            [,1]      [,2]
#> [1,] 0.03299602 0.0000000
#> [2,] 0.00000000 0.2177475
confint(fit)
#>         2.5%    97.5%
#> mu  4.824788 5.536835
#> var 2.385016 4.214188
```
