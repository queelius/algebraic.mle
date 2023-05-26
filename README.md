
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package: `algebraic.mle`

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/queelius/algebraic.mle/blob/master/LICENSE)
<!-- badges: end -->

`algebraic.mle` is an R package that provides an algebra over Maximum
Likelihood Estimators (MLEs). These estimators possess many desirable,
well-defined statistical properties which the package helps you
manipulate and utilize.

## Installation

`algebraic.mle` can be installed from GitHub by using the devtools
package in R:

``` r
install.packages("devtools")
devtools::install_github("queelius/algebraic.mle")
```

## Purpose

The primary focus of this package is the likelihood and the
log-likelihood, fundamental statistical concepts for a parametric model.
They are integral to both Bayesian and frequentist statistics, as well
as for those who prioritize likelihood. The algebraic.mle package
enables easy handling of MLEs, which, under certain conditions and
assumptions (such as independence and identical distribution (iid) for a
sample), present numerous advantages, including asymptotic normality and
being the uniformly minimum variance unbiased estimator of theta.

## API Overview

The main object in the `algebraic.mle` package is the `mle` object,
which represents a fitted model. The package provides a number of
generic methods designed for `mle` objects. A comprehensive list of
functions is available in the [function
reference](https://queelius.github.io/algebraic.mle/reference/index.html)
for `algebraic.mle`.

## Example

Here is an example of fitting a conditional exponential model to some
data using `algebraic.mle`. The true DGP is given by `Y | x ~
EXP(rate(x))` where `rate(x) = exp(b0 + b1*x)`, and we do not care how
`x` is distributed.

Let’s fit a conditional exponential model to some data. In this model,
`Y | x ~ EXP(rate(x))` where `rate(x) = exp(b0 + b1*x)`. First, let’s
define the DGP (data generating process):

``` r
b0 <- -.1
b1 <- 0.5

dgp <- function(n, x) {
    rate <- exp(b0 + b1 * x)
    rexp(n, rate) + rnorm(n, 0, 1e-3)
}
```

Let’s generate some date:

``` r
n <- 75
set.seed(1231)
df <- data.frame(x = rep(NA, n), y = rep(NA, n))
for (i in 1:n) {
    x <- runif(1, -10, 10)
    y <- dgp(n = 1, x = x)
    df[i, ] <- c(x, y)
}
```

Now, we define two functions, `resp`, `rate`, and `loglik` function
which will be used to define the model.

``` r
resp <- function(df) df$y
rate <- function(df, beta) exp(beta[1] + beta[2] * df$x)
loglik <- function(df, resp, rate) {
  function(beta) sum(dexp(x = resp(df), rate = rate(df, beta), log = TRUE))
}
```

Let’s fit the model. We’ll use the `optim` function in `stats` to fit
the model and then wrap it into an `mle` object using `mle_numerical`.

``` r
library(algebraic.mle)

# initial guess for the parameters
par0 <- c(0, 0)
names(par0) <- c("b0", "b1")

sol <- mle_numerical(optim(par = par0,
    fn = loglik(df, resp, rate),
    control = list(fnscale = -1),
    hessian = TRUE))
summary(sol)
#> Maximum likelihood estimator of type mle_numerical is normally distributed.
#> The estimates of the parameters are given by:
#>         b0         b1 
#> -0.2253626  0.4560893 
#> The standard error is  0.1167634 0.02145606 .
#> The asymptotic 95% confidence interval of the parameters are given by:
#>          2.5%       97.5%
#> b0 -0.4174213 -0.03330395
#> b1  0.4207972  0.49138139
#> The MSE of the estimator is  0.01409405 .
#> The log-likelihood is  -119.6977 .
#> The AIC is  243.3954 .
```

Let’s plot it:

``` r
# plot the x-y points from the data frame
plot(df$x,df$y)

# now overlay a plot of the conditional mean
x <- seq(-10, 10, .1)
b0.hat <- point(sol)[1]
b1.hat <- point(sol)[2]
y.hat <- 1/exp(b0.hat + b1.hat*x)
y <- 1/exp(b0 + b1*x)
lines(x, y, col = "green", lwd = 10)
lines(x, y.hat, col = "blue", lwd = 10)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

### Hypothesis test and model selection

Let’s test the hypothesis that `b0 = 0` using a likelihood ratio test.
We can use the LRT because this null model is a special case (nested) of
the full model. The null model is `Y | x ~ EXP(rate(x))` where `rate(x)
= exp(b1*x)`, while the full model is `Y | x ~ EXP(rate(x))` where
`rate(x) = exp(b0 + b1*x)`.

``` r
# construct null model where b1 = 0
rateb0_zero <- function(df, b1) exp(b1 * df$x)

# initial guess for the parameters
b1 <- 0
names(b1) <- "b1"

# fit the model under the null hypothesis
sol_null_b0 <- mle_numerical(optim(par = b1,
    fn = loglik(df, resp, rateb0_zero),
    control = list(fnscale = -1),
    hessian = TRUE,
    method = "BFGS"))
summary(sol_null_b0)
#> Maximum likelihood estimator of type mle_numerical is normally distributed.
#> The estimates of the parameters are given by:
#>        b1 
#> 0.4617093 
#> The standard error is  0.01899941 .
#> The asymptotic 95% confidence interval of the parameters are given by:
#>         2.5%     97.5%
#> b1 0.4304581 0.4929605
#> The MSE of the estimator is  0.0003609774 .
#> The log-likelihood is  -121.7164 .
#> The AIC is  245.4328 .
```

Let’s compute the likelihood ratio test statistic and p-value:

``` r
(lrt <- -2 * (loglike(sol_null_b0) - loglike(sol)))
#> [1] 4.037435
pchisq(lrt, df = 1, lower.tail = FALSE) # compute the p-value
#> [1] 0.04450142
```

We see that the `p > 0.05`, so the data appears to be *compatible* with
the null hypothesis `b0 = 0`.

If we wanted to do model selection, we could use the AIC:

``` r
aic(sol)
#> [1] 243.3954
aic(sol_null_b0)
#> [1] 245.4328
```

Since the AIC of the null model is greater than the AIC of the full
model, we conclude that the data is “more” compatible with the full
model.

So, even though the data appears to be compatible with the null
hypothesis according to the LRT, the AIC suggests that the full model is
better. We actually know the DGP and both models are reasonable
approximations, but neither of course models the DGP exactly.

“All models are wrong, but some are useful.” - George Box

Eventually, if we have a sufficiently large sample, any model that is
not the DGP can be discarded, but reality is so complex that we will
never have a large enough sample and we will never be able to come up
with a model that is exactly the DGP.

Let’s do another test, ![b1
= 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b1%20%3D%200
"b1 = 0"), i.e., it’s an unconditional exponential model, or just a
standard exponential distribution.

``` r
rate_b1_zero <- function(df, b0) exp(b0)
# fit the model under the null hypothesis
sol_null_b1 <- mle_numerical(optim(par = 0,
    fn = loglik(df, resp, rate_b1_zero),
    control = list(fnscale = -1),
    hessian = TRUE,
    method = "BFGS"))
(lrt.null_b1 <- -2 * (loglike(sol_null_b1) - loglike(sol)))
#> [1] 285.0265
pchisq(lrt.null_b1, df = 1, lower.tail = FALSE) # compute the p-value
#> [1] 6.029289e-64
```

This has a `p`-value of essentially zero, so we reject the null
hypothesis that `b1 = 0`.

You can see tutorials for more examples of using the package in the
[vignettes](https://queelius.github.io/algebraic.mle/articles/index.html).
