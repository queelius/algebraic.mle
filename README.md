
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package: `algebraic.mle`

<!-- badges: start -->

<!-- badges: end -->

An algebra over maximum likelihood estimators (MLE).

MLEs have many desirable, well-defined statistical properties. We define
an algebra over MLEs.

## Installation

You can install `algebraic.mle` from
[GitHub](https://github.com/queelius/algebraic.mle) with:

``` r
install.packages("devtools")
devtools::install_github("queelius/algebraic.mle")
```

## Philosophy

The likelihood (and the log-likelihood) is the most important statistical
object for a parametric model. It's important for the Bayesians and the
frequentists, and its certainly imporatnt for the likelihood-ists (?).
When we assume iid for a sample, the MLE enjoys a host of benefits, from
asymptotic normality (if we look at the log-likelihood, it's just a sum
of log densities, `loglike(theta|data) = sum(log(pmodel(data, theta))`)
to, under the right conditions, being the uniform minimum variance unbaised
estimator of `theta` given the information in `data`.
Since it is asymptotically normal, we can also do a number of algebraic
operations that produce other normal distributions.

## API

The object representing a fitted model is a type of `mle` object, the
maximum likelihood estimator of the model with respect to observed data.

The API mostly consists of generic methods with implementations for
various `mle` type objects. For a full list of functions, see the
[function
reference](https://queelius.github.io/algebraic.mle/reference/index.html)
for `algebraic.mle`.

Let’s fit a conditional exponential model to some data. In this model,
`Y | x ~ EXP(rate(x))` where `rate(x) = exp(b0 + b1*x)`. First, let’s
observe some data from the DGP (data generating process):

``` r
n <- 200
b0 <- -5
b1 <- .5
df <- data.frame(x = rep(NA, n), y = rep(NA, n))
for (i in 1:n) {
    x <- runif(1, -10, 10)
    y <- rexp(n, rate = exp(b0 + b1 * x))
    df[i, ] <- c(x, y)
}
```

Now, we define three functions, `resp`, `rate`, and `loglik`, which
defines our conditional model.

``` r
resp <- function(df) df$y
rate <- function(df, beta) exp(beta[1] + beta[2] * df$x)
loglike <- function(df, resp, rate) {
  function(beta) sum(dexp(x = resp(df), rate = rate(df, beta), log = TRUE))
}
```

Let’s fit the model to the data in `df`. We’ll use the `optim` function in
`stats` to fit the model and then wrap it into an `mle` object using
`mle_numerical`.

``` r
library(algebraic.mle)
sol <- mle_numerical(optim(par = c(0, 0),
    fn = loglike(df, resp, rate),
    control = list(fnscale = -1),
    hessian = TRUE))
```

We have fit the model. We can do a lot of things with this model, for instance,
we can print out summary info:
``` r
summary(sol)
#> Maximum likelihood estimator of type mle_numerical is normally distributed.
#> The estimates of the parameters are given by:
#> [1] -5.0366809  0.4883971
#> The standard error is  0.07082793 0.01182874 .
#> The asymptotic 95% confidence interval of the parameters are given by:
#>              2.5%      97.5%
#> param1 -5.1531825 -4.9201794
#> param2  0.4689405  0.5078536
#> The MSE of the estimator is  0.005156515 .
#> The log-likelihood is  -1171.433 .
#> The AIC is  2346.866 .
```

See the API reference for more options.

Let’s plot the data, the true expectation of the
DGP (green), and the estimate of the DGP (red):

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

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

You can see tutorials for more examples of using the package in the
[vignettes](https://queelius.github.io/algebraic.mle/articles/index.html).
These vignettes are very simple illustrations, in many ways simpler than
what we showed above, but they reveal more of the API and how to use it.
