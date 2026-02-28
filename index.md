# algebraic.mle

The maximum likelihood estimator is a **technology**: under regularity
conditions, any MLE is asymptotically normal with variance given by the
inverse Fisher information. `algebraic.mle` exploits that structure by
defining an **algebra** over MLEs — you can compose them, combine them,
transform them, and convert them to distribution objects, all while
propagating uncertainty automatically.

## Installation

Install from CRAN:

``` r
install.packages("algebraic.mle")
```

Or the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("queelius/algebraic.mle")
```

``` r
library(algebraic.mle)
library(algebraic.dist)
```

## The algebra

### `joint()` — compose independent experiments

Two labs estimate different parameters from independent experiments.
[`joint()`](https://queelius.github.io/algebraic.mle/reference/joint.md)
concatenates the parameter vectors and produces a block-diagonal
covariance:

``` r
fit_rate <- mle(theta.hat = c(lambda = 2.1), sigma = matrix(0.04), nobs = 50L)
fit_shape <- mle(
  theta.hat = c(k = 1.5, s = 3.2),
  sigma = matrix(c(0.10, 0.02,
                    0.02, 0.30), 2, 2),
  nobs = 100L
)

j <- joint(fit_rate, fit_shape)
params(j)
#> lambda      k      s 
#>    2.1    1.5    3.2
vcov(j)
#>      [,1] [,2] [,3]
#> [1,] 0.04 0.00 0.00
#> [2,] 0.00 0.10 0.02
#> [3,] 0.00 0.02 0.30
```

The off-diagonal blocks are zero because the experiments are
independent.

### `combine()` — pool repeated estimates

Three labs each estimate the same rate parameter.
[`combine()`](https://queelius.github.io/algebraic.mle/reference/combine.md)
uses inverse-variance weighting to produce the minimum-variance unbiased
combination:

``` r
fit1 <- mle(theta.hat = c(lambda = 2.1), sigma = matrix(0.04), nobs = 50L)
fit2 <- mle(theta.hat = c(lambda = 1.9), sigma = matrix(0.02), nobs = 100L)
fit3 <- mle(theta.hat = c(lambda = 2.0), sigma = matrix(0.03), nobs = 70L)

comb <- combine(fit1, fit2, fit3)
params(comb)
#>   lambda 
#> 1.976923
se(comb)             # smaller than any individual SE
#> [1] 0.09607689
c(se(fit1), se(fit2), se(fit3))
#> [1] 0.2000000 0.1414214 0.1732051
```

### `rmap()` — transform via the delta method

The MLE of any smooth function of the parameters is the function applied
to the MLE (invariance property).
[`rmap()`](https://queelius.github.io/algebraic.dist/reference/rmap.html)
propagates uncertainty through the Jacobian:

``` r
# MLE for an exponential rate
fit_exp <- mle(theta.hat = c(lambda = 2.0), sigma = matrix(0.08), nobs = 50L)

# Transform rate to mean lifetime: E[T] = 1/lambda
mean_life <- rmap(fit_exp, function(theta) c(mean_life = 1 / theta[1]),
                  method = "delta")
params(mean_life)
#> mean_life 
#>       0.5
se(mean_life)
#> [1] 0.07071068
confint(mean_life)
#>                2.5%     97.5%
#> mean_life 0.3614096 0.6385904
```

### `as_dist()` — bridge to distribution algebra

Convert an MLE to its asymptotic normal distribution, then use the full
distribution algebra from `algebraic.dist`:

``` r
d1 <- as_dist(fit1)   # Normal(2.1, 0.04)
d2 <- as_dist(fit2)   # Normal(1.9, 0.02)

# Sum of independent normals is normal
d_sum <- d1 + d2
d_sum
#> Normal distribution (mu = 4, var = 0.06)
```

## Inference

All MLE objects support a common interface for statistical inference:

``` r
fit <- mle(
  theta.hat = c(mu = 5.1, sigma2 = 3.8),
  sigma = diag(c(0.04, 0.29)),
  loglike = -154.3,
  nobs = 100L
)

confint(fit)
#>            2.5%    97.5%
#> mu     4.708007 5.491993
#> sigma2 2.744527 4.855473
se(fit)
#> [1] 0.2000000 0.5385165
aic(fit)
#> [1] 312.6
```

Draw samples from the asymptotic distribution:

``` r
samp <- sampler(fit)
head(samp(5))
#>          [,1]     [,2]
#> [1,] 4.881999 3.464467
#> [2,] 5.199934 3.735020
#> [3,] 5.067908 4.302101
#> [4,] 4.964794 3.924047
#> [5,] 5.059629 2.571551
```

## Creating MLE objects

### `mle()` — direct construction

When you know the point estimate and variance-covariance (e.g., from
analytical results), construct directly:

``` r
fit <- mle(
  theta.hat = c(mu = 5.0),
  sigma = matrix(0.04),
  loglike = -120.5,
  nobs = 100L
)
summary(fit)
#> Maximum likelihood estimator of type mle is normally distributed.
#> The estimates of the parameters are given by:
#> mu 
#>  5 
#> The standard error is  0.2 .
#> The asymptotic 95% confidence interval of the parameters are given by:
#>        2.5%    97.5%
#> mu 4.608007 5.391993
#> The MSE of the estimator is  0.04 .
#> The log-likelihood is  -120.5 .
#> The AIC is  243 .
```

### `mle_numerical()` — wrap `optim()` output

Fit a conditional exponential model `Y | x ~ Exp(rate(x))` where
`rate(x) = exp(b0 + b1 * x)`:

``` r
# Simulate data from the model
set.seed(1231)
n <- 75
b0_true <- -0.1
b1_true <-  0.5
df <- data.frame(
  x = runif(n, -10, 10),
  y = NA
)
df$y <- rexp(n, rate = exp(b0_true + b1_true * df$x))

# Define log-likelihood
loglik <- function(beta) {
  rate <- exp(beta[1] + beta[2] * df$x)
  sum(dexp(df$y, rate = rate, log = TRUE))
}

# Fit and wrap
par0 <- c(b0 = 0, b1 = 0)
sol <- mle_numerical(optim(
  par = par0,
  fn = loglik,
  control = list(fnscale = -1),
  hessian = TRUE
))
summary(sol)
#> Maximum likelihood estimator of type mle_numerical is normally distributed.
#> The estimates of the parameters are given by:
#>         b0         b1 
#> -0.1425359  0.5226732 
#> The standard error is  0.1154688 0.01986773 .
#> The asymptotic 95% confidence interval of the parameters are given by:
#>          2.5%      97.5%
#> b0 -0.3688506 0.08377888
#> b1  0.4837332 0.56161328
#> The MSE of the individual components in a multivariate estimator is:
#>              [,1]         [,2]
#> [1,] 1.333305e-02 4.876544e-05
#> [2,] 4.876544e-05 3.947267e-04
#> The log-likelihood is  -90.58435 .
#> The AIC is  185.1687 .
confint(sol)
#>          2.5%      97.5%
#> b0 -0.3688506 0.08377888
#> b1  0.4837332 0.56161328
```

### `mle_boot()` — bootstrap MLE

When sample size is small or regularity conditions are in doubt, use the
bootstrap:

``` r
set.seed(42)
x <- rexp(30, rate = 2)
boot_fn <- function(data, idx) 1 / mean(data[idx])
b <- boot::boot(data = x, statistic = boot_fn, R = 999)
fit_boot <- mle_boot(b)
params(fit_boot)
#> [1] 2.019028
confint(fit_boot, type = "perc")
#>            2.5%    97.5%
#> param1 1.367237 3.175378
```

## Links

- [Package documentation](https://queelius.github.io/algebraic.mle/)
- [Vignettes](https://queelius.github.io/algebraic.mle/articles/index.html)
- [Function
  reference](https://queelius.github.io/algebraic.mle/reference/index.html)
- [GitHub repository](https://github.com/queelius/algebraic.mle)
