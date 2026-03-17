# The Algebra of MLEs

## Introduction

The maximum likelihood estimator (MLE) is a technology. Under regularity
conditions, any MLE is asymptotically normal:
$\widehat{\theta} \sim N\left( \theta,I(\theta)^{- 1} \right)$, where
$I(\theta)$ is the Fisher information matrix. This fact is remarkable
because it holds regardless of the underlying model – exponential,
Weibull, normal, or anything else. The MLE reduces the problem of
inference to a mean vector and a covariance matrix.

`algebraic.mle` is the algebra that follows from this fact. It does not
find MLEs – it takes them as input and provides operations for
composition, transformation, and inference. Given two independent
experiments,
[`joint()`](https://queelius.github.io/algebraic.mle/reference/joint.md)
composes their MLEs into a single joint estimator with block-diagonal
covariance. Given three labs estimating the same quantity,
[`combine()`](https://queelius.github.io/algebraic.mle/reference/combine.md)
produces the optimal inverse-variance weighted estimate. Given a
transformation $g(\theta)$,
[`rmap()`](https://queelius.github.io/algebraic.dist/reference/rmap.html)
propagates uncertainty through the delta method. And
[`as_dist()`](https://queelius.github.io/algebraic.dist/reference/as_dist.html)
bridges to the distribution algebra in `algebraic.dist`, where normal
distributions can be further composed.

This vignette walks through each operation with concrete examples,
culminating in a full pipeline from independent experiments to system
reliability inference.

``` r
library(algebraic.mle)
```

## Creating MLEs

The package provides three constructors for wrapping estimation results
into the MLE algebra. Each produces an object that supports the full
interface:
[`params()`](https://queelius.github.io/algebraic.dist/reference/params.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`se()`](https://queelius.github.io/algebraic.mle/reference/se.md),
[`aic()`](https://queelius.github.io/algebraic.mle/reference/aic.md),
and the algebraic operations.

### Direct construction with `mle()`

When you know the point estimate and its variance-covariance matrix
(e.g., from analytical results or published tables), construct an MLE
directly:

``` r
# Exponential rate from a published study: lambda_hat = 0.5, se = 0.07, n = 200
fit_pub <- mle(
  theta.hat = c(lambda = 0.5),
  sigma = matrix(0.07^2),
  nobs = 200L
)
params(fit_pub)
#> lambda 
#>    0.5
se(fit_pub)
#> [1] 0.07
```

### From numerical optimization with `mle_numerical()`

The most common path: maximize a log-likelihood with
[`optim()`](https://rdrr.io/r/stats/optim.html), then wrap the result.
The Hessian at the optimum gives the variance-covariance matrix.

``` r
set.seed(42)
x <- rexp(80, rate = 2)

loglik <- function(rate) {
  if (rate <= 0) return(-Inf)
  sum(dexp(x, rate = rate, log = TRUE))
}

result <- optim(
  par = c(lambda = 1),
  fn = loglik,
  method = "Brent", lower = 0.01, upper = 20,
  hessian = TRUE,
  control = list(fnscale = -1)
)

fit_num <- mle_numerical(result, options = list(nobs = length(x)))
params(fit_num)
#> [1] 1.687
se(fit_num)
#> [1] 0.1886
```

### From bootstrap with `mle_boot()`

When the sample is small or regularity conditions are questionable,
bootstrap the MLE and wrap the `boot` object:

``` r
set.seed(42)
y <- rexp(30, rate = 2)

boot_result <- boot::boot(
  data = y,
  statistic = function(d, i) 1 / mean(d[i]),
  R = 999
)

fit_boot <- mle_boot(boot_result)
params(fit_boot)
#> [1] 2.019
se(fit_boot)
#> [1] 0.4629
confint(fit_boot, type = "perc")
#>         2.5% 97.5%
#> param1 1.367 3.175
```

## Composing Independent MLEs

When independent experiments measure different parameters of a system,
[`joint()`](https://queelius.github.io/algebraic.mle/reference/joint.md)
composes their MLEs into a single joint estimator. The result has a
block-diagonal variance-covariance matrix because the experiments are
independent – there is no cross-covariance between their parameter
estimates.

**Motivating example**: A series system has two components. Lab A tests
component 1 and estimates its failure rate $\lambda_{1}$. Lab B tests
component 2 and estimates $\lambda_{2}$. Neither lab knows about the
other component, but we need the joint parameter vector
$\left( \lambda_{1},\lambda_{2} \right)$ for system-level inference.

``` r
# Lab A: component 1 failure rate
fit_A <- mle(
  theta.hat = c(lambda1 = 0.02),
  sigma = matrix(0.001^2),
  nobs = 150L
)

# Lab B: component 2 failure rate
fit_B <- mle(
  theta.hat = c(lambda2 = 0.05),
  sigma = matrix(0.003^2),
  nobs = 80L
)

# Joint MLE: block-diagonal covariance
fit_joint <- joint(fit_A, fit_B)
params(fit_joint)
#> lambda1 lambda2 
#>    0.02    0.05
vcov(fit_joint)
#>       [,1]  [,2]
#> [1,] 1e-06 0e+00
#> [2,] 0e+00 9e-06
```

The off-diagonal entries are zero – the hallmark of independence. The
joint MLE supports the full interface:

``` r
confint(fit_joint)
#>            2.5%   97.5%
#> lambda1 0.01804 0.02196
#> lambda2 0.04412 0.05588
se(fit_joint)
#> [1] 0.001 0.003
```

Use
[`marginal()`](https://queelius.github.io/algebraic.dist/reference/marginal.html)
to recover individual component estimates from the joint:

``` r
# Recover component 2 parameters
fit_B_recovered <- marginal(fit_joint, 2)
params(fit_B_recovered)
#> lambda2 
#>    0.05
se(fit_B_recovered)
#> [1] 0.003
```

[`joint()`](https://queelius.github.io/algebraic.mle/reference/joint.md)
extends to any number of independent MLEs:

``` r
fit_C <- mle(
  theta.hat = c(lambda3 = 0.01),
  sigma = matrix(0.0005^2),
  nobs = 200L
)

fit_system <- joint(fit_A, fit_B, fit_C)
params(fit_system)
#> lambda1 lambda2 lambda3 
#>    0.02    0.05    0.01
```

## Combining Repeated Estimates

When multiple independent experiments estimate the **same** parameter,
[`combine()`](https://queelius.github.io/algebraic.mle/reference/combine.md)
produces the optimal estimator via inverse-variance weighting. The
combined estimate weights more precise estimates more heavily, and the
combined variance is always less than any individual variance.

**Motivating example**: Three labs independently estimate the failure
rate $\lambda$ of the same component type.

``` r
lab1 <- mle(theta.hat = c(lambda = 0.050), sigma = matrix(0.004^2), nobs = 50L)
lab2 <- mle(theta.hat = c(lambda = 0.048), sigma = matrix(0.002^2), nobs = 200L)
lab3 <- mle(theta.hat = c(lambda = 0.053), sigma = matrix(0.003^2), nobs = 100L)

combined <- combine(lab1, lab2, lab3)
params(combined)
#>  lambda 
#> 0.04961
se(combined)
#> [1] 0.001536
```

The combined standard error is smaller than any individual standard
error:

``` r
cat("Lab SEs:     ", se(lab1), se(lab2), se(lab3), "\n")
#> Lab SEs:      0.004 0.002 0.003
cat("Combined SE: ", se(combined), "\n")
#> Combined SE:  0.001536
```

The combined estimate is pulled toward the more precise estimates. Lab 2
has the smallest variance, so it receives the most weight.

When all estimates have equal precision,
[`combine()`](https://queelius.github.io/algebraic.mle/reference/combine.md)
reduces to the simple average:

``` r
eq1 <- mle(theta.hat = c(mu = 10.0), sigma = matrix(1))
eq2 <- mle(theta.hat = c(mu = 12.0), sigma = matrix(1))
eq3 <- mle(theta.hat = c(mu = 11.0), sigma = matrix(1))

eq_combined <- combine(eq1, eq2, eq3)
params(eq_combined)  # (10 + 12 + 11) / 3 = 11
#> mu 
#> 11
```

## Transformations via Invariance

By the invariance property of the MLE, if $\widehat{\theta}$ is the MLE
of $\theta$, then $g\left( \widehat{\theta} \right)$ is the MLE of
$g(\theta)$ for any function $g$. The
[`rmap()`](https://queelius.github.io/algebraic.dist/reference/rmap.html)
function computes this transformation and propagates uncertainty using
the delta method:
$$\text{Var}\left( g\left( \widehat{\theta} \right) \right) \approx J\,\text{Var}\left( \widehat{\theta} \right)\, J^{\top},$$
where $J$ is the Jacobian of $g$ evaluated at $\widehat{\theta}$.

### Univariate example: rate to mean lifetime

An exponential component has rate $\lambda$. The mean lifetime is
$\mu = 1/\lambda$. Transform the MLE:

``` r
fit_rate <- mle(
  theta.hat = c(lambda = 0.02),
  sigma = matrix(0.001^2),
  nobs = 150L
)

# g: rate -> mean lifetime
fit_lifetime <- rmap(fit_rate, function(theta) 1 / theta[1], method = "delta")
params(fit_lifetime)
#> [1] 50
se(fit_lifetime)
#> [1] 2.5
confint(fit_lifetime)
#>        2.5% 97.5%
#> param1 45.1  54.9
```

### Multivariate example: component rates to system reliability

For a series system with independent exponential components, the system
reliability at time $t$ is
$$R(t) = \exp\left( - \lambda_{1}t \right) \cdot \exp\left( - \lambda_{2}t \right) = \exp\left( - \left( \lambda_{1} + \lambda_{2} \right)t \right).$$

Starting from the joint MLE of $\left( \lambda_{1},\lambda_{2} \right)$:

``` r
fit_rates <- joint(
  mle(theta.hat = c(lambda1 = 0.02), sigma = matrix(0.001^2), nobs = 150L),
  mle(theta.hat = c(lambda2 = 0.05), sigma = matrix(0.003^2), nobs = 80L)
)

# System reliability at t = 100 hours
t0 <- 100
R_mle <- rmap(fit_rates,
  function(theta) exp(-(theta[1] + theta[2]) * t0),
  method = "delta"
)

params(R_mle)
#> [1] 0.0009119
se(R_mle)
#> [1] 0.0002884
confint(R_mle)
#>             2.5%    97.5%
#> param1 0.0003467 0.001477
```

The delta method automatically computes the Jacobian via numerical
differentiation (using
[`numDeriv::jacobian`](https://rdrr.io/pkg/numDeriv/man/jacobian.html))
and propagates the covariance.

## Bridging to Distribution Algebra

[`as_dist()`](https://queelius.github.io/algebraic.dist/reference/as_dist.html)
converts an MLE to its asymptotic normal (or multivariate normal)
distribution, bridging into the `algebraic.dist` package where
distributions can be further composed.

``` r
fit1 <- mle(theta.hat = c(lambda = 0.02), sigma = matrix(0.001^2))
d1 <- as_dist(fit1)
d1
#> Normal distribution (mu = 0.02, var = 1e-06)
```

This is useful when you want to reason about the sampling distribution
of the MLE as a distribution object. For instance, you can compute
probabilities that the true parameter exceeds a threshold:

``` r
# P(lambda > 0.025) under the asymptotic distribution
1 - cdf(d1)(0.025)
#> [1] 2.867e-07
```

For bootstrap MLEs,
[`as_dist()`](https://queelius.github.io/algebraic.dist/reference/as_dist.html)
returns an empirical distribution built from the bootstrap replicates:

``` r
d_boot <- as_dist(fit_boot)
d_boot
#> Empirical distribution (999 observations, 1 dimensions)
```

## Full Pipeline

Here is an end-to-end example tying together all the algebraic
operations.

**Scenario**: A series system has two independent components with
exponential lifetimes. Separate accelerated life tests estimate each
component’s failure rate. We want to infer the system reliability at a
mission time of $t = 500$ hours, including a confidence interval.

``` r
set.seed(123)

# --- Step 1: Independent MLEs from separate experiments ---

# Component 1: 200 observed lifetimes
x1 <- rexp(200, rate = 0.003)
loglik1 <- function(rate) {
  if (rate <= 0) return(-Inf)
  sum(dexp(x1, rate = rate, log = TRUE))
}
fit1 <- mle_numerical(
  optim(par = c(lambda1 = 0.001), fn = loglik1,
        method = "Brent", lower = 1e-6, upper = 1,
        hessian = TRUE, control = list(fnscale = -1)),
  options = list(nobs = length(x1))
)

# Component 2: 120 observed lifetimes
x2 <- rexp(120, rate = 0.008)
loglik2 <- function(rate) {
  if (rate <= 0) return(-Inf)
  sum(dexp(x2, rate = rate, log = TRUE))
}
fit2 <- mle_numerical(
  optim(par = c(lambda2 = 0.001), fn = loglik2,
        method = "Brent", lower = 1e-6, upper = 1,
        hessian = TRUE, control = list(fnscale = -1)),
  options = list(nobs = length(x2))
)

cat("Component 1 rate:", params(fit1), "  SE:", se(fit1), "\n")
#> Component 1 rate: 0.002978   SE: 0.0001827
cat("Component 2 rate:", params(fit2), "  SE:", se(fit2), "\n")
#> Component 2 rate: 0.007983   SE: 0.0007171
```

``` r
# --- Step 2: Compose into joint parameter space ---
fit_system <- joint(fit1, fit2)
cat("Joint parameters:", params(fit_system), "\n")
#> Joint parameters: 0.002978 0.007983
cat("Joint vcov:\n")
#> Joint vcov:
vcov(fit_system)
#>           [,1]      [,2]
#> [1,] 3.336e-08 0.000e+00
#> [2,] 0.000e+00 5.142e-07
```

``` r
# --- Step 3: Transform to system reliability at t = 500 ---
mission_time <- 500
R_system <- rmap(fit_system,
  function(theta) exp(-(theta[1] + theta[2]) * mission_time),
  method = "delta"
)

cat("System reliability R(500):", params(R_system), "\n")
#> System reliability R(500): 0.004167
cat("SE of R(500):             ", se(R_system), "\n")
#> SE of R(500):              0.001542
```

``` r
# --- Step 4: Inference ---
cat("95% CI for R(500):\n")
#> 95% CI for R(500):
confint(R_system)
#>            2.5%    97.5%
#> param1 0.001145 0.007188
```

``` r
# --- Step 5: Bridge to distribution algebra ---
R_dist <- as_dist(R_system)
cat("Asymptotic distribution of R(500):", "\n")
#> Asymptotic distribution of R(500):
R_dist
#> Normal distribution (mu = 0.00416658, var = 2.3765e-06)

# Probability that system reliability exceeds 90%
cat("P(R(500) > 0.90):", 1 - cdf(R_dist)(0.90), "\n")
#> P(R(500) > 0.90): 0
```

The pipeline reads as a chain of algebraic operations: two independent
MLEs are composed via
[`joint()`](https://queelius.github.io/algebraic.mle/reference/joint.md),
transformed to the quantity of interest via
[`rmap()`](https://queelius.github.io/algebraic.dist/reference/rmap.html),
and then reasoned about as a distribution via
[`as_dist()`](https://queelius.github.io/algebraic.dist/reference/as_dist.html).
Each step preserves the uncertainty structure inherited from the
original experiments.
