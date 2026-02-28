## ── density.mle ──────────────────────────────────────────────────────────

test_that("density.mle univariate returns correct PDF", {
  fit <- mle(theta.hat = c(mu = 0), sigma = matrix(1))
  f <- density(fit)

  # PDF of N(0,1) at 0 is 1/sqrt(2*pi) ~ 0.3989
  expect_equal(f(0), dnorm(0), tolerance = 1e-10)
  expect_equal(f(1), dnorm(1), tolerance = 1e-10)
})

test_that("density.mle multivariate returns correct PDF", {
  fit <- mle(theta.hat = c(a = 0, b = 0), sigma = diag(2))
  f <- density(fit)

  # PDF of MVN(0, I) at origin
  expect_equal(f(c(0, 0)), mvtnorm::dmvnorm(c(0, 0)), tolerance = 1e-10)
})

test_that("density.mle errors without vcov", {
  fit <- mle(theta.hat = c(x = 1))
  expect_error(density(fit), "no variance-covariance")
})

## ── cdf.mle ─────────────────────────────────────────────────────────────

test_that("cdf.mle univariate returns correct CDF", {
  fit <- mle(theta.hat = c(mu = 0), sigma = matrix(1))
  F <- cdf(fit)

  expect_equal(F(0), 0.5, tolerance = 1e-10)
  expect_equal(F(qnorm(0.975)), 0.975, tolerance = 1e-6)
})

test_that("cdf.mle multivariate returns correct CDF", {
  fit <- mle(theta.hat = c(a = 0, b = 0), sigma = diag(2))
  F <- cdf(fit)

  # CDF at (Inf, Inf) should be ~1
  val <- F(c(10, 10))
  expect_true(val > 0.999)
})

## ── inv_cdf.mle ─────────────────────────────────────────────────────────

test_that("inv_cdf.mle returns correct quantiles", {
  fit <- mle(theta.hat = c(mu = 5), sigma = matrix(4))
  Q <- inv_cdf(fit)

  expect_equal(Q(0.5), 5, tolerance = 1e-10)
  # Q(0.975) should be 5 + 1.96*2 = 8.92
  expect_equal(Q(0.975), qnorm(0.975, mean = 5, sd = 2), tolerance = 1e-6)
})

## ── sup.mle ─────────────────────────────────────────────────────────────

test_that("sup.mle returns interval for normal", {
  fit <- mle(theta.hat = c(x = 0), sigma = matrix(1))
  s <- sup(fit)

  expect_true(inherits(s, "interval"))
})

## ── dim.mle ─────────────────────────────────────────────────────────────

test_that("dim.mle returns nparams", {
  fit1 <- mle(theta.hat = c(x = 1))
  fit3 <- mle(theta.hat = c(a = 1, b = 2, c = 3))

  expect_equal(dim(fit1), 1)
  expect_equal(dim(fit3), 3)
})

## ── mean.mle ────────────────────────────────────────────────────────────

test_that("mean.mle returns params", {
  fit <- mle(theta.hat = c(a = 1, b = 2))
  expect_equal(mean(fit), params(fit))
})

## ── conditional.mle ─────────────────────────────────────────────────────

test_that("conditional.mle closed-form Schur complement", {
  sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  fit <- mle(theta.hat = c(a = 0, b = 0), sigma = sigma)

  # Condition on a=1, get distribution of b|a=1
  d <- conditional(fit, given_indices = 1, given_values = 1)

  # Schur complement: mu_cond = 0 + 0.5*1*(1-0) = 0.5
  # var_cond = 1 - 0.5*1*0.5 = 0.75
  expect_true(algebraic.dist::is_normal(d))
  expect_equal(mean(d), 0.5, tolerance = 1e-10)
  expect_equal(vcov(d), 0.75, tolerance = 1e-10)
})

test_that("conditional.mle with 3 params", {
  sigma <- diag(3)
  sigma[1, 2] <- sigma[2, 1] <- 0.3
  sigma[1, 3] <- sigma[3, 1] <- 0.1
  sigma[2, 3] <- sigma[3, 2] <- 0.2
  fit <- mle(theta.hat = c(a = 1, b = 2, c = 3), sigma = sigma)

  # Condition on c=5, get distribution of (a,b)|c=5
  d <- conditional(fit, given_indices = 3, given_values = 5)

  expect_true(algebraic.dist::is_mvn(d))
  expect_equal(length(mean(d)), 2)
})

## ── density.mle_boot ────────────────────────────────────────────────────

test_that("density.mle_boot returns empirical PMF", {
  set.seed(42)
  x <- rexp(50, rate = 2)
  rate_mle <- function(data, indices) 1 / mean(data[indices])
  boot_result <- boot::boot(data = x, statistic = rate_mle, R = 100)
  fit <- mle_boot(boot_result)

  f <- density(fit)
  expect_true(is.function(f))
})

## ── dim.mle_boot ────────────────────────────────────────────────────────

test_that("dim.mle_boot returns nparams", {
  set.seed(42)
  x <- rexp(50, rate = 2)
  rate_mle <- function(data, indices) 1 / mean(data[indices])
  boot_result <- boot::boot(data = x, statistic = rate_mle, R = 100)
  fit <- mle_boot(boot_result)

  expect_equal(dim(fit), 1)
})

## ── mean.mle_boot ───────────────────────────────────────────────────────

test_that("mean.mle_boot returns empirical mean of replicates", {
  set.seed(42)
  x <- rexp(50, rate = 2)
  rate_mle <- function(data, indices) 1 / mean(data[indices])
  boot_result <- boot::boot(data = x, statistic = rate_mle, R = 500)
  fit <- mle_boot(boot_result)

  # mean.mle_boot should be mean of replicates (includes bias)
  expect_equal(mean(fit), mean(boot_result$t), tolerance = 1e-10)
})

test_that("mean.mle_boot multivariate returns colMeans", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)
  stat_fn <- function(data, indices) {
    d <- data[indices]
    c(mu = mean(d), var = var(d))
  }
  boot_result <- boot::boot(data = x, statistic = stat_fn, R = 200)
  fit <- mle_boot(boot_result)

  expect_equal(mean(fit), colMeans(boot_result$t), tolerance = 1e-10)
})
