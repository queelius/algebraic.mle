## ── as_dist.mle (univariate) ─────────────────────────────────────────────

test_that("as_dist.mle univariate returns normal", {
  fit <- mle(theta.hat = c(mu = 5), sigma = matrix(0.25))
  d <- as_dist(fit)

  expect_true(algebraic.dist::is_normal(d))
  expect_equal(mean(d), 5)
  expect_equal(vcov(d), 0.25)
})

test_that("as_dist.mle scalar sigma handled", {
  fit <- mle(theta.hat = c(x = 3), sigma = 0.1)
  d <- as_dist(fit)

  expect_true(algebraic.dist::is_normal(d))
  expect_equal(mean(d), 3)
  expect_equal(vcov(d), 0.1)
})

## ── as_dist.mle (multivariate) ──────────────────────────────────────────

test_that("as_dist.mle multivariate returns mvn", {
  theta <- c(a = 1, b = 2)
  sigma <- matrix(c(0.1, 0.02, 0.02, 0.2), 2, 2)
  fit <- mle(theta.hat = theta, sigma = sigma)
  d <- as_dist(fit)

  expect_true(algebraic.dist::is_mvn(d))
  expect_equal(mean(d), c(1, 2))
  expect_equal(vcov(d), sigma)
})

## ── as_dist.mle errors ──────────────────────────────────────────────────

test_that("as_dist.mle errors without vcov", {
  fit <- mle(theta.hat = c(x = 1))
  expect_error(as_dist(fit), "no variance-covariance")
})

## ── as_dist.mle_boot ────────────────────────────────────────────────────

test_that("as_dist.mle_boot returns empirical_dist", {
  set.seed(42)
  x <- rexp(50, rate = 2)
  rate_mle <- function(data, indices) 1 / mean(data[indices])
  boot_result <- boot::boot(data = x, statistic = rate_mle, R = 100)
  fit <- mle_boot(boot_result)

  d <- as_dist(fit)
  expect_true(algebraic.dist::is_empirical_dist(d))
  expect_equal(nrow(algebraic.dist::obs(d)), 100)
})

test_that("as_dist.mle_boot multivariate returns empirical_dist", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)
  stat_fn <- function(data, indices) {
    d <- data[indices]
    c(mu = mean(d), var = var(d))
  }
  boot_result <- boot::boot(data = x, statistic = stat_fn, R = 200)
  fit <- mle_boot(boot_result)

  d <- as_dist(fit)
  expect_true(algebraic.dist::is_empirical_dist(d))
  expect_equal(dim(d), 2)
  expect_equal(nrow(algebraic.dist::obs(d)), 200)
})

## ── Algebra participation via as_dist ───────────────────────────────────

test_that("as_dist enables distribution algebra (Normal + Normal)", {
  fit1 <- mle(theta.hat = c(x = 3), sigma = matrix(1))
  fit2 <- mle(theta.hat = c(y = 4), sigma = matrix(2))

  d1 <- as_dist(fit1)
  d2 <- as_dist(fit2)
  d_sum <- d1 + d2

  # Normal + Normal = Normal with mean = 3+4, var = 1+2
  expect_true(algebraic.dist::is_normal(d_sum))
  expect_equal(mean(d_sum), 7)
  expect_equal(vcov(d_sum), 3)
})

test_that("as_dist enables scalar multiplication", {
  fit <- mle(theta.hat = c(x = 5), sigma = matrix(4))
  d <- as_dist(fit)

  d_scaled <- 2 * d
  expect_true(algebraic.dist::is_normal(d_scaled))
  expect_equal(mean(d_scaled), 10)
  expect_equal(vcov(d_scaled), 16)
})
