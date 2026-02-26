## ── Constructor ───────────────────────────────────────────────────────────

test_that("mle_boot wraps boot object correctly", {
  set.seed(42)
  x <- rexp(50, rate = 2)

  rate_mle <- function(data, indices) {
    d <- data[indices]
    1 / mean(d)
  }

  b <- boot::boot(data = x, statistic = rate_mle, R = 200)
  fit <- mle_boot(b)

  expect_true(is_mle(fit))
  expect_true(is_mle_boot(fit))
  expect_true(inherits(fit, "mle_boot"))
})

test_that("is_mle_boot returns FALSE for non-boot objects", {
  fit <- mle(theta.hat = 5, sigma = matrix(0.1))
  expect_false(is_mle_boot(fit))
})

## ── Accessors ─────────────────────────────────────────────────────────────

test_that("params.mle_boot returns original estimates", {
  set.seed(42)
  x <- rexp(50, rate = 2)

  rate_mle <- function(data, indices) {
    d <- data[indices]
    1 / mean(d)
  }

  b <- boot::boot(data = x, statistic = rate_mle, R = 200)
  fit <- mle_boot(b)

  expect_equal(params(fit), b$t0)
})

test_that("nobs.mle_boot returns data length for vector", {
  set.seed(42)
  x <- rexp(30, rate = 1)

  stat <- function(data, i) mean(data[i])
  b <- boot::boot(data = x, statistic = stat, R = 100)
  fit <- mle_boot(b)

  expect_equal(nobs(fit), 30)
})

test_that("obs.mle_boot returns original data", {
  set.seed(42)
  x <- rexp(20, rate = 1)

  stat <- function(data, i) mean(data[i])
  b <- boot::boot(data = x, statistic = stat, R = 100)
  fit <- mle_boot(b)

  expect_equal(obs(fit), x)
})

## ── vcov and se ───────────────────────────────────────────────────────────

test_that("vcov.mle_boot returns covariance of replicates", {
  set.seed(42)
  x <- rexp(50, rate = 2)

  rate_mle <- function(data, indices) {
    d <- data[indices]
    1 / mean(d)
  }

  b <- boot::boot(data = x, statistic = rate_mle, R = 500)
  fit <- mle_boot(b)

  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_true(all(diag(V) > 0))
})

## ── bias ──────────────────────────────────────────────────────────────────

test_that("bias.mle_boot estimates bias from replicates", {
  set.seed(42)
  x <- rexp(50, rate = 2)

  rate_mle <- function(data, indices) {
    d <- data[indices]
    1 / mean(d)
  }

  b <- boot::boot(data = x, statistic = rate_mle, R = 500)
  fit <- mle_boot(b)

  b_est <- bias(fit)
  expect_length(b_est, 1)
  # Bias should be relatively small compared to the estimate
  expect_true(abs(b_est) < abs(params(fit)))
})

test_that("bias.mle_boot warns when theta is provided", {
  set.seed(42)
  x <- rexp(20, rate = 1)

  stat <- function(data, i) mean(data[i])
  b <- boot::boot(data = x, statistic = stat, R = 100)
  fit <- mle_boot(b)

  expect_warning(bias(fit, theta = 1), "ignore")
})

## ── confint ───────────────────────────────────────────────────────────────

test_that("confint.mle_boot returns bootstrap CIs", {
  set.seed(42)
  x <- rexp(50, rate = 2)

  rate_mle <- function(data, indices) {
    d <- data[indices]
    1 / mean(d)
  }

  b <- boot::boot(data = x, statistic = rate_mle, R = 500)
  fit <- mle_boot(b)

  ci <- confint(fit, type = "perc")
  expect_equal(nrow(ci), 1)
  expect_equal(ncol(ci), 2)
  expect_true(ci[1, 1] < params(fit))
  expect_true(ci[1, 2] > params(fit))
})

## ── sampler ───────────────────────────────────────────────────────────────

test_that("sampler.mle_boot resamples from replicates", {
  set.seed(42)
  x <- rexp(50, rate = 2)

  rate_mle <- function(data, indices) {
    d <- data[indices]
    1 / mean(d)
  }

  b <- boot::boot(data = x, statistic = rate_mle, R = 200)
  fit <- mle_boot(b)

  samp <- sampler(fit)
  draws <- samp(100)
  expect_length(draws, 100)
  # All draws should be values from the bootstrap replicates
  expect_true(all(draws %in% b$t))
})

## ── mse ───────────────────────────────────────────────────────────────────

test_that("mse.mle_boot incorporates bias", {
  set.seed(42)
  x <- rexp(50, rate = 2)

  rate_mle <- function(data, indices) {
    d <- data[indices]
    1 / mean(d)
  }

  b <- boot::boot(data = x, statistic = rate_mle, R = 500)
  fit <- mle_boot(b)

  m <- mse(fit)
  expect_true(is.matrix(m))
  # MSE >= variance (MSE = var + bias^2)
  expect_true(all(m >= vcov(fit) - 1e-10))
})

## ── coef dispatches through inheritance ───────────────────────────────────

test_that("coef works on mle_boot via inheritance", {
  set.seed(42)
  x <- rexp(30, rate = 1)

  stat <- function(data, i) mean(data[i])
  b <- boot::boot(data = x, statistic = stat, R = 100)
  fit <- mle_boot(b)

  expect_equal(coef(fit), params(fit))
})
