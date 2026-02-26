## ── Delta method ──────────────────────────────────────────────────────────

test_that("rmap delta method transforms point estimate", {
  fit <- mle(
    theta.hat = c(mu = 4, var = 9),
    sigma = diag(c(0.1, 0.5)),
    nobs = 100L
  )

  g <- function(theta) sqrt(theta[2])  # sd = sqrt(var)
  sd_mle <- rmap(fit, g, method = "delta")

  expect_true(is_mle(sd_mle))
  expect_true(inherits(sd_mle, "rmap_mle"))
  expect_equal(as.numeric(params(sd_mle)), 3, tolerance = 1e-6)
})

test_that("rmap delta method propagates variance correctly", {
  # For g(x) = 2*x, Var(g(x)) = 4 * Var(x)
  fit <- mle(theta.hat = c(x = 5), sigma = matrix(0.25), nobs = 50L)
  g <- function(theta) 2 * theta
  g_mle <- rmap(fit, g, method = "delta")

  expect_equal(as.numeric(vcov(g_mle)), 4 * 0.25, tolerance = 1e-4)
})

test_that("rmap delta method preserves nobs", {
  fit <- mle(theta.hat = c(x = 1), sigma = matrix(1), nobs = 42L)
  g_mle <- rmap(fit, function(t) t^2, method = "delta")
  expect_equal(nobs(g_mle), 42L)
})

## ── Monte Carlo method ────────────────────────────────────────────────────

test_that("rmap MC method transforms point estimate", {
  fit <- mle(
    theta.hat = c(mu = 4, var = 9),
    sigma = diag(c(0.1, 0.5)),
    nobs = 100L
  )

  g <- function(theta) sqrt(theta[2])
  set.seed(42)
  sd_mle <- rmap(fit, g, method = "mc", n = 5000)

  expect_true(is_mle(sd_mle))
  expect_equal(as.numeric(params(sd_mle)), 3, tolerance = 1e-6)
})

test_that("rmap MC and delta give similar results", {
  fit <- mle(
    theta.hat = c(x = 10),
    sigma = matrix(0.5),
    nobs = 100L
  )
  g <- function(theta) log(theta)

  delta_mle <- rmap(fit, g, method = "delta")
  set.seed(42)
  mc_mle <- rmap(fit, g, method = "mc", n = 10000)

  # Point estimates should match (both = g(theta))
  expect_equal(params(delta_mle), params(mc_mle))

  # Variances should be similar (not exact due to MC noise)
  expect_equal(as.numeric(vcov(delta_mle)), as.numeric(vcov(mc_mle)),
               tolerance = 0.05)
})

## ── n parameter accepts numeric ───────────────────────────────────────────

test_that("rmap accepts numeric n (not just integer)", {
  fit <- mle(theta.hat = c(x = 5), sigma = matrix(0.1), nobs = 50L)
  # This should NOT error with numeric n
  set.seed(1)
  result <- rmap(fit, function(t) t^2, method = "mc", n = 100)
  expect_true(is_mle(result))
})

## ── Multi-dimensional output ──────────────────────────────────────────────

test_that("rmap handles vector-valued transformation", {
  fit <- mle(theta.hat = c(a = 1, b = 2), sigma = diag(c(0.1, 0.2)))

  g <- function(theta) c(theta[1] + theta[2], theta[1] * theta[2])
  g_mle <- rmap(fit, g, method = "delta")

  expect_equal(length(params(g_mle)), 2)
  expect_equal(as.numeric(params(g_mle)), c(3, 2), tolerance = 1e-6)
  expect_true(is.matrix(vcov(g_mle)))
  expect_equal(dim(vcov(g_mle)), c(2, 2))
})
