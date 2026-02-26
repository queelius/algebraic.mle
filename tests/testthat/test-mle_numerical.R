## ── Constructor from optim output ──────────────────────────────────────────

test_that("mle_numerical wraps optim result correctly", {
  set.seed(42)
  x <- rexp(200, rate = 2)

  loglik <- function(rate) {
    if (rate <= 0) return(-Inf)
    sum(dexp(x, rate = rate, log = TRUE))
  }

  result <- optim(
    par = 1, fn = loglik, method = "Brent",
    lower = 0.01, upper = 20,
    hessian = TRUE, control = list(fnscale = -1)
  )

  fit <- mle_numerical(result, options = list(nobs = length(x)))

  expect_true(is_mle(fit))
  expect_true(inherits(fit, "mle_numerical"))
  expect_true(fit$converged)

  # MLE of rate should be close to 1/mean(x)
  expect_equal(params(fit), 1 / mean(x), tolerance = 0.01)
  expect_equal(nobs(fit), 200L)
})

test_that("mle_numerical computes vcov from hessian", {
  set.seed(42)
  x <- rexp(200, rate = 2)

  loglik <- function(rate) {
    if (rate <= 0) return(-Inf)
    sum(dexp(x, rate = rate, log = TRUE))
  }

  result <- optim(
    par = 1, fn = loglik, method = "Brent",
    lower = 0.01, upper = 20,
    hessian = TRUE, control = list(fnscale = -1)
  )

  fit <- mle_numerical(result)

  # vcov should be -ginv(hessian)
  expect_false(is.null(vcov(fit)))
  expect_true(all(diag(vcov(fit)) > 0))  # variances positive
})

test_that("mle_numerical stores log-likelihood value", {
  result <- list(
    par = c(mu = 5),
    value = -120.5,
    convergence = 0,
    hessian = matrix(-10)
  )

  fit <- mle_numerical(result)
  expect_equal(loglik_val(fit), -120.5)
})

test_that("mle_numerical stores observed FIM", {
  H <- matrix(-10)
  result <- list(par = c(x = 1), value = -50, convergence = 0, hessian = H)

  fit <- mle_numerical(result)
  # info = -hessian
  expect_equal(observed_fim(fit), -H)
})

## ── Base R generics on mle_numerical ──────────────────────────────────────

test_that("coef and logLik work on mle_numerical", {
  result <- list(
    par = c(a = 1, b = 2),
    value = -80,
    convergence = 0,
    hessian = -diag(c(10, 20))
  )
  fit <- mle_numerical(result, options = list(nobs = 50L))

  expect_equal(coef(fit), c(a = 1, b = 2))

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), -80)
  expect_equal(attr(ll, "df"), 2)

  expect_equal(AIC(fit), -2 * (-80) + 2 * 2)
})

## ── Multivariate fit ──────────────────────────────────────────────────────

test_that("mle_numerical works with multivariate optim", {
  set.seed(42)
  x <- rnorm(100, mean = 3, sd = 2)

  loglik <- function(par) {
    mu <- par[1]; sigma <- par[2]
    if (sigma <= 0) return(-Inf)
    sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
  }

  result <- optim(
    par = c(mu = 0, sigma = 1), fn = loglik,
    hessian = TRUE, control = list(fnscale = -1)
  )

  fit <- mle_numerical(result)

  expect_equal(length(params(fit)), 2)
  expect_equal(names(params(fit)), c("mu", "sigma"))
  expect_true(is.matrix(vcov(fit)))
  expect_equal(dim(vcov(fit)), c(2, 2))
})
