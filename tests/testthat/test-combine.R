## -- combine: basic same-parameter combination --------------------------------

test_that("combine two univariate MLEs gives optimal weighting", {
  # Exact inverse-variance weighting
  fit1 <- mle(theta.hat = c(mu = 10), sigma = matrix(0.04),
              info = matrix(25), nobs = 50L)
  fit2 <- mle(theta.hat = c(mu = 11), sigma = matrix(0.02),
              info = matrix(50), nobs = 100L)

  comb <- combine(fit1, fit2)

  expect_true(is_mle(comb))
  # Combined info = 25 + 50 = 75
  expect_equal(observed_fim(comb), matrix(75))
  # Combined theta = (75)^{-1} * (25*10 + 50*11) = 800/75
  expect_equal(as.numeric(params(comb)), 800 / 75, tolerance = 1e-10)
  # nobs = 150
  expect_equal(nobs(comb), 150L)
})

test_that("combine multivariate MLEs", {
  sigma1 <- matrix(c(0.1, 0.01, 0.01, 0.2), 2, 2)
  sigma2 <- matrix(c(0.05, 0.005, 0.005, 0.1), 2, 2)
  info1 <- solve(sigma1)
  info2 <- solve(sigma2)

  fit1 <- mle(theta.hat = c(a = 1, b = 2), sigma = sigma1,
              info = info1, nobs = 50L)
  fit2 <- mle(theta.hat = c(a = 1.1, b = 2.1), sigma = sigma2,
              info = info2, nobs = 80L)

  comb <- combine(fit1, fit2)

  expect_equal(nparams(comb), 2)
  expect_equal(nobs(comb), 130L)
  # Combined variance should be less than either individual
  expect_true(all(diag(vcov(comb)) < diag(sigma1)))
  expect_true(all(diag(vcov(comb)) < diag(sigma2)))
})

## -- combine: varargs and list ------------------------------------------------

test_that("combine works with 3+ varargs", {
  fit1 <- mle(theta.hat = c(x = 10), sigma = matrix(1), info = matrix(1), nobs = 10L)
  fit2 <- mle(theta.hat = c(x = 11), sigma = matrix(1), info = matrix(1), nobs = 10L)
  fit3 <- mle(theta.hat = c(x = 12), sigma = matrix(1), info = matrix(1), nobs = 10L)

  comb <- combine(fit1, fit2, fit3)
  # Equal weights: theta = (10+11+12)/3 = 11
  expect_equal(as.numeric(params(comb)), 11, tolerance = 1e-10)
  expect_equal(nobs(comb), 30L)
})

test_that("combine accepts a list", {
  mles <- lapply(1:4, function(i) {
    mle(theta.hat = c(x = i), sigma = matrix(1), info = matrix(1), nobs = 10L)
  })

  comb <- combine(mles)
  expect_equal(as.numeric(params(comb)), 2.5, tolerance = 1e-10)
})

## -- combine: FIM fallback from vcov ------------------------------------------

test_that("combine falls back to ginv(vcov) when FIM unavailable", {
  fit1 <- mle(theta.hat = c(x = 10), sigma = matrix(0.04), nobs = 50L)
  fit2 <- mle(theta.hat = c(x = 11), sigma = matrix(0.02), nobs = 100L)

  # Neither has info field, but both have sigma
  expect_null(observed_fim(fit1))

  comb <- combine(fit1, fit2)
  expect_true(!is.null(vcov(comb)))
  expect_true(!is.null(observed_fim(comb)))
})

## -- combine: passthrough and validation --------------------------------------

test_that("combine of single MLE returns it unchanged", {
  fit <- mle(theta.hat = c(x = 5), sigma = matrix(0.1), nobs = 50L)
  expect_identical(combine(fit), fit)
})

test_that("combine errors on non-mle input", {
  fit <- mle(theta.hat = c(x = 1), sigma = matrix(0.1))
  expect_error(combine(fit, 42), "mle")
})

test_that("combine errors when no vcov or FIM available", {
  fit1 <- mle(theta.hat = c(x = 1))
  fit2 <- mle(theta.hat = c(x = 2))
  expect_error(combine(fit1, fit2), "variance-covariance|information")
})

## -- combine: matches mle_weighted -------------------------------------------

test_that("combine matches mle_weighted for FIM-available inputs", {
  set.seed(123)
  make_mle <- function(mu, n) {
    s2 <- 4
    mle(theta.hat = mu, sigma = matrix(s2 / n),
        info = matrix(n / s2), nobs = n)
  }
  fit1 <- make_mle(10, 50L)
  fit2 <- make_mle(11, 30L)
  fit3 <- make_mle(9.5, 70L)

  comb <- combine(fit1, fit2, fit3)
  weighted <- mle_weighted(list(fit1, fit2, fit3))

  expect_equal(as.numeric(params(comb)), as.numeric(params(weighted)),
               tolerance = 1e-10)
  expect_equal(as.numeric(vcov(comb)), as.numeric(vcov(weighted)),
               tolerance = 1e-10)
})
