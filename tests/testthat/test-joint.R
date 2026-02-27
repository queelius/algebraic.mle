## -- joint: basic two-MLE composition ----------------------------------------

test_that("joint of two univariate MLEs produces correct structure", {
  fit1 <- mle(theta.hat = c(lambda = 2.1), sigma = matrix(0.04),
              loglike = -30, info = matrix(25), nobs = 50L,
              score = c(lambda = 0.001))
  fit2 <- mle(theta.hat = c(k = 1.5), sigma = matrix(0.1),
              loglike = -45, info = matrix(10), nobs = 100L,
              score = c(k = -0.002))

  j <- joint(fit1, fit2)

  expect_true(is_mle(j))
  expect_equal(params(j), c(lambda = 2.1, k = 1.5))
  expect_equal(nparams(j), 2)

  # Block-diagonal vcov
  expected_sigma <- matrix(0, 2, 2)
  expected_sigma[1, 1] <- 0.04
  expected_sigma[2, 2] <- 0.1
  expect_equal(vcov(j), expected_sigma)

  # Additive log-likelihood
  expect_equal(loglik_val(j), -75)

  # Block-diagonal FIM
  expected_info <- matrix(0, 2, 2)
  expected_info[1, 1] <- 25
  expected_info[2, 2] <- 10
  expect_equal(observed_fim(j), expected_info)

  # Concatenated score
  expect_equal(score_val(j), c(lambda = 0.001, k = -0.002))

  # nobs is NULL (different experiments)
  expect_null(nobs(j))
})

test_that("joint of univariate and multivariate MLE", {
  fit1 <- mle(theta.hat = c(a = 1), sigma = matrix(0.5))
  fit2 <- mle(theta.hat = c(b = 2, c = 3),
              sigma = matrix(c(0.1, 0.02, 0.02, 0.3), 2, 2))

  j <- joint(fit1, fit2)

  expect_equal(params(j), c(a = 1, b = 2, c = 3))
  expect_equal(nparams(j), 3)
  expect_equal(vcov(j)[1, 1], 0.5)
  expect_equal(vcov(j)[2:3, 2:3], matrix(c(0.1, 0.02, 0.02, 0.3), 2, 2))
  expect_equal(vcov(j)[1, 2:3], c(0, 0))
  expect_equal(vcov(j)[2:3, 1], c(0, 0))
})

test_that("joint of three MLEs (variadic)", {
  fit1 <- mle(theta.hat = c(a = 1), sigma = matrix(0.1))
  fit2 <- mle(theta.hat = c(b = 2), sigma = matrix(0.2))
  fit3 <- mle(theta.hat = c(c = 3), sigma = matrix(0.3))

  j <- joint(fit1, fit2, fit3)

  expect_equal(params(j), c(a = 1, b = 2, c = 3))
  expect_equal(diag(vcov(j)), c(0.1, 0.2, 0.3))
  # Off-diagonals are zero
  V <- vcov(j)
  diag(V) <- 0
  expect_true(all(V == 0))
})

## -- joint: NULL field handling -----------------------------------------------

test_that("joint with one missing loglike sets loglike to NULL", {
  fit1 <- mle(theta.hat = c(a = 1), sigma = matrix(0.1), loglike = -10)
  fit2 <- mle(theta.hat = c(b = 2), sigma = matrix(0.2))

  j <- joint(fit1, fit2)
  expect_null(loglik_val(j))
})

test_that("joint with one missing vcov errors", {
  fit1 <- mle(theta.hat = c(a = 1), sigma = matrix(0.1))
  fit2 <- mle(theta.hat = c(b = 2))

  expect_error(joint(fit1, fit2), "variance-covariance")
})

## -- joint: validation --------------------------------------------------------

test_that("joint errors on overlapping parameter names", {
  fit1 <- mle(theta.hat = c(a = 1, b = 2), sigma = diag(2))
  fit2 <- mle(theta.hat = c(b = 3, c = 4), sigma = diag(2))

  expect_error(joint(fit1, fit2), "disjoint|overlap")
})

test_that("joint errors on non-mle input", {
  fit1 <- mle(theta.hat = c(a = 1), sigma = matrix(0.1))
  expect_error(joint(fit1, "not_an_mle"), "mle")
})

test_that("joint errors with fewer than 2 inputs", {
  fit1 <- mle(theta.hat = c(a = 1), sigma = matrix(0.1))
  expect_error(joint(fit1), "at least 2")
})

## -- joint: algebraic closure -------------------------------------------------

test_that("marginal recovers component from joint", {
  fit1 <- mle(theta.hat = c(a = 1), sigma = matrix(0.5))
  fit2 <- mle(theta.hat = c(b = 2, c = 3),
              sigma = matrix(c(0.1, 0.02, 0.02, 0.3), 2, 2))

  j <- joint(fit1, fit2)
  m <- marginal(j, 2:3)

  expect_equal(params(m), c(b = 2, c = 3))
  expect_equal(vcov(m), matrix(c(0.1, 0.02, 0.02, 0.3), 2, 2))
})

test_that("as_dist works on joint MLE", {
  fit1 <- mle(theta.hat = c(a = 1), sigma = matrix(0.5))
  fit2 <- mle(theta.hat = c(b = 2), sigma = matrix(0.3))

  j <- joint(fit1, fit2)
  d <- as_dist(j)

  expect_true(algebraic.dist::is_mvn(d))
  expect_equal(mean(d), c(1, 2))
})

test_that("rmap works on joint MLE (delta method)", {
  fit1 <- mle(theta.hat = c(a = 2), sigma = matrix(0.1))
  fit2 <- mle(theta.hat = c(b = 3), sigma = matrix(0.2))

  j <- joint(fit1, fit2)
  # Transform: product a*b
  g <- function(p) c(product = p[1] * p[2])
  transformed <- rmap(j, g, method = "delta")

  expect_equal(params(transformed), c(product = 6), tolerance = 1e-6)
  expect_true(!is.null(vcov(transformed)))
})
