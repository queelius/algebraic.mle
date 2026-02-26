## ── Basic weighting ───────────────────────────────────────────────────────

test_that("mle_weighted combines two estimates", {
  fit1 <- mle(theta.hat = c(mu = 10), sigma = matrix(1),
              info = matrix(1), nobs = 50L)
  fit2 <- mle(theta.hat = c(mu = 12), sigma = matrix(0.5),
              info = matrix(2), nobs = 100L)

  combined <- mle_weighted(list(fit1, fit2))

  expect_true(is_mle(combined))
  expect_true(inherits(combined, "mle_weighted"))

  # Weighted mean: (1*10 + 2*12) / (1+2) = 34/3
  expect_equal(as.numeric(params(combined)), 34 / 3, tolerance = 1e-10)

  # Combined info = 1 + 2 = 3, vcov = 1/3
  expect_equal(as.numeric(vcov(combined)), 1 / 3, tolerance = 1e-10)

  # Combined nobs = 50 + 100 = 150
  expect_equal(nobs(combined), 150L)
})

## ── Single element ────────────────────────────────────────────────────────

test_that("mle_weighted with single-element list returns it", {
  fit <- mle(theta.hat = c(x = 5), sigma = matrix(1), info = matrix(1))
  combined <- mle_weighted(list(fit))
  expect_equal(params(combined), params(fit))
})

test_that("mle_weighted with single mle returns it directly", {
  fit <- mle(theta.hat = c(x = 5), sigma = matrix(1))
  combined <- mle_weighted(fit)
  expect_equal(params(combined), params(fit))
})

## ── Multivariate ──────────────────────────────────────────────────────────

test_that("mle_weighted works with multivariate parameters", {
  info1 <- diag(c(2, 4))
  info2 <- diag(c(3, 6))
  fit1 <- mle(theta.hat = c(a = 1, b = 2), sigma = solve(info1),
              info = info1, nobs = 30L)
  fit2 <- mle(theta.hat = c(a = 3, b = 4), sigma = solve(info2),
              info = info2, nobs = 70L)

  combined <- mle_weighted(list(fit1, fit2))

  expect_equal(length(params(combined)), 2)
  expect_true(is.matrix(vcov(combined)))
  expect_equal(nobs(combined), 100L)

  # Combined info = info1 + info2
  expected_info <- info1 + info2
  expect_equal(observed_fim(combined), expected_info, tolerance = 1e-10)
})

## ── Error handling ────────────────────────────────────────────────────────

test_that("mle_weighted errors on NULL input", {
  expect_error(mle_weighted(NULL), "null")
})

test_that("mle_weighted errors on non-mle elements", {
  expect_error(mle_weighted(list(1, 2, 3)), "not all elements")
})

test_that("mle_weighted errors on non-list input", {
  expect_error(mle_weighted("foo"), "not a list")
})
