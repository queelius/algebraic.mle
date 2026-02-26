## ── Constructor and accessors ──────────────────────────────────────────────

test_that("mle constructor stores all fields", {
  theta <- c(mu = 5, var = 4)
  sigma <- diag(c(0.04, 0.32))
  fit <- mle(
    theta.hat = theta,
    loglike = -150.3,
    score = c(0.001, -0.002),
    sigma = sigma,
    info = solve(sigma),
    obs = 1:100,
    nobs = 100L
  )

  expect_true(is_mle(fit))
  expect_equal(params(fit), theta)
  expect_equal(vcov(fit), sigma)
  expect_equal(loglik_val(fit), -150.3)
  expect_equal(score_val(fit), c(0.001, -0.002))
  expect_equal(observed_fim(fit), solve(sigma))
  expect_equal(obs(fit), 1:100)
  expect_equal(nobs(fit), 100L)
})

test_that("mle constructor handles NULL fields", {
  fit <- mle(theta.hat = c(x = 1))

  expect_true(is_mle(fit))
  expect_equal(params(fit), c(x = 1))
  expect_null(vcov(fit))
  expect_null(loglik_val(fit))
  expect_null(score_val(fit))
  expect_null(observed_fim(fit))
  expect_null(obs(fit))
  expect_null(nobs(fit))
})

test_that("is_mle returns FALSE for non-mle objects", {
  expect_false(is_mle(list(a = 1)))
  expect_false(is_mle(42))
  expect_false(is_mle("hello"))
})

## ── coef and logLik (base R generics) ─────────────────────────────────────

test_that("coef.mle returns same as params", {
  fit <- mle(theta.hat = c(a = 1, b = 2), sigma = diag(2))
  expect_equal(coef(fit), params(fit))
})

test_that("logLik.mle returns proper logLik object", {
  fit <- mle(theta.hat = c(a = 1, b = 2), sigma = diag(2),
             loglike = -50, nobs = 100L)
  ll <- logLik(fit)

  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), -50)
  expect_equal(attr(ll, "df"), 2)
  expect_equal(attr(ll, "nobs"), 100L)
})

test_that("logLik.mle returns NULL when loglike is NULL", {
  fit <- mle(theta.hat = c(x = 1))
  expect_null(logLik(fit))
})

test_that("AIC and BIC work via logLik", {
  fit <- mle(theta.hat = c(a = 1, b = 2), sigma = diag(2),
             loglike = -50, nobs = 100L)
  # AIC = -2*loglik + 2*k
  expect_equal(AIC(fit), -2 * (-50) + 2 * 2)
  # BIC = -2*loglik + k*log(n)
  expect_equal(BIC(fit), -2 * (-50) + 2 * log(100))
})

## ── nparams ───────────────────────────────────────────────────────────────

test_that("nparams returns correct count", {
  fit1 <- mle(theta.hat = c(x = 1))
  fit2 <- mle(theta.hat = c(a = 1, b = 2, c = 3))

  expect_equal(nparams(fit1), 1)
  expect_equal(nparams(fit2), 3)
})

## ── aic ───────────────────────────────────────────────────────────────────

test_that("aic computes correctly", {
  fit <- mle(theta.hat = c(a = 1, b = 2), loglike = -50, sigma = diag(2))
  # AIC = -2 * loglike + 2 * k
  expect_equal(aic(fit), -2 * (-50) + 2 * 2)
})

## ── se ────────────────────────────────────────────────────────────────────

test_that("se returns sqrt of diagonal of vcov", {
  sigma <- matrix(c(4, 1, 1, 9), nrow = 2)
  fit <- mle(theta.hat = c(a = 0, b = 0), sigma = sigma)

  expect_equal(se(fit), c(2, 3))
})

test_that("se returns NULL when vcov is NULL", {
  fit <- mle(theta.hat = c(x = 1))
  expect_null(se(fit))
})

test_that("se with scalar vcov", {
  fit <- mle(theta.hat = c(x = 1), sigma = 4)
  expect_equal(se(fit), 2)
})

## ── confint ───────────────────────────────────────────────────────────────

test_that("confint returns correct structure", {
  fit <- mle(theta.hat = c(mu = 5), sigma = matrix(0.25), nobs = 100L)
  ci <- confint(fit)

  expect_equal(nrow(ci), 1)
  expect_equal(ncol(ci), 2)
  expect_equal(rownames(ci), "mu")
  expect_true(ci[1, 1] < 5)
  expect_true(ci[1, 2] > 5)
})

test_that("confint respects level parameter", {
  fit <- mle(theta.hat = c(x = 0), sigma = matrix(1), nobs = 50L)
  ci_90 <- confint(fit, level = 0.90)
  ci_99 <- confint(fit, level = 0.99)

  expect_true((ci_99[1, 2] - ci_99[1, 1]) > (ci_90[1, 2] - ci_90[1, 1]))
})

test_that("confint errors without vcov", {
  fit <- mle(theta.hat = c(x = 1))
  expect_error(confint(fit), "No variance-covariance matrix")
})

## ── bias ──────────────────────────────────────────────────────────────────

test_that("bias.mle returns zero vector (asymptotic)", {
  fit <- mle(theta.hat = c(a = 1, b = 2), sigma = diag(2))
  expect_equal(bias(fit), c(0, 0))
})

## ── mse ───────────────────────────────────────────────────────────────────

test_that("mse.mle equals vcov when bias is zero", {
  sigma <- diag(c(0.1, 0.2))
  fit <- mle(theta.hat = c(a = 1, b = 2), sigma = sigma)
  expect_equal(mse(fit), sigma)
})

test_that("mse.mle with scalar parameter", {
  fit <- mle(theta.hat = c(x = 5), sigma = 0.25)
  expect_equal(mse(fit), 0.25)
})

## ── orthogonal ────────────────────────────────────────────────────────────

test_that("orthogonal detects diagonal FIM", {
  info <- diag(c(10, 20))
  fit <- mle(theta.hat = c(a = 1, b = 2), sigma = solve(info), info = info)
  orth <- orthogonal(fit)

  # Off-diagonals should be TRUE (near zero), diagonals FALSE (not near zero)
  expect_true(orth[1, 2])
  expect_true(orth[2, 1])
  expect_false(orth[1, 1])
  expect_false(orth[2, 2])
})

test_that("orthogonal returns NULL when no FIM", {
  fit <- mle(theta.hat = c(x = 1))
  expect_null(orthogonal(fit))
})

## ── marginal ──────────────────────────────────────────────────────────────

test_that("marginal extracts correct subset", {
  theta <- c(a = 1, b = 2, c = 3)
  sigma <- diag(c(0.1, 0.2, 0.3))
  fit <- mle(theta.hat = theta, sigma = sigma, nobs = 50L)

  m <- marginal(fit, c(1, 3))
  expect_equal(params(m), c(a = 1, c = 3))
  expect_equal(vcov(m), sigma[c(1, 3), c(1, 3)])
  expect_equal(nobs(m), 50L)
})

test_that("marginal errors on empty or invalid indices", {
  fit <- mle(theta.hat = c(a = 1, b = 2), sigma = diag(2))
  expect_error(marginal(fit, integer(0)), "non-empty")
  expect_error(marginal(fit, c(-1)), "\\[1, dim")
  expect_error(marginal(fit, c(3)), "\\[1, dim")
})

## ── sampler ───────────────────────────────────────────────────────────────

test_that("sampler.mle univariate returns numeric vector", {
  fit <- mle(theta.hat = c(x = 5), sigma = 0.1)
  samp <- sampler(fit)
  set.seed(1)
  draws <- samp(100)

  expect_length(draws, 100)
  expect_true(is.numeric(draws))
})

test_that("sampler.mle multivariate returns matrix", {
  fit <- mle(theta.hat = c(a = 1, b = 2), sigma = diag(c(0.1, 0.2)))
  samp <- sampler(fit)
  set.seed(1)
  draws <- samp(50)

  expect_equal(nrow(draws), 50)
  expect_equal(ncol(draws), 2)
})

## ── summary and print ────────────────────────────────────────────────────

test_that("summary returns summary_mle object", {
  fit <- mle(theta.hat = c(x = 5), sigma = matrix(0.1), loglike = -20)
  s <- summary(fit)
  expect_s3_class(s, "summary_mle")
})

test_that("print does not error", {
  fit <- mle(theta.hat = c(x = 5), sigma = matrix(0.1), loglike = -20)
  expect_output(print(fit), "Maximum likelihood estimator")
})

## ── superclasses ──────────────────────────────────────────────────────────

test_that("superclasses are prepended to class vector", {
  fit <- mle(theta.hat = c(x = 1), superclasses = "my_custom_mle")
  expect_equal(class(fit), c("my_custom_mle", "mle"))
  expect_true(is_mle(fit))
})
