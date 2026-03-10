## -- stats4 compatibility regression test -------------------------------------
## The S3 class was renamed from "mle" to "mle_fit" specifically to avoid
## collision with stats4::mle (S4). This test verifies the fix holds.

test_that("all generics work when stats4 is loaded", {
  skip_if_not_installed("stats4")
  requireNamespace("stats4", quietly = TRUE)

  fit <- mle(theta.hat = c(x = 5), sigma = matrix(0.1),
             loglike = -20, nobs = 100L)

  expect_no_error(summary(fit))
  expect_no_error(vcov(fit))
  expect_no_error(confint(fit))
  expect_no_error(nobs(fit))
  expect_no_error(AIC(fit))
  expect_no_error(BIC(fit))
  expect_no_error(logLik(fit))
  expect_no_error(coef(fit))
})
