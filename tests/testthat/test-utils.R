test_that("confint_from_sigma works with matrix input", {
  theta <- c(mu = 5, sigma2 = 4)
  sigma <- diag(c(0.25, 0.5))
  ci <- confint_from_sigma(sigma, theta, level = 0.95)

  expect_equal(nrow(ci), 2)
  expect_equal(ncol(ci), 2)
  expect_equal(rownames(ci), c("mu", "sigma2"))

  # CI should be centered on theta
  midpoints <- (ci[, 1] + ci[, 2]) / 2
  expect_equal(midpoints, theta, tolerance = 1e-10)

  # Lower < theta < upper

  expect_true(all(ci[, 1] < theta))
  expect_true(all(ci[, 2] > theta))
})

test_that("confint_from_sigma works with vector input", {
  theta <- c(a = 1, b = 2)
  sigma_vec <- c(0.1, 0.2)
  ci <- confint_from_sigma(sigma_vec, theta, level = 0.95)

  expect_equal(nrow(ci), 2)
  expect_equal(rownames(ci), c("a", "b"))
})

test_that("confint_from_sigma respects level parameter", {
  theta <- c(x = 0)
  sigma <- matrix(1)

  ci_95 <- confint_from_sigma(sigma, theta, level = 0.95)
  ci_99 <- confint_from_sigma(sigma, theta, level = 0.99)

  width_95 <- ci_95[1, 2] - ci_95[1, 1]
  width_99 <- ci_99[1, 2] - ci_99[1, 1]
  expect_true(width_99 > width_95)
})

test_that("confint_from_sigma uses unnamed params correctly", {
  theta <- c(1, 2, 3)
  sigma <- diag(c(0.1, 0.2, 0.3))
  ci <- confint_from_sigma(sigma, theta)

  expect_equal(rownames(ci), c("param1", "param2", "param3"))
})
