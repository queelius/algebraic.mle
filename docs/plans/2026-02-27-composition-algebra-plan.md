# MLE Composition Algebra Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add `joint()` and `combine()` composition operators to `algebraic.mle`, rewrite documentation around "MLE as technology" narrative, bump to v1.2.0.

**Architecture:** Two new generics (`joint`, `combine`) with `mle` methods in separate files. `joint()` builds block-diagonal covariance from independent MLEs with disjoint parameters. `combine()` does inverse-variance weighting of same-parameter MLEs (cleaner API over existing `mle_weighted`). Documentation reframed around the idea that the MLE is a technology you exploit once you have an estimator.

**Tech Stack:** R, S3 generics, testthat 3, roxygen2, knitr/rmarkdown for vignettes.

**Design doc:** `docs/plans/2026-02-27-composition-algebra-design.md`

---

### Task 1: `joint()` — Tests

**Files:**
- Create: `tests/testthat/test-joint.R`

**Step 1: Write failing tests**

```r
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
```

**Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "joint")'`
Expected: FAIL — `joint` function does not exist

**Step 3: Commit failing tests**

```bash
git add tests/testthat/test-joint.R
git commit -m "test: add failing tests for joint() composition operator"
```

---

### Task 2: `joint()` — Implementation

**Files:**
- Create: `R/joint.R`
- Modify: `R/algebraic.mle.R` (add export for `joint`)

**Step 1: Write `R/joint.R`**

```r
#' Compose independent MLEs into a joint MLE.
#'
#' Given two or more independent MLEs with disjoint parameter sets, produces
#' a joint MLE with block-diagonal variance-covariance structure.
#'
#' The joint MLE has:
#' \itemize{
#'   \item \code{theta.hat}: concatenation of all parameter vectors
#'   \item \code{sigma}: block-diagonal from individual vcov matrices
#'   \item \code{loglike}: sum of log-likelihoods (when all available)
#'   \item \code{info}: block-diagonal from individual FIMs (when all available)
#'   \item \code{score}: concatenation of score vectors (when all available)
#'   \item \code{nobs}: NULL (different experiments have no shared sample size)
#' }
#'
#' @param x An \code{mle} object.
#' @param ... Additional \code{mle} objects to join.
#' @return An \code{mle} object representing the joint MLE.
#' @examples
#' # Two independent experiments
#' fit_rate <- mle(theta.hat = c(lambda = 2.1), sigma = matrix(0.04), nobs = 50L)
#' fit_shape <- mle(theta.hat = c(k = 1.5, s = 3.2),
#'                  sigma = matrix(c(0.1, 0.02, 0.02, 0.3), 2, 2), nobs = 100L)
#'
#' # Joint MLE: 3 params, block-diagonal covariance
#' j <- joint(fit_rate, fit_shape)
#' params(j)   # c(lambda = 2.1, k = 1.5, s = 3.2)
#' vcov(j)     # 3x3 block-diagonal
#'
#' # Existing algebra works on the joint:
#' marginal(j, 2:3)   # recover shape params
#' as_dist(j)          # MVN for distribution algebra
#' @export
joint <- function(x, ...) {
    UseMethod("joint", x)
}

#' @rdname joint
#' @importFrom algebraic.dist params nparams
#' @importFrom stats vcov nobs
#' @importFrom MASS ginv
#' @export
joint.mle <- function(x, ...) {
    mles <- c(list(x), list(...))

    if (length(mles) < 2L) {
        stop("joint() requires at least 2 mle objects.")
    }
    if (!all(vapply(mles, is_mle, logical(1)))) {
        stop("All arguments to joint() must be mle objects.")
    }

    # Validate disjoint parameter names
    all_names <- unlist(lapply(mles, function(m) names(params(m))))
    if (anyDuplicated(all_names)) {
        dups <- all_names[duplicated(all_names)]
        stop("Parameter names must be disjoint. Overlapping: ",
             paste(unique(dups), collapse = ", "))
    }

    # All must have vcov
    vcovs <- lapply(mles, vcov)
    if (any(vapply(vcovs, is.null, logical(1)))) {
        stop("All mle objects must have a variance-covariance matrix for joint().")
    }

    # Concatenate parameters
    theta_joint <- unlist(lapply(mles, params))

    # Build block-diagonal vcov
    dims <- vapply(mles, nparams, integer(1))
    p <- sum(dims)
    sigma_joint <- matrix(0, p, p)
    offset <- 0L
    for (i in seq_along(mles)) {
        idx <- offset + seq_len(dims[i])
        V <- vcovs[[i]]
        if (!is.matrix(V)) V <- matrix(V, 1, 1)
        sigma_joint[idx, idx] <- V
        offset <- offset + dims[i]
    }

    # Sum log-likelihoods (NULL if any missing)
    loglikes <- lapply(mles, loglik_val)
    loglike_joint <- if (any(vapply(loglikes, is.null, logical(1)))) {
        NULL
    } else {
        sum(unlist(loglikes))
    }

    # Block-diagonal FIM (NULL if any missing)
    fims <- lapply(mles, observed_fim)
    info_joint <- if (any(vapply(fims, is.null, logical(1)))) {
        NULL
    } else {
        info <- matrix(0, p, p)
        offset <- 0L
        for (i in seq_along(mles)) {
            idx <- offset + seq_len(dims[i])
            I <- fims[[i]]
            if (!is.matrix(I)) I <- matrix(I, 1, 1)
            info[idx, idx] <- I
            offset <- offset + dims[i]
        }
        info
    }

    # Concatenate scores (NULL if any missing)
    scores <- lapply(mles, score_val)
    score_joint <- if (any(vapply(scores, is.null, logical(1)))) {
        NULL
    } else {
        unlist(scores)
    }

    mle(theta.hat = theta_joint,
        loglike = loglike_joint,
        score = score_joint,
        sigma = sigma_joint,
        info = info_joint,
        obs = NULL,
        nobs = NULL)
}
```

**Step 2: Add export to `R/algebraic.mle.R`**

After the existing generics (around line 177), add nothing — `joint` defines its own generic and is exported via roxygen `@export`. No re-export needed since `joint` is native to this package.

**Step 3: Regenerate NAMESPACE**

Run: `Rscript -e 'devtools::document()'`

**Step 4: Run tests**

Run: `Rscript -e 'devtools::test(filter = "joint")'`
Expected: All PASS

**Step 5: Run full test suite for regressions**

Run: `Rscript -e 'devtools::test()'`
Expected: All 182+ tests PASS

**Step 6: Commit**

```bash
git add R/joint.R NAMESPACE man/
git commit -m "feat: add joint() for block-diagonal composition of independent MLEs"
```

---

### Task 3: `combine()` — Tests

**Files:**
- Create: `tests/testthat/test-combine.R`

**Step 1: Write failing tests**

```r
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
```

**Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "combine")'`
Expected: FAIL — `combine` function does not exist

**Step 3: Commit failing tests**

```bash
git add tests/testthat/test-combine.R
git commit -m "test: add failing tests for combine() same-parameter MLE composition"
```

---

### Task 4: `combine()` — Implementation

**Files:**
- Create: `R/combine.R`

**Step 1: Write `R/combine.R`**

```r
#' Combine independent MLEs for the same parameter.
#'
#' Given multiple independent MLEs that estimate the same parameter \eqn{\theta},
#' produces an optimally weighted combination using inverse-variance (Fisher
#' information) weighting.
#'
#' The combined estimator has:
#' \itemize{
#'   \item \code{theta.hat}: \eqn{(\sum I_i)^{-1} \sum I_i \hat\theta_i}
#'   \item \code{sigma}: \eqn{(\sum I_i)^{-1}}
#'   \item \code{info}: \eqn{\sum I_i}
#'   \item \code{nobs}: sum of individual sample sizes
#' }
#'
#' When the Fisher information matrix is not directly available but the
#' variance-covariance matrix is, the FIM is computed as \code{ginv(vcov)}.
#'
#' For the legacy interface that accepts a list, see \code{\link{mle_weighted}}.
#'
#' @param x An \code{mle} object, or a list of \code{mle} objects.
#' @param ... Additional \code{mle} objects to combine.
#' @return An \code{mle} object representing the optimally weighted combination.
#' @seealso \code{\link{mle_weighted}}, \code{\link{joint}}
#' @examples
#' # Three independent estimates of the same rate
#' fit1 <- mle(theta.hat = c(lambda = 2.1), sigma = matrix(0.04), nobs = 50L)
#' fit2 <- mle(theta.hat = c(lambda = 1.9), sigma = matrix(0.02), nobs = 100L)
#' fit3 <- mle(theta.hat = c(lambda = 2.0), sigma = matrix(0.03), nobs = 70L)
#'
#' comb <- combine(fit1, fit2, fit3)
#' params(comb)
#' se(comb)
#' @export
combine <- function(x, ...) {
    UseMethod("combine", x)
}

#' @rdname combine
#' @export
combine.list <- function(x, ...) {
    if (length(list(...)) > 0L) {
        stop("When passing a list, do not pass additional arguments.")
    }
    combine_mles(x)
}

#' @rdname combine
#' @importFrom MASS ginv
#' @importFrom algebraic.dist params nparams
#' @importFrom stats vcov nobs
#' @export
combine.mle <- function(x, ...) {
    dots <- list(...)
    if (length(dots) == 0L) return(x)
    combine_mles(c(list(x), dots))
}

#' Internal workhorse for combine.
#' @param mles A list of mle objects.
#' @return An mle object.
#' @keywords internal
combine_mles <- function(mles) {
    if (!all(vapply(mles, is_mle, logical(1)))) {
        stop("All arguments to combine() must be mle objects.")
    }
    if (length(mles) == 1L) return(mles[[1L]])

    # Get FIMs, falling back to ginv(vcov) when needed
    fims <- lapply(mles, function(m) {
        I <- observed_fim(m)
        if (!is.null(I)) return(if (is.matrix(I)) I else matrix(I, 1, 1))
        V <- vcov(m)
        if (is.null(V)) {
            stop("combine() requires either a Fisher information matrix or ",
                 "variance-covariance matrix for each mle object.")
        }
        if (!is.matrix(V)) V <- matrix(V, 1, 1)
        ginv(V)
    })

    info_combined <- Reduce(`+`, fims)
    sigma_combined <- ginv(info_combined)
    theta_combined <- as.vector(
        sigma_combined %*%
        Reduce(`+`, Map(`%*%`, fims, lapply(mles, function(m) {
            p <- params(m)
            matrix(p, ncol = 1)
        })))
    )
    names(theta_combined) <- names(params(mles[[1L]]))

    # Sum nobs (use 0L for NULL, then convert back)
    nobs_vals <- vapply(mles, function(m) {
        n <- nobs(m)
        if (is.null(n)) NA_integer_ else as.integer(n)
    }, integer(1))
    nobs_combined <- if (anyNA(nobs_vals)) NULL else sum(nobs_vals)

    mle(theta.hat = theta_combined,
        loglike = NULL,
        score = NULL,
        sigma = sigma_combined,
        info = info_combined,
        obs = NULL,
        nobs = nobs_combined)
}
```

**Step 2: Regenerate NAMESPACE**

Run: `Rscript -e 'devtools::document()'`

**Step 3: Run combine tests**

Run: `Rscript -e 'devtools::test(filter = "combine")'`
Expected: All PASS

**Step 4: Run full test suite**

Run: `Rscript -e 'devtools::test()'`
Expected: All tests PASS (no regressions)

**Step 5: Commit**

```bash
git add R/combine.R NAMESPACE man/
git commit -m "feat: add combine() for optimal weighting of same-parameter MLEs"
```

---

### Task 5: Version Bump + DESCRIPTION Rewrite

**Files:**
- Modify: `DESCRIPTION`

**Step 1: Bump version and rewrite Description**

Change `Version:` from `1.1.0` to `1.2.0`.

Rewrite the `Description:` field to:

```
Description: The maximum likelihood estimator (MLE) is a technology: under
    regularity conditions, any MLE is asymptotically normal with variance given
    by the inverse Fisher information. This package exploits that structure by
    defining an algebra over MLEs. Compose independent estimators into joint
    MLEs via block-diagonal covariance ('joint'), optimally combine repeated
    estimates via inverse-variance weighting ('combine'), propagate
    transformations via the delta method ('rmap'), and bridge to distribution
    algebra via conversion to normal or multivariate normal objects ('as_dist').
    Supports asymptotic ('mle', 'mle_numerical') and bootstrap ('mle_boot')
    estimators with a unified interface for inference: confidence intervals,
    standard errors, AIC, Fisher information, and predictive intervals.
    For background on maximum likelihood estimation, see Casella and Berger
    (2002, ISBN:978-0534243128). For the delta method and variance estimation,
    see Lehmann and Casella (1998, ISBN:978-0387985022).
```

**Step 2: Verify**

Run: `Rscript -e 'devtools::check(args = "--no-vignettes")' 2>&1 | tail -5`
Expected: 0 errors, 0 warnings

**Step 3: Commit**

```bash
git add DESCRIPTION
git commit -m "docs: bump to v1.2.0, rewrite Description with MLE-as-technology framing"
```

---

### Task 6: README.Rmd Rewrite

**Files:**
- Modify: `README.Rmd`

**Step 1: Rewrite README.Rmd**

Replace the current content with a restructured version organized around:

1. **Opening paragraph** — "MLE as technology" framing (the MLE is asymptotically normal, this package is the algebra that exploits that)
2. **Installation** — same as current
3. **The Algebra** — quick showcase of the four operations:
   - `joint()` — compose independent MLEs (new)
   - `combine()` — optimally weight same-parameter MLEs (new)
   - `rmap()` — transform via invariance / delta method (existing)
   - `as_dist()` — bridge to distribution algebra (new from earlier work)
4. **Inference** — quick showcase of `confint`, `se`, `aic`, `sampler` (existing, keep brief)
5. **Creating MLE Objects** — `mle()`, `mle_numerical()`, `mle_boot()` (reuse some existing content)
6. **Links** — vignettes, pkgdown site

Keep the conditional exponential fitting example but move it to section 5 as a motivating example for `mle_numerical()`. The top of the README should lead with the algebra, not model fitting.

**Step 2: Knit to verify**

Run: `Rscript -e 'rmarkdown::render("README.Rmd")'`
Expected: Renders without error

**Step 3: Commit**

```bash
git add README.Rmd README.md
git commit -m "docs: rewrite README around MLE-as-technology narrative"
```

---

### Task 7: New Vignette — "The Algebra of MLEs"

**Files:**
- Create: `vignettes/mle-algebra.Rmd`

**Step 1: Write vignette**

Structure:
1. **Introduction** — The MLE is a technology. Under regularity conditions, theta_hat ~ MVN(theta, I(theta)^{-1}). This package is the algebra that follows from this fact.
2. **Creating MLEs** — Brief: `mle()` for known quantities, `mle_numerical()` from optimization, `mle_boot()` from bootstrap.
3. **Composing Independent MLEs** — `joint()` example: two independent experiments estimating different parameters of a system. Show block-diagonal structure. Then use `marginal()` to slice back.
4. **Combining Repeated Estimates** — `combine()` example: three labs estimate the same rate parameter. Show inverse-variance weighting, variance reduction.
5. **Transformations via Invariance** — `rmap()` example: transform rate to mean lifetime. Show delta method variance propagation.
6. **Bridging to Distribution Algebra** — `as_dist()` example: convert to normal/MVN, use algebraic.dist operators (addition, scalar multiplication, conditioning).
7. **Full Pipeline** — Component MLEs → `joint()` → `rmap()` to system reliability → `as_dist()` for distribution algebra → `confint()` for inference.

**Step 2: Build vignette to verify**

Run: `Rscript -e 'rmarkdown::render("vignettes/mle-algebra.Rmd")'`
Expected: Renders without error

**Step 3: Commit**

```bash
git add vignettes/mle-algebra.Rmd
git commit -m "docs: add 'The Algebra of MLEs' vignette"
```

---

### Task 8: Final Verification

**Step 1: Full test suite**

Run: `Rscript -e 'devtools::test()'`
Expected: All tests PASS (should be ~210+)

**Step 2: Coverage check**

Run: `Rscript -e 'cov <- covr::package_coverage(); cat("Coverage:", covr::percent_coverage(cov), "%\n")'`
Expected: >= 80%

**Step 3: R CMD check**

Run: `Rscript -e 'devtools::check()'`
Expected: 0 errors, 0 warnings

**Step 4: Verify backwards compatibility**

Run: `Rscript -e 'library(algebraic.mle); fit <- mle_weighted(list(mle(theta.hat=1, sigma=matrix(1), info=matrix(1), nobs=10L), mle(theta.hat=2, sigma=matrix(1), info=matrix(1), nobs=10L))); cat("mle_weighted still works:", params(fit), "\n")'`
Expected: Prints combined parameter value (1.5)

**Step 5: Commit any remaining changes**

---

### Task Summary

| Task | Component | Type |
|------|-----------|------|
| 1 | joint() tests | TDD: failing tests |
| 2 | joint() implementation | TDD: make tests pass |
| 3 | combine() tests | TDD: failing tests |
| 4 | combine() implementation | TDD: make tests pass |
| 5 | DESCRIPTION rewrite + version bump | Documentation |
| 6 | README.Rmd rewrite | Documentation |
| 7 | mle-algebra.Rmd vignette | Documentation |
| 8 | Final verification | Verification |
