# algebraic.mle v2.0.0 Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Rename S3 class `"mle"` to `"mle_fit"` to eliminate stats4 collision, reclaim stats generics, and clean up code quality.

**Architecture:** Mechanical class rename across all R source, tests, and docs. Remove `mle_weighted` (replaced by `combine()`). Remove custom `aic`/`bic`/`loglik_val` generics; implement proper `logLik.mle_fit` so `AIC()`/`BIC()` work via `stats::AIC.default`. Add `coef.mle_fit` as alias for `params()`.

**Tech Stack:** R, testthat, roxygen2, devtools

---

### Task 1: Delete mle_weighted and fixing/ directory

Reduce surface area before the rename.

**Files:**
- Delete: `R/mle_weighted.R`
- Delete: `tests/testthat/test-mle_weighted.R`
- Delete: `fixing/mle_weighted_slow.R`
- Delete: `fixing/hyp_tests.R`
- Modify: `R/algebraic.mle.R` (remove any mle_weighted re-export if present)

**Step 1:** Delete the files

```bash
rm R/mle_weighted.R tests/testthat/test-mle_weighted.R
rm -rf fixing/
```

**Step 2:** Run `devtools::document()` to regenerate NAMESPACE (removes mle_weighted exports)

```r
devtools::document()
```

**Step 3:** Run tests — expect mle_weighted tests gone, all others pass

```r
devtools::test()
```

**Step 4:** Commit

```
feat: remove mle_weighted (replaced by combine())
```

---

### Task 2: Class rename in R source files

Mechanical find-and-replace. Constructors keep their names; only the class strings and S3 method names change.

**Files:**
- Modify: `R/mle.R`
- Modify: `R/mle_numerical.R`
- Modify: `R/mle_boot.R`
- Modify: `R/mle_rmap.R`
- Modify: `R/combine.R`
- Modify: `R/joint.R`
- Modify: `R/as_dist.R`
- Modify: `R/algebraic.mle.R`

**Rename map (class strings in `class =` / `inherits()` / `structure()`):**

| Old string | New string | Files |
|-----------|-----------|-------|
| `"mle"` (as class) | `"mle_fit"` | mle.R:57, mle_boot.R:39, combine.R:52, joint.R:43 |
| `"mle_numerical"` | `"mle_fit_numerical"` | mle_numerical.R:54 |
| `"mle_boot"` | `"mle_fit_boot"` | mle_boot.R:39, mle_boot.R:53 |
| `"rmap_mle"` | `"mle_fit_rmap"` | mle_rmap.R:73 |
| `"summary_mle"` | `"summary_mle_fit"` | mle.R:293 |

**Rename map (S3 method function names):**

All `*.mle` methods in `R/mle.R` become `*.mle_fit`:
- `print.mle` → `print.mle_fit`
- `vcov.mle` → `vcov.mle_fit`
- `params.mle` → `params.mle_fit`
- `nparams.mle` → `nparams.mle_fit`
- `aic.mle` → (will be removed in Task 4)
- `bic.mle` → (will be removed in Task 4)
- `nobs.mle` → `nobs.mle_fit`
- `obs.mle` → `obs.mle_fit`
- `loglik_val.mle` → (will be removed in Task 4)
- `confint.mle` → `confint.mle_fit`
- `sampler.mle` → `sampler.mle_fit`
- `mse.mle` → `mse.mle_fit`
- `observed_fim.mle` → `observed_fim.mle_fit`
- `summary.mle` → `summary.mle_fit`
- `print.summary_mle` → `print.summary_mle_fit`
- `se.mle` → `se.mle_fit`
- `orthogonal.mle` → `orthogonal.mle_fit`
- `score_val.mle` → `score_val.mle_fit`
- `bias.mle` → `bias.mle_fit`
- `pred.mle` → `pred.mle_fit`
- `expectation.mle` → `expectation.mle_fit`
- `marginal.mle` → `marginal.mle_fit`
- `density.mle` → `density.mle_fit`
- `cdf.mle` → `cdf.mle_fit`
- `inv_cdf.mle` → `inv_cdf.mle_fit`
- `sup.mle` → `sup.mle_fit`
- `dim.mle` → `dim.mle_fit`
- `mean.mle` → `mean.mle_fit`
- `conditional.mle` → `conditional.mle_fit`
- `combine.mle` → `combine.mle_fit`
- `joint.mle` → `joint.mle_fit`
- `as_dist.mle` → `as_dist.mle_fit`
- `rmap.mle` → `rmap.mle_fit`

All `*.mle_boot` methods in `R/mle_boot.R` become `*.mle_fit_boot`:
- `params.mle_boot` → `params.mle_fit_boot`
- `nparams.mle_boot` → `nparams.mle_fit_boot`
- `nobs.mle_boot` → `nobs.mle_fit_boot`
- `obs.mle_boot` → `obs.mle_fit_boot`
- `vcov.mle_boot` → `vcov.mle_fit_boot`
- `mse.mle_boot` → `mse.mle_fit_boot`
- `bias.mle_boot` → `bias.mle_fit_boot`
- `sampler.mle_boot` → `sampler.mle_fit_boot`
- `confint.mle_boot` → `confint.mle_fit_boot`
- `density.mle_boot` → `density.mle_fit_boot`
- `dim.mle_boot` → `dim.mle_fit_boot`
- `mean.mle_boot` → `mean.mle_fit_boot`
- `as_dist.mle_boot` → `as_dist.mle_fit_boot`

**Predicate functions (update inherits target, keep function name):**
- `is_mle(x)`: `inherits(x, "mle")` → `inherits(x, "mle_fit")`
- `is_mle_boot(x)`: `inherits(x, "mle_boot")` → `inherits(x, "mle_fit_boot")`

**Roxygen `\code{}` references:** Update all `\code{mle}` → `\code{mle_fit}` etc. in doc comments.

**Error messages:** Update `"mle objects"` → `"mle_fit objects"` in combine.R, joint.R.

**Step 1:** Apply all renames using editor

**Step 2:** Run `devtools::document()` to regenerate NAMESPACE

**Step 3:** Do NOT run tests yet (Task 3 updates tests first)

**Step 4:** Commit

```
refactor: rename S3 class "mle" to "mle_fit" (stats4 collision fix)
```

---

### Task 3: Update tests for class rename

Mechanical find-and-replace in test files.

**Files:**
- Modify: `tests/testthat/test-mle.R`
- Modify: `tests/testthat/test-mle_numerical.R`
- Modify: `tests/testthat/test-mle_boot.R`
- Modify: `tests/testthat/test-mle_rmap.R`
- Modify: `tests/testthat/test-combine.R`
- Modify: `tests/testthat/test-joint.R`
- Modify: `tests/testthat/test-as_dist.R`
- Modify: `tests/testthat/test-dist_methods.R`

**Renames in tests:**
- `"summary_mle"` → `"summary_mle_fit"` (test-mle.R:210)
- `c("my_custom_mle", "mle")` → `c("my_custom_mle", "mle_fit")` (test-mle.R:222)
- `inherits(fit, "mle_numerical")` → `inherits(fit, "mle_fit_numerical")` (test-mle_numerical.R:21)
- `inherits(fit, "mle_boot")` → `inherits(fit, "mle_fit_boot")` (test-mle_boot.R:17)
- `inherits(sd_mle, "rmap_mle")` → `inherits(sd_mle, "mle_fit_rmap")` (test-mle_rmap.R:14)
- All `is_mle()` / `is_mle_boot()` calls stay as-is (function names unchanged)

**Also in test-mle.R:**
- Remove aic/bic tests (will be replaced in Task 5)
- Keep everything else

**Step 1:** Apply renames

**Step 2:** Run tests — all should pass

```r
devtools::test()
```

**Step 3:** Commit

```
test: update class assertions for mle_fit rename
```

---

### Task 4: Remove old generics, implement logLik and coef

Replace custom generics with standard R generics.

**Files:**
- Modify: `R/algebraic.mle.R` — remove `aic()`, `bic()`, `loglik_val()` generics
- Modify: `R/mle.R` — remove `aic.mle_fit`, `bic.mle_fit`, `loglik_val.mle_fit`; add `logLik.mle_fit`, `coef.mle_fit`
- Modify: `R/mle_boot.R` — check if loglik_val referenced
- Modify: `R/as_dist.R` — check if loglik_val referenced
- Modify: `R/combine.R` — check references
- Modify: `tests/testthat/test-mle.R` — replace aic/bic tests with AIC/BIC/logLik/coef tests

**Step 1: Write failing tests**

Add to `tests/testthat/test-mle.R`:

```r
## -- logLik -------------------------------------------------------------------

test_that("logLik returns proper logLik object", {
  fit <- mle(theta.hat = c(a = 1, b = 2), sigma = diag(2),
             loglike = -50, nobs = 100L)
  ll <- logLik(fit)

  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), -50)
  expect_equal(attr(ll, "df"), 2)
  expect_equal(attr(ll, "nobs"), 100L)
})

test_that("logLik returns NA logLik when loglike is NULL", {
  fit <- mle(theta.hat = c(x = 1))
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.na(as.numeric(ll)))
})

## -- AIC/BIC via logLik -------------------------------------------------------

test_that("AIC computes correctly via logLik", {
  fit <- mle(theta.hat = c(a = 1, b = 2), sigma = diag(2),
             loglike = -50, nobs = 100L)
  # AIC = -2 * loglik + 2 * k
  expect_equal(AIC(fit), -2 * (-50) + 2 * 2)
})

test_that("BIC computes correctly via logLik", {
  fit <- mle(theta.hat = c(a = 1, b = 2), sigma = diag(2),
             loglike = -50, nobs = 100L)
  # BIC = -2 * loglik + k * log(n)
  expect_equal(BIC(fit), -2 * (-50) + 2 * log(100))
})

## -- coef ---------------------------------------------------------------------

test_that("coef delegates to params", {
  fit <- mle(theta.hat = c(mu = 5, var = 4), sigma = diag(2))
  expect_equal(coef(fit), params(fit))
  expect_equal(names(coef(fit)), c("mu", "var"))
})
```

**Step 2: Run tests — expect failures**

```r
devtools::test(filter = "mle$")
```

**Step 3: Remove old generics from `R/algebraic.mle.R`**

Remove `aic`, `bic`, `loglik_val` generic definitions and their roxygen.
Keep `observed_fim`, `mse`, `bias`, `score_val`, `se`, `orthogonal`, `pred`.

**Step 4: Remove old methods from `R/mle.R`**

Remove `aic.mle_fit`, `bic.mle_fit`, `loglik_val.mle_fit`.

Add:

```r
#' Log-likelihood of an \code{mle_fit} object.
#'
#' Returns a \code{"logLik"} object with \code{df} (number of parameters) and
#' \code{nobs} attributes, enabling \code{AIC()} and \code{BIC()} via
#' \code{stats::AIC.default}.
#'
#' @param object the \code{mle_fit} object
#' @param ... additional arguments (not used)
#' @return A \code{"logLik"} object.
#' @importFrom stats logLik
#' @export
logLik.mle_fit <- function(object, ...) {
    val <- object$loglike
    if (is.null(val)) val <- NA_real_
    structure(val,
              df = nparams(object),
              nobs = nobs(object),
              class = "logLik")
}

#' Extract coefficients from an \code{mle_fit} object.
#'
#' Delegates to \code{\link[algebraic.dist]{params}()}.
#'
#' @param object the \code{mle_fit} object
#' @param ... additional arguments (not used)
#' @return Named numeric vector of parameter estimates.
#' @importFrom stats coef
#' @export
coef.mle_fit <- function(object, ...) {
    params(object)
}
```

**Step 5: Update internal references**

In `R/mle.R`, `print.summary_mle_fit` currently calls `aic(x$x)`. Change to `AIC(x$x)`.

Search all R/ files for remaining calls to `aic(`, `bic(`, `loglik_val(` and update:
- `loglik_val(x)` → `object$loglike` (internal) or `logLik(x)` (if going through generic)
- In `print.summary_mle_fit`: `loglik_val(x$x)` → `as.numeric(logLik(x$x))`

**Step 6: Run `devtools::document()`**

**Step 7: Run tests — all should pass**

```r
devtools::test()
```

**Step 8: Commit**

```
feat: replace aic/bic/loglik_val with standard logLik/AIC/BIC/coef
```

---

### Task 5: Add stats4 collision regression test

**Files:**
- Create: `tests/testthat/test-stats4-compat.R`

**Step 1: Write the test**

```r
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
```

**Step 2: Run the test with stats4**

```r
devtools::test(filter = "stats4")
```

**Step 3: Commit**

```
test: add stats4 collision regression test
```

---

### Task 6: Code quality fixes

**Files:**
- Modify: `R/mle.R` — seq_len, vectorize pred.mle_fit
- Modify: `R/mle_boot.R` — seq_len in confint loop
- Modify: `R/mle_rmap.R` — seq_len, coerce n to integer, vectorize MC loop
- Modify: `R/joint.R` — seq_len
- Modify: `R/algebraic.mle.R` — modernize to `"_PACKAGE"` sentinel

**Step 1: seq_len / seq_along replacements**

In `R/mle.R`:
- `for (i in 2:R)` → `for (i in seq_len(R)[-1])` (pred.mle_fit, line ~460)
- `for (j in 1:p)` → `for (j in seq_len(p))` (pred.mle_fit, line ~465)

In `R/mle_boot.R`:
- `for (j in 1:p)` → `for (j in seq_len(p))` (confint.mle_fit_boot, line ~206)

In `R/mle_rmap.R`:
- `for (i in 1:n)` → `for (i in seq_len(n))` (rmap.mle_fit, line ~57)
- Add `n <- as.integer(n)` after the stopifnot validation

**Step 2: Vectorize PI computation in pred.mle_fit**

Replace the `for (j in 1:p)` loop with:

```r
PI <- t(apply(data, 2, function(dj) {
    c(mean(dj), quantile(x = dj, probs = c(alpha / 2, 1 - alpha / 2)))
}))
colnames(PI) <- c("mean", "lower", "upper")
```

**Step 3: Modernize package docs in algebraic.mle.R**

Replace:
```r
#' @docType package
#' @name algebraic.mle
NULL
#> NULL
```

With:
```r
#' @keywords internal
"_PACKAGE"
```

**Step 4: Run tests**

```r
devtools::test()
```

**Step 5: Commit**

```
refactor: code quality — seq_len, vectorize loops, modernize package docs
```

---

### Task 7: Update vignettes and README

**Files:**
- Modify: `vignettes/statistics.Rmd` — `mle_weighted()` → `combine()`
- Modify: `vignettes/fitting-common-dist.Rmd` — `mle_weighted()` → `combine()`
- Modify: `vignettes/mle-algebra.Rmd` — check for aic/bic/loglik_val references
- Modify: `README.Rmd` — update any class name references in prose

**Key changes:**
- All `mle_weighted(list(...))` → `combine(...)`
- Any `aic(fit)` → `AIC(fit)`, `bic(fit)` → `BIC(fit)`
- Any `loglik_val(fit)` → `as.numeric(logLik(fit))` or just `logLik(fit)`
- Prose references to class `"mle"` → `"mle_fit"`

**Step 1:** Apply changes to each vignette

**Step 2:** Rebuild README

```r
rmarkdown::render("README.Rmd")
```

**Step 3:** Commit

```
docs: update vignettes and README for v2.0 class rename
```

---

### Task 8: Update downstream packages

**Files (likelihood.model):**
- Modify: `/home/spinoza/github/rlang/likelihood.model/R/core-fisher_mle.R`
  - Line 86: `class = c("fisher_mle", "mle")` → `c("fisher_mle", "mle_fit")`
  - Line 435: `class = c("fisher_boot", "fisher_mle", "mle", "boot")` → `c("fisher_boot", "fisher_mle", "mle_fit", "boot")`
  - Remove: `aic.fisher_mle`, `bic.fisher_mle`, `bic` generic, `loglik_val.fisher_mle`
  - Update re-export block: drop `aic`, `bic`, `loglik_val` from `@importFrom algebraic.mle`
- Modify: `/home/spinoza/github/rlang/likelihood.model/DESCRIPTION`
  - Line 18: `algebraic.mle` → `algebraic.mle (>= 2.0.0)`

**Files (compositional.mle):**
- Modify: `/home/spinoza/github/rlang/compositional.mle/DESCRIPTION`
  - Line 21: `algebraic.mle` → `algebraic.mle (>= 2.0.0)`

**Step 1:** Apply likelihood.model changes

**Step 2:** Run likelihood.model tests

```r
setwd("/home/spinoza/github/rlang/likelihood.model")
devtools::document()
devtools::test()
```

**Step 3:** Apply compositional.mle version pin

**Step 4:** Commit each package separately

```
# In likelihood.model:
feat: update for algebraic.mle v2.0 class rename

# In compositional.mle:
chore: pin algebraic.mle >= 2.0.0
```

---

### Task 9: Version bump, NEWS.md, final check

**Files:**
- Modify: `DESCRIPTION` — version 1.2.0 → 2.0.0
- Modify: `NEWS.md` — add v2.0.0 entry
- Modify: `CITATION.cff` — update version if present

**Step 1:** Bump version in DESCRIPTION

**Step 2:** Write NEWS.md entry

```markdown
# algebraic.mle 2.0.0

## Breaking changes

* S3 class renamed from `"mle"` to `"mle_fit"` to resolve collision with
  `stats4::mle`. Subclasses follow: `mle_fit_numerical`, `mle_fit_boot`,
  `mle_fit_rmap`. Constructor function names (`mle()`, `mle_numerical()`,
  `mle_boot()`) are unchanged.

* Removed `aic()`, `bic()`, and `loglik_val()` generics. Use `AIC()`,
  `BIC()`, and `logLik()` from the `stats` package instead. `logLik()`
  returns a proper `"logLik"` object with `df` and `nobs` attributes.

* Added `coef.mle_fit()` method delegating to `params()`.

* Removed `mle_weighted()` constructor and `mle_fit_weighted` class.
  Use `combine()` instead, which has the same inverse-variance weighting
  with better error handling.

## Improvements

* Code quality: replaced `1:n` with `seq_len()`, vectorized loops in
  `pred()` and `rmap()`, modernized package documentation.

* Cleaned up `fixing/` directory (removed dead experimental code).
```

**Step 3:** Run `devtools::document()`

**Step 4:** Run full R CMD check

```r
devtools::check()
```

Expected: 0 errors, 0 warnings, 0 notes (or only the "New submission" note).

**Step 5:** Verify stats4 compatibility with installed package

```r
devtools::install()
# In a fresh R session:
library(stats4)
library(algebraic.mle)
fit <- mle(theta.hat = c(x = 5), sigma = matrix(0.1), loglike = -20, nobs = 100L)
summary(fit)
AIC(fit)
BIC(fit)
logLik(fit)
coef(fit)
vcov(fit)
confint(fit)
nobs(fit)
```

**Step 6:** Commit

```
release: algebraic.mle v2.0.0
```
