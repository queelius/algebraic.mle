# algebraic.mle v2.0.0 Breaking Change Design

Date: 2026-03-08

## Problem

The S3 class `"mle"` collides with `stats4::mle` (S4 class). When stats4 is
loaded, S4 dispatch intercepts calls to `AIC()`, `BIC()`, `logLik()`, `coef()`,
`summary()`, `vcov()`, `confint()`, and `nobs()` — all fail with
`Error in object@min: no applicable method for '@'`.

stats4 ships with R. Any user with both packages loaded hits this.
CRAN rejected v1.2.0 for this reason.

## Solution: Class Rename + Stats Generics Reclaim

### Class Rename

| Old | New |
|-----|-----|
| `mle` | `mle_fit` |
| `mle_numerical` | `mle_fit_numerical` |
| `mle_boot` | `mle_fit_boot` |
| `rmap_mle` | `mle_fit_rmap` |
| `summary_mle` | `summary_mle_fit` |
| `mle_weighted` | **removed** (replaced by `combine()`) |

Constructor function names stay the same: `mle()`, `mle_numerical()`, etc.
Predicates `is_mle()`, `is_mle_boot()` updated internally.

### Reclaim Stats Generics

Implement `logLik.mle_fit` returning a proper `"logLik"` object with `df` and
`nobs` attributes. `AIC()` and `BIC()` work for free via `stats::AIC.default`.

Add `coef.mle_fit` delegating to `params()`.

Remove custom generics: `aic()`, `bic()`, `loglik_val()`.

`params()` remains the primary accessor (ecosystem generic from algebraic.dist).

### Parameter Naming

Keep current convention: `x` for our own generics, `object` where forced by
stats generics (`vcov`, `nobs`, `confint`, `summary`, `logLik`, `coef`).
Already consistent across ecosystem.

### Code Quality

- Replace all `1:n` / `1:p` with `seq_len()` / `seq_along()`
- Vectorize loops in `pred()`, `mle_rmap.R`, `mle_boot.R`
- Modernize package docs: `"_PACKAGE"` sentinel
- Delete `fixing/mle_weighted_slow.R` and `fixing/hyp_tests.R`
- Delete `R/mle_weighted.R`
- Coerce `n` to integer in `rmap()`

### Tests

- Rename all class references in existing tests
- Add `logLik`/`AIC`/`BIC`/`coef` tests
- Add stats4 collision regression test
- Remove `mle_weighted` and `aic`/`bic` tests

### Downstream: likelihood.model

- Class vectors: `"mle"` -> `"mle_fit"` (2 lines in fisher_mle.R)
- Remove: `aic.fisher_mle`, `bic.fisher_mle`, `bic` generic, `loglik_val.fisher_mle`
  (all superseded by `logLik.fisher_mle` + stats defaults)
- Update re-exports: drop `aic`, `bic`, `loglik_val` from `@importFrom algebraic.mle`
- DESCRIPTION: `algebraic.mle (>= 2.0.0)`

### Downstream: compositional.mle

- 0 code changes (uses `mle_numerical()` constructor, class handled internally)
- DESCRIPTION: `algebraic.mle (>= 2.0.0)`

### Release Order

1. algebraic.mle 2.0.0
2. likelihood.model (next version, pins >= 2.0.0)
3. compositional.mle (next version, pins >= 2.0.0)

### Version

1.2.0 -> 2.0.0
